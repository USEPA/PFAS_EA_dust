### Script to process Exposure Assessment input files from ATSDR
#   Output: a file where each row is a respondent, including serum concentrations, house dust concentrations, sex, age and weight from questionnaire

import os
import pandas as pd
import numpy as np
import math
import scipy.stats as stats

data_dir = os.path.abspath('data_read_only')
output = pd.DataFrame()
PFAS = ['PFOA', 'PFOS', 'PFNA', 'PFDA', 'PFUnA', 'PFHxS', 'MeFOSAA']
PFAS_serum_convert = {'Sb-PFOA':'PFOA', 'n-PFOA':'PFOA', 'Sm-PFOS':'PFOS', 'n-PFOS':'PFOS', 'PFNA':'PFNA',
                      'PFDEA':'PFDA', 'PFUA':'PFUnA', 'PFHXS':'PFHxS', 'ME-PFOSA-ACOH':'MeFOSAA'}

# Loop through sites A-E
for site in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']:
    # Load basic respondent info from questionnaires
    try:
        adult = pd.read_excel(os.path.join(data_dir, 'Site{}_AdultQxNonPII_toEPA.xlsx'.format(site)), usecols='A:J')
        child = pd.read_excel(os.path.join(data_dir, 'Site{}_ChildQxNonPII_toEPA.xlsx'.format(site)), usecols='A:J')
    except:
        adult = pd.read_excel(os.path.join(data_dir, 'Site{}_AdultQx_toEPA.xlsx'.format(site)), usecols='A:J')
        child = pd.read_excel(os.path.join(data_dir, 'Site{}_ChildQx_toEPA.xlsx'.format(site)), usecols='A:J')
    if site == 'D':
        adult = adult.rename(columns = {'Age':'AAge'})  # fix a typo in original data column labels
    adult.columns = adult.columns.map(lambda x: x[1:].lower())
    adult['qtype'] = 'Adult'
    child.columns = child.columns.map(lambda x: x[1:].lower())
    child['qtype'] = 'Child'
    subjects = adult.append(child)
    subjects = subjects.set_index('respondentid')
    print('Site {} number of subjects: {}'.format(site,len(subjects)))

    # load and process serum data
    serum = pd.read_excel(os.path.join(data_dir, 'Site{}_serum_results_toEPA.xlsx'.format(site)))
    serum.columns = serum.columns.str.lower()
    serum['sampleid'] = serum['sampleid'].str.replace(' ','') # remove any spaces in IDs
    serum['respondentid'] = serum.sampleid.str.slice(4,8)
    serum['serum_ng/ml'] = serum['level'].copy()
    serum.loc[serum['serum_ng/ml'].isin(['<LOD','ND']), 'serum_ng/ml'] = (0.1 / math.sqrt(2)) # Set <LOD values to LOD/sqrt(2). All LODs for serum are 0.1 ng/ml

    #standardize PFAS abbreviations for serum
    serum = serum.replace({"abbreviation": PFAS_serum_convert})
    serum['serum_col'] = serum['abbreviation'] + '_serum_ng/ml'

    #add together PFOS and PFOA isomers to calculate total PFOS and PFOA in serum
    serum = serum[['respondentid', 'serum_ng/ml', 'serum_col']].copy()
    serum['serum_ng/ml'] = serum['serum_ng/ml'].astype('float') 
    serum_summed = serum.groupby(['respondentid','serum_col'],as_index=False).sum()
    
    #convert serum table to wide format
    serum_summed = serum_summed.pivot(index='respondentid', columns='serum_col', values='serum_ng/ml')

    #join the subjects with their serum data
    df = subjects.join(serum_summed, how='inner') # drops subjects without serum data
    print('Site {} number of subjects with serum data: {}'.format(site,len(df)))

    # read in and clean dust data
    dust= pd.read_excel(os.path.join(data_dir, 'Site{}_dust_as received_toEPA.xlsx'.format(site)),sheet_name=1, header=None)
    dust = dust.loc[dust[0].isin(PFAS+['CLIENT_ID','UNITS'])]  # drop PFAS not included in blood serum
    for i in range(1,len(dust.columns)-1,3):  # fix sample labels
        dust.iloc[0,(i+1)] = dust.iloc[0,i]
        dust.iloc[0,(i+2)] = dust.iloc[0,i]
    dust = dust.transpose()
    dust.columns = dust.iloc[0,:]
    dust = dust.drop(dust.index[0])  # transpose dataframe
    dust = dust.loc[~dust['CLIENT_ID'].str.contains('DEQ|Blank|Spiked')]  # drop blanks
    dust['CLIENT_ID'] = dust['CLIENT_ID'].str.replace(' ','')  # clean up labels
    dust['UNITS'] = dust['UNITS'].str.replace(' \(as_received weight basis\)','')  # clean up labels
    dust['UNITS'] = dust['UNITS'].str.replace(' \(as received weight basis\)','')  # clean up labels
    dust = dust.set_index('CLIENT_ID')
    dust_levels = dust[dust.UNITS == 'ng/g']
    dust_RLs = dust[dust.UNITS == 'ng/g (RL)'].drop(columns='UNITS')
    dust_levels = dust_levels.fillna(dust_RLs/math.sqrt(2))  # replace non-detects with reporting limit / sqrt(2)
    dust_levels['env_id'] = dust_levels.index.str.slice(4,8)
    dust_levels = dust_levels.set_index('env_id').drop(columns=['UNITS'])
    dust_levels.columns += '_dust_ng/g'
    dust_levels= dust_levels.astype('float')

    # link dustIDs to household IDs
    env_to_house= pd.read_excel(os.path.join(data_dir, 'DustID & HouseholdID.xlsx'),sheet_name='Site {}'.format(site), skiprows=2).set_index('Sample ID')
    env_to_house['env_id_trimmed'] = env_to_house.index.str.slice(4,8)
    env_to_house.columns = ['householdid', 'env_id']
    dust_levels = dust_levels.merge(env_to_house, on='env_id')
    dust_levels = dust_levels.groupby(['householdid', 'env_id']).agg(stats.gmean).reset_index('env_id')  # if n>1 per house, take geometric mean
    
    # join subjects with their household dust levels
    df = df.merge(dust_levels, left_on='householdid', right_index=True, how='inner')  # drops all subjects with no dust data
    df['site'] = site
    print('Site {} number of subjects with serum+dust data: {}\n'.format(site,len(df)))

    # append to master dataframe
    if site == 'A':
        df_all_sites = df
        df_all_sites_MRLs = dust_RLs
    else:
        df_all_sites = pd.concat([df_all_sites,df])
        df_all_sites_MRLs = pd.concat([df_all_sites_MRLs, dust_RLs])

print('Min and max MRLs per analyte: \n{}\n'.format(df_all_sites_MRLs.agg([min,'median',max])))   
print('Total number of subjects with serum+dust data: {}\n '.format(len(df_all_sites)))

output_path = 'output/data_all_sites.csv'
df_all_sites.to_csv(output_path)

print('Processed data output to: {}'.format(output_path))
print('data_processing.py script completed.')