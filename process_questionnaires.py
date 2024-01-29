### Script to process Exposure Assessment input files from ATSDR
#   Output: a file where each row is a respondent, including serum concentrations, house dust concentrations, questionnaire data

import os
import pandas as pd
import numpy as np
import math
import scipy.stats as stats

data_dir = os.path.abspath('data_read_only')
output = pd.DataFrame()

print('Beginning process_questionnaires.py...')
# Loop through all sites and append the questionnaire data for all respondents
for site in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']:
    # Load basic respondent info from questionnaires
    try:
        adult = pd.read_excel(os.path.join(data_dir, 'Site{}_AdultQxNonPII_toEPA.xlsx'.format(site)))#, usecols='A:J')
        child = pd.read_excel(os.path.join(data_dir, 'Site{}_ChildQxNonPII_toEPA.xlsx'.format(site)))#, usecols='A:J')
    except:
        adult = pd.read_excel(os.path.join(data_dir, 'Site{}_AdultQx_toEPA.xlsx'.format(site)))#, usecols='A:J')
        child = pd.read_excel(os.path.join(data_dir, 'Site{}_ChildQx_toEPA.xlsx'.format(site)))#, usecols='A:J')
    #correct column name discrepencies if present
    adult = adult.rename(columns={'Age': 'Aage', 'AAge': 'Aage', 'ARespondentID': 'respondentid'})
    child = child.rename(columns={'CRespondentID': 'respondentid'})
    # append to master dataframe
    if site == 'A':
        qdata_all_sites = adult
        qdata_all_sites_child = child
    else:
        qdata_all_sites = pd.concat([qdata_all_sites,adult])
        qdata_all_sites_child = pd.concat([qdata_all_sites_child,child])
qdata_all_sites = qdata_all_sites.set_index('respondentid')
qdata_all_sites_child = qdata_all_sites_child.set_index('respondentid')

qs_adults_path = 'output/questionnaire_data_all_sites_adults.csv'
qs_children_path = 'output/questionnaire_data_all_sites_children.csv'
qdata_all_sites.to_csv(qs_adults_path)
qdata_all_sites_child.to_csv(qs_children_path)
print('Adult questionnaires for all sites written to: {}'.format(qs_adults_path))
print('Child questionnaires for all sites written to: {}'.format(qs_children_path))

# Load dust+serum data and join questionnaire data to each participant
# Append full questionnaire data to our dust+serum dataset
df_all_sites = pd.read_csv('output/data_all_sites.csv')
df_complete = df_all_sites.merge(qdata_all_sites.iloc[:,9:],on='respondentid', how='left')
df_complete = df_complete.merge(qdata_all_sites_child.iloc[:,9:],on='respondentid', how='left', suffixes=('_adult', '_child'))
df_complete = df_complete.set_index('respondentid')

#standardize binary columns
print('Standardizing column formats between sites...')
binary_cols = [col for col in df_complete.columns if df_complete[col].astype(str).isin(['TRUE', 'FALSE', 'True', 'False', 'true', 'false', 'No', 'Yes', '0', '1', '0.0', '1.0', 'nan', 'NaN',
                                                                                       "Don't know", 'Refused to answer']).all()]
for col in binary_cols:
    #print(col)
    df_complete[col] = df_complete[col].replace({"Don't know": np.nan, 'Refused to answer': np.nan, '0': False, '1': True, '0.0': False, '1.0': True, 0.0: False,
                                                 1.0: True, 0: False, 1: True, 'Yes': True, 'yes': True, 'No': False, 'no': False})

# combine duplicate adult and child columns
print('combining duplicate adult and child columns...')
combine_adult_child_cols = {'CQ3_AAorAN': 'AQ2_AAorAN', 'CQ3_Asian': 'AQ2_Asian', 'CQ3_Black': 'AQ2_Black', 'CQ3_Hawaiian': 'AQ2_Hawaiian', 'CQ3_White':'AQ2_White', 'CQ4_Months':'AQ3_Months', 'CQ4_Years':'AQ3_Years', 'CQ5_FullTimeRes':'AQ4_FullTimeRes', 'CQ5_Days':'AQ4_Days', 'CQ5_Weeks':'AQ4_Weeks', 'CQ5_Months':'AQ4_Months', 'CQ6_WaterHomeCups':'AQ12_WaterHomeCups', 'CQ6_NoResponse':'AQ12_NoResponse', 'CQ7_SoilFreq':'AQ19_SoilFreq', 'CQ9_FruitVeg':'AQ21_FruitVegFreq', 'CQ10_Fish':'AQ22_FishFreq', 'CQ11_Milk':'AQ23_MilkFreq', 'CQ17_Blood':'AQ30_Blood', 'PublicWaterSupply_child':'PublicWaterSupply_adult'}
df_complete_1 = df_complete.copy()
for c_col, a_col in combine_adult_child_cols.items():
    df_complete[a_col] = df_complete[a_col].combine_first(df_complete[c_col])
    df_complete.drop(c_col, axis=1, inplace=True)

# writing outputs
out_path = 'output/data_all_sites_questionnaire.csv'
df_complete.to_csv(out_path)
print('Dust + serum + full questionnaire data for all sites written to: {}'.format(out_path))
print('Completed process_questionnaires.py.\n')

