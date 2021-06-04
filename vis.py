import os
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.preprocessing import StandardScaler


if __name__ == '__main__':
	log_score = True

	# Load Fst data
	with open('data/preprocessed/Fst_bl10000.json') as json_file:
		Fst_dict = json.load(json_file)

	# Load study scores
	save_dir = 'data/results'

	height_studies = [
		'GCST005908',   #EUR1
		'GCST000817',   #EUR2
		'GCST000611',   #JAP1
		'GCST008839',   #JAP2
		'GCST001263',   #AFR1
		'GCST001290',	#AFR2
		'GCST002702',   #EAS1
		'GCST004212',   #GBR1
		'GCST008163',   #EUR3
		'GCST002647',	#EUR4
		'GCST001885',   #CHB1
		'GCST001956',   #EUR5
		'GCST000372',	#EUR6
		'GCST000174',	#EUR7
		'GCST000175',	#EUR8
		'GCST000176',	#EUR9
	]

	schizophrenia_studies = [
		'GCST006803',   #EUR1_S
		'GCST003880',   #CHB1_S
		'GCST001242',	#EUR2_S
		'GCST007205',	#JPT1_S
		'GCST001757',	#EUR3_S
		'GCST002149',	#EUR4_S
		'GCST001657',	#EUR5_S
		'GCST009336',   #EAS1_S
		'GCST001565',	#EUR6_S
	]

	bmi_studies = [
		'GCST009004',	#EUR1_B
		'GCST009003',	#EUR1B_B
		'GCST009001',	#EUR1c_B
		'GCST002461',	#EAS1_B
		'GCST006369',	#EAS2_B
		'GCST002021',	#EUR2_B
		'GCST005951',	#EUR3_B
		'GCST008158',	#EUR4_B
		'GCST006802',	#GBR1_B
	]

	t2diabetes_studies = [
		'GCST008114',	#AFR1_T2D
		'GCST007847',	#JPT1_T2D
		'GCST010118',	#EAS1_T2D
		'GCST009379',	#EUR1_T2D
		'GCST006801',	#GBR1_T2D
		'GCST004773',	#EUR2_T2D
		'GCST005047',	#EUR3_T2D
		'GCST002560',	#ASW1_T2D
		'GCST002128',	#EAS2_T2D
		'GCST001351',	#EAS3_T2D
		'GCST90013693',	#JPT2_T2D
		'GCST007517',	#EUR4_T2D
		'GCST010556',   #ASW2_T2D
		'GCST000712',	#EUR5_T2D
	]

	panc_cancer_studies = [
		'GCST005786',	#AFR1_PC
		'GCST006085',    #EUR1_PC
	]

	
	beta_score_studies = height_studies + bmi_studies
	or_score_studies = schizophrenia_studies + t2diabetes_studies + panc_cancer_studies

	studies = schizophrenia_studies

	study_to_sample_pop = {
		'GCST005908': 'EUR',    #EUR1
		'GCST000817': 'EUR',    #EUR2
		'GCST000611': 'JPT',    #JAP1
		'GCST008839': 'JPT',    #JAP2
		'GCST001263': 'AFR',    #AFR1
		'GCST001290': 'AFR',	#AFR2
		'GCST002702': 'EAS',    #EAS1
		'GCST004212': 'GBR',	#GBR1
		'GCST008163': 'EUR',    #EUR3
		'GCST002647': 'EUR',    #EUR4
		'GCST001885': 'CHB',    #CHB1
		'GCST001956': 'EUR',    #EUR5
		'GCST000372': 'EUR',	#EUR6
		'GCST000174': 'EUR',	#EUR7
		'GCST000175': 'EUR',	#EUR8
		'GCST000176': 'EUR',	#EUR9
		'GCST006803': 'EUR',    #EUR1_S
		'GCST003880': 'CHB',    #CHB1_S
		'GCST001242': 'EUR',	#EUR2_S
		'GCST007205': 'JPT',	#JPT1_S
		'GCST001757': 'EUR',	#EUR3_S
		'GCST002149': 'EUR',	#EUR4_S
		'GCST001657': 'EUR',	#EUR5_S
		'GCST009004': 'EUR',	#EUR1_B
		'GCST009001': 'EUR',	#EUR1c_B
		'GCST009003': 'EUR',	#EUR1B_B
		'GCST002461': 'EAS',	#EAS1_B
		'GCST006369': 'EAS',	#EAS2_B
		'GCST002021': 'EUR',	#EUR2_B
		'GCST005951': 'EUR',	#EUR3_B
		'GCST008158': 'EUR',	#EUR4_B
		'GCST006802': 'GBR',	#GBR1_B
		'GCST008114': 'AFR',	#AFR1_T2D
		'GCST007847': 'JPT',	#JPT1_T2D
		'GCST010118': 'EAS',	#EAS1_T2D
		'GCST009379': 'EUR',	#EUR1_T2D
		'GCST006801': 'GBR',	#GBR1_T2D
		'GCST004773': 'EUR',	#EUR2_T2D
		'GCST005047': 'EUR',	#EUR3_T2D
		'GCST002560': 'ASW',	#ASW1_T2D
		'GCST002128': 'EAS',	#EAS2_T2D
		'GCST001351': 'EAS',	#EAS3_T2D
		'GCST90013693': 'JPT',	#JPT2_T2D
		'GCST007517': 'EUR',	#EUR4_T2D
		'GCST010556': 'ASW',    #ASW2_T2D
		'GCST000712': 'EUR',	#EUR5_T2D
		'GCST005786': 'AFR',	#AFR1_PC
		'GCST006085': 'EUR',    #EUR1_PC
		'GCST009336': 'EAS',    #EAS1_S
		'GCST001565': 'EUR',	#EUR6_S
	}

	# load data and add scaled score
	all_dfs = []
	std_scaler = StandardScaler()

	for s in studies:
		df = pd.read_csv(os.path.join(save_dir, '{}_scores.csv'.format(s)))

		if log_score:
			df['scores'] = np.log(df.scores)

		df['scaled_score'] = std_scaler.fit_transform(df.scores.values.reshape(-1,1))
		df['gender_scaled_score'] = 0
		df.loc[df.gender == 'male', 'gender_scaled_score'] = std_scaler.fit_transform(
			df.scores[df.gender == 'male'].values.reshape(-1,1)
		)
		df.loc[df.gender == 'female', 'gender_scaled_score'] = std_scaler.fit_transform(
			df.scores[df.gender == 'female'].values.reshape(-1,1)
		)

		df['sample_pop'] = study_to_sample_pop[s]
		df['study'] = df[['STUDY ACCESSION', 'sample_pop']].agg(' '.join, axis=1)

		all_dfs.append(df)

	df = pd.concat(all_dfs)

	with sns.plotting_context(context='notebook', font_scale=1.):

		# For log odds scores
		sns.displot(data=df, x='scaled_score', hue='super_pop', 
					col='study', kind='kde', 
					facet_kws={'margin_titles':True})
		plt.show()

		# sns.displot(data=df, x='gender_scaled_score', hue='super_pop', 
		#             row='gender', col='study', kind='kde', 
		#             facet_kws={'margin_titles':True})
		# plt.show()

		# sns.displot(data=df, x='gender_scaled_score', hue='super_pop', 
		#             row='sample_pop', col='gender', kind='kde',
		#             facet_kws={'margin_titles':True})
		# plt.show()

		# sns.displot(data=df, x='gender_scaled_score', hue='sample_pop', 
		#             row='gender', col='super_pop', kind='kde',
		#             facet_kws={'margin_titles':True})
		# plt.show()

		# sns.displot(data=df, x='gender_scaled_score', hue='study', 
		#     row='gender', col='super_pop', kind='kde',
		#     facet_kws={'margin_titles':True})
		# plt.show()

	# show north/southern europe height difference
