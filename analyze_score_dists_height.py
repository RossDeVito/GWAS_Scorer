import os
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
import scikit_posthocs as sp
from statsmodels.stats.multicomp import MultiComparison
from statsmodels.stats.multitest import multipletests

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
		'GCST004212',   #GBR1
		'GCST008163',   #EUR3
		'GCST002647',	#EUR4
		'GCST001956',   #EUR5
		'GCST000372',	#EUR6
		'GCST000174',	#EUR7
		'GCST000175',	#EUR8
		'GCST000176',	#EUR9
		'GCST000611',   #JAP1
		'GCST008839',   #JAP2
		'GCST001263',   #AFR1
		'GCST001290',	#AFR2
		'GCST002702',   #EAS1
		'GCST001885',   #CHB1
	]

	height_multi_pop_studies = [
		'GCST008053',	#Non-EUR
		'GCST008904',	#AFR,AMR,EAS,EUR
	]

	height_studies_1 = [
		'GCST001956',   #EUR5
		'GCST001885',   #CHB1
	]

	height_studies_EUR = [
		'GCST005908',   #EUR1
		'GCST000817',   #EUR2
		'GCST004212',   #GBR1
		'GCST008163',   #EUR3
		'GCST002647',	#EUR4
		'GCST001956',   #EUR5
		'GCST000372',	#EUR6
		'GCST000174',	#EUR7
		'GCST000175',	#EUR8
		'GCST000176',	#EUR9
	]

	schizophrenia_studies = [
		'GCST006803',   #EUR1_S
		'GCST001242',	#EUR2_S
		'GCST001757',	#EUR3_S
		'GCST002149',	#EUR4_S
		'GCST001657',	#EUR5_S
		'GCST001565',	#EUR6_S
		'GCST009336',   #EAS1_S
		'GCST007205',	#JPT1_S
		'GCST003880',   #CHB1_S
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

	studies = height_studies + height_multi_pop_studies

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
		'GCST008053': 'Non-EUR',	#Non-EUR
		'GCST008904': 'AFR,AMR,EAS,EUR',	#AFR,AMR,EAS,EUR
	}

	study_to_sample_super_pop = {
		'GCST005908': 'EUR',    #EUR1
		'GCST000817': 'EUR',    #EUR2
		'GCST000611': 'EAS',    #JAP1
		'GCST008839': 'EAS',    #JAP2
		'GCST001263': 'AFR',    #AFR1
		'GCST001290': 'AFR',	#AFR2
		'GCST002702': 'EAS',    #EAS1
		'GCST004212': 'EUR',	#GBR1
		'GCST008163': 'EUR',    #EUR3
		'GCST002647': 'EUR',    #EUR4
		'GCST001885': 'EAS',    #CHB1
		'GCST001956': 'EUR',    #EUR5
		'GCST000372': 'EUR',	#EUR6
		'GCST000174': 'EUR',	#EUR7
		'GCST000175': 'EUR',	#EUR8
		'GCST000176': 'EUR',	#EUR9
		'GCST006803': 'EUR',    #EUR1_S
		'GCST003880': 'EAS',    #CHB1_S
		'GCST001242': 'EUR',	#EUR2_S
		'GCST007205': 'EAS',	#JPT1_S
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
		'GCST006802': 'EUR',	#GBR1_B
		'GCST008114': 'AFR',	#AFR1_T2D
		'GCST007847': 'EAS',	#JPT1_T2D
		'GCST010118': 'EAS',	#EAS1_T2D
		'GCST009379': 'EUR',	#EUR1_T2D
		'GCST006801': 'EUR',	#GBR1_T2D
		'GCST004773': 'EUR',	#EUR2_T2D
		'GCST005047': 'EUR',	#EUR3_T2D
		'GCST002560': 'AFR',	#ASW1_T2D
		'GCST002128': 'EAS',	#EAS2_T2D
		'GCST001351': 'EAS',	#EAS3_T2D
		'GCST90013693': 'EAS',	#JPT2_T2D
		'GCST007517': 'EUR',	#EUR4_T2D
		'GCST010556': 'AFR',    #ASW2_T2D
		'GCST000712': 'EUR',	#EUR5_T2D
		'GCST005786': 'AFR',	#AFR1_PC
		'GCST006085': 'EUR',    #EUR1_PC
		'GCST009336': 'EAS',    #EAS1_S
		'GCST001565': 'EUR',	#EUR6_S
		'GCST008053': 'Non-EUR',	#Non-EUR
		'GCST008904': 'AFR,AMR,EAS,EUR',	#AFR,AMR,EAS,EUR
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
		df['sample_super_pop'] = study_to_sample_super_pop[s]
		df['study'] = df[['STUDY ACCESSION', 'sample_pop']].agg(' '.join, axis=1)

		all_dfs.append(df)

	df = pd.concat(all_dfs)

	# Remove samples with null pop labels
	df = df.dropna(subset=['pop', 'super_pop'])

	# Add Fst
	df['pop_Fst'] = df.apply(lambda s: 
		0. if (',' in s.sample_pop) or ('-' in s.sample_pop) else Fst_dict[ s['pop'] ][ s['sample_pop'] ], 
		axis=1
	)
	df.loc[df['pop_Fst'] < 0, 'pop_Fst'] = 0.
	df['super_pop_Fst'] = df.apply(lambda s: 
		0. if (',' in s.sample_pop) or ('-' in s.sample_pop) else Fst_dict[ s['super_pop'] ][ s['sample_pop'] ], 
		axis=1
	)
	df.loc[df['super_pop_Fst'] < 0, 'super_pop_Fst'] = 0.

	# Add in study variance by raw and std score
	df['pop_scaled_score_variance'] = df.groupby(['STUDY ACCESSION', 'pop']).scaled_score.transform('var')
	df['pop_scaled_score_variance_bg'] = df.groupby(['STUDY ACCESSION', 'pop', 'gender']).scaled_score.transform('var')


	# Box plots
	gdf = df.groupby(['pop', 'STUDY ACCESSION', 'sample_super_pop', 'super_pop']).mean().reset_index()
	# gdf.loc[gdf['sample_pop'] != 'EUR', 'sample_pop'] = 'EAS'
	sns.catplot(x="super_pop", y="pop_scaled_score_variance",
                hue="sample_super_pop",
                data=gdf, kind="box")
	plt.show()

	# Correlations
	study_res = {
		'STUDY ACCESSION': [],
		'study': [],
		'sample_pop': [],
		'n_loci_used': [],
		'B-F stat': [],
		'B-F p-val': [],
		'Fst & variance corr': [],
		'Fst & variance corr p-val': [],
		'M-W stat': [],
		'M-W p-val': [],
	}

	# # Across all studies using super pops
	# pop_scores_df = df.groupby('super_pop')['scaled_score'].apply(list).reset_index()
	# pop_labels = pop_scores_df.loc[:, 'super_pop'].values
	# pop_scores = pop_scores_df.loc[:, 'scaled_score'].values
	# pop_scores = [np.array(score_list) for score_list in pop_scores]

	# # Kruskal–Wallis test
	# print(stats.kruskal(*pop_scores))
	# # Conover-Iman posthoc
	# r = sp.posthoc_conover(df, val_col='scaled_score', group_col='super_pop', p_adjust='bonferroni')
	# sp.sign_plot(
	# 	r, 
	# 	xticklabels=True, 
	# 	yticklabels=True, 
	# 	cmap = ['1', '#9BCE8D', '#A43D47', '#EE856D', '#F5C0A3']
	# )
	# plt.show()

	# # Across all studies with EUR sample pop using super pops
	# pop_scores_df = df[df.sample_pop == 'EUR'].groupby('super_pop')['scaled_score'].apply(list).reset_index()
	# pop_labels = pop_scores_df.loc[:, 'super_pop'].values
	# pop_scores = pop_scores_df.loc[:, 'scaled_score'].values
	# pop_scores = [np.array(score_list) for score_list in pop_scores]

	# # Kruskal–Wallis test
	# print(stats.kruskal(*pop_scores))
	# # Conover-Iman posthoc
	# r = sp.posthoc_conover(df, val_col='scaled_score', group_col='super_pop', p_adjust='bonferroni')
	# sp.sign_plot(
	# 	r, 
	# 	xticklabels=True, 
	# 	yticklabels=True, 
	# 	cmap = ['1', '#9BCE8D', '#A43D47', '#EE856D', '#F5C0A3']
	# )
	# plt.show()

	# # Across all studies with non-EUR sample pop using super pops
	# pop_scores_df = df[df.sample_pop != 'EUR'].groupby('super_pop')['scaled_score'].apply(list).reset_index()
	# pop_labels = pop_scores_df.loc[:, 'super_pop'].values
	# pop_scores = pop_scores_df.loc[:, 'scaled_score'].values
	# pop_scores = [np.array(score_list) for score_list in pop_scores]

	# # Kruskal–Wallis test
	# print(stats.kruskal(*pop_scores))
	# # Conover-Iman posthoc
	# r = sp.posthoc_conover(df, val_col='scaled_score', group_col='super_pop', p_adjust='bonferroni')
	# sp.sign_plot(
	# 	r, 
	# 	xticklabels=True, 
	# 	yticklabels=True, 
	# 	cmap = ['1', '#9BCE8D', '#A43D47', '#EE856D', '#F5C0A3']
	# )
	# plt.show()

	# By study
	n_unique = len(df['STUDY ACCESSION'].unique())
	fig, axs = plt.subplots(3, n_unique)
	fig.suptitle("Height GWAS Scores")

	for i, study in enumerate(df['STUDY ACCESSION'].unique()):
		study_res['STUDY ACCESSION'].append(study)
		study_df = df[df['STUDY ACCESSION'] == study]
		study_res['study'].append(study_df.iloc[0].study)
		study_res['sample_pop'].append(study_df.iloc[0].sample_pop)
		study_res['n_loci_used'].append(study_df.iloc[0].n_loci_used)

		do_leg = False
		if i == n_unique-1:
			do_leg = True

		eur_hight_df = study_df[study_df['pop'].isin(['IBS', 'TSI', 'FIN', 'CEU', 'GBR'])]
		# Plot distributions
		sns.kdeplot(
			data=eur_hight_df, 
			x='scaled_score', 
			hue='pop',
			common_norm=False,
			ax=axs[0, i],
			legend=do_leg
		)

		if i == n_unique-1:
			legend = axs[0, i].get_legend()
			labels = (x.get_text() for x in legend.get_texts())
			axs[0, i].legend(
				legend.legendHandles,
				labels,
				bbox_to_anchor=(1.1,0.5), 
				loc="center left",
				borderaxespad=0	
			)

		# Using super pops
		pop_scores_df = study_df.groupby('super_pop')['scaled_score'].apply(list).reset_index()
		pop_labels = pop_scores_df.loc[:, 'super_pop'].values
		pop_scores = pop_scores_df.loc[:, 'scaled_score'].values
		pop_scores = [np.array(score_list) for score_list in pop_scores]

		axs[0, i].set_title(study_df.iloc[0].study.replace(' ', '\n'))
		# axs[1, i].set_title("K-W p-val = {:.3}".format(p))

		# Pairwise variance differences
		s, p = stats.levene(*pop_scores)
		study_res['B-F stat'].append(s)
		study_res['B-F p-val'].append(p)

		mc = MultiComparison(data=study_df.scaled_score, groups=study_df.super_pop)
		r = mc.allpairtest(stats.levene, alpha=.01, method='holm')
	
		df1 = pd.DataFrame(r[-1])
		df2 = df1.copy()
		df2['group1'] = df1['group2']
		df2['group2'] = df1['group1']
		r_df = pd.concat((df1,df2)).pivot(index='group1', columns='group2', values='pval_corr').fillna(1.)
		sp.sign_plot(
			r_df, 
			xticklabels=True, 
			yticklabels=True, 
			ax=axs[1, i], 
			cbar_ax_bbox=[0.905, 0.4, 0.02, 0.15], 
			cmap=['1', '#9BCE8D', '#A43D47', '#EE856D', '#F5C0A3']
		)
		# axs[1, i].set_title("B-F p-val = {:.3}".format(p))
		axs[1, i].set_xlabel(None)

		# Plot relationship with Fst
		sns.scatterplot(data=study_df, x='pop_Fst', y='pop_scaled_score_variance',
					hue='super_pop', ax=axs[2, i], legend=False)

		pop_df = study_df.groupby(['pop']).mean().reset_index()
		r, p = stats.pearsonr(pop_df.pop_Fst, pop_df.pop_scaled_score_variance)
		study_res['Fst & variance corr'].append(r)
		study_res['Fst & variance corr p-val'].append(p)

		axs[2, i].set_xlabel("Fst w/ sample pop")
		# axs[2, i].set_title("Pearson r = {:.3}\np-val = {:.3}".format(r, p))

		# Height cline
		eur_hight_df = study_df[study_df['pop'].isin(['IBS', 'TSI', 'FIN', 'CEU', 'GBR'])]
		south_heights = eur_hight_df[eur_hight_df['pop'].isin(['IBS', 'TSI'])].scaled_score.values
		north_heights = eur_hight_df[eur_hight_df['pop'].isin(['FIN', 'CEU', 'GBR'])].scaled_score.values

		s, p = stats.mannwhitneyu(south_heights, north_heights, alternative='less')
		study_res['M-W stat'].append(s)
		study_res['M-W p-val'].append(p)

		if i > 0:
			axs[0, i].set_ylabel(None)
			axs[1, i].set_ylabel(None)
			axs[2, i].set_ylabel(None)
		else:
			axs[1, i].set_ylabel("Equality of Super\nPop. Variances")
			axs[2, i].set_ylabel("Standardized Score\nVariance")

	all_studies_df = pd.DataFrame(study_res)
	all_studies_df['B-F corrected p-val'] = multipletests(
		all_studies_df['B-F p-val'],
		alpha=.01,
		method='holm')[1]
	all_studies_df['M-W corrected p-val'] = multipletests(
		all_studies_df['M-W p-val'],
		alpha=.01,
		method='holm')[1]
	all_studies_df['Fst & variance corr corrected p-val'] = multipletests(
		all_studies_df['Fst & variance corr p-val'],
		alpha=.01,
		method='holm')[1]
	print(all_studies_df)
	plt.show()