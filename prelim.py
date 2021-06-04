import numpy as np
import pandas as pd

import allel


if __name__ == '__main__':
	# height
	study_path = 'data/GWAS_Catalog/gwas-catalog-v1.0.3-studies-r2021-05-05.tsv'
	assoc_path = 'data/GWAS_Catalog/gwas_catalog_v1.0.2-associations_e100_r2021-05-05.tsv'
	ancest_path = 'data/GWAS_Catalog/gwas-catalog-v1.0.3-ancestries-r2021-05-05.tsv'

	all_study_df = pd.read_csv(study_path, delimiter='\t')
	assoc_df = pd.read_csv(assoc_path, delimiter='\t')
	ancest_df = pd.read_csv(ancest_path, delimiter='\t')

	ancest_df['n_pop'] = ancest_df.groupby(['STUDY ACCESSION','STAGE']).transform('count').max(level=0).DATE

	pd.set_option('display.max_rows', 500)
	# print(study_df['DISEASE/TRAIT'].value_counts().head(400))

	hs_df = all_study_df[all_study_df['DISEASE/TRAIT'] == 'Height']

	prev_cols = ['DATE', 'JOURNAL',
		'LINK', 'STUDY', 'DISEASE/TRAIT', 'INITIAL SAMPLE SIZE',
		'REPLICATION SAMPLE SIZE',
		'ASSOCIATION COUNT', 'MAPPED_TRAIT', 'MAPPED_TRAIT_URI',
		'STUDY ACCESSION']

	prev_cols_min = ['DATE', 'DISEASE/TRAIT', 'INITIAL SAMPLE SIZE',
		'REPLICATION SAMPLE SIZE',
		'ASSOCIATION COUNT', 'MAPPED_TRAIT']

	hs_df = hs_df[prev_cols]
	hs_df_min = hs_df[prev_cols_min]

	joined_cols = ['DATE_x', 'JOURNAL', 'LINK', 'STUDY', 'DISEASE/TRAIT',
       'INITIAL SAMPLE SIZE', 'REPLICATION SAMPLE SIZE', 'ASSOCIATION COUNT',
       'MAPPED_TRAIT', 'MAPPED_TRAIT_URI', 'STUDY ACCESSION', 'PUBMED ID',
       'FIRST AUTHOR', 'INITIAL SAMPLE DESCRIPTION',
       'REPLICATION SAMPLE DESCRIPTION', 'STAGE', 'NUMBER OF INDIVIDUALS',
       'BROAD ANCESTRAL CATEGORY', 'COUNTRY OF ORIGIN',
       'COUNTRY OF RECRUITMENT', 'n_pop']
	
	joined_cols_min = ['DATE_x', 'LINK', 'STUDY', 'DISEASE/TRAIT',
       'INITIAL SAMPLE SIZE', 'REPLICATION SAMPLE SIZE', 'ASSOCIATION COUNT',
       'MAPPED_TRAIT', 'STUDY ACCESSION', 'INITIAL SAMPLE DESCRIPTION',
       'REPLICATION SAMPLE DESCRIPTION', 'STAGE', 'NUMBER OF INDIVIDUALS',
       'BROAD ANCESTRAL CATEGORY', 'COUNTRY OF ORIGIN',
       'COUNTRY OF RECRUITMENT', 'n_pop']

	all_jdf = pd.merge(all_study_df, ancest_df, how='left', on='STUDY ACCESSION')
	jdf = pd.merge(hs_df, ancest_df, how='left', on='STUDY ACCESSION')
	mjdf = jdf[joined_cols_min]

	one_samp_pop_df = mjdf[mjdf.n_pop == 1]

	is_euro = one_samp_pop_df['BROAD ANCESTRAL CATEGORY'] == 'European'
	edf = one_samp_pop_df[is_euro]
	nedf = one_samp_pop_df[~is_euro]

	# get joined data for study
	study_id = 'GCST001956'
	study_df = jdf[jdf['STUDY ACCESSION'] == study_id]
	ss_df = assoc_df[assoc_df['STUDY ACCESSION'] == study_id].iloc[:, 11:-5]
	ss_df = ss_df[ss_df['OR or BETA'].notnull()]

	# # open vcf
	# chromosome = 19
	# path = 'data/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr{}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz'.format(chromosome)
	
	# # callset = allel.read_vcf(path, '*', samples=['HG00096'])
	# print('Loading site dataframe...')
	# vcf_df = allel.vcf_to_dataframe(path, '*')
	# # callset = allel.read_vcf(path, samples=['HG00096'])
	# print('Saving vcf as npz...')
	# allel.vcf_to_npz(path, 'data/preprocessed/chr{}.npz'.format(chromosome), 
	# 					fields='*', overwrite=True)
	# print('Loading npz...')
	# callset = np.load('data/preprocessed/chr{}.npz'.format(chromosome), 
	# 					allow_pickle=True)