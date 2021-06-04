import h5py
import os

import numpy as np
import pandas as pd

import allel


def get_gwas_scores(study_ids, assoc_path, hdf5_dir, panel_path, 
					save_dir=None, odds_ratio_score=False):
	"""For each study, generates score for all samples and saves results as csv.
	
	Generates scores for all samples in hdf5_dir data. When save_dir
	is not None, results are saved in save_dir with the format 
	{study_id}_scores.csv.

	Args:
		study_ids (string or list-like): String identifier for a study or
			a list of these strings. Each study ID should correspond to
			a value in the 'STUDY ACCESSION' column of assoc_path's
			resulting table.
		assoc_path: Path to tsv containing studies' resulting set of sites
			and weights. Tab seperated value file should follow format of
			NHGRI-EBI GWAS Catalog 'all associations' data.
		hdf5_dir: Directory containing preprocessed hdf5 data (see 
			preprocess.py).
		panel_path: Path to tab seperated value file containing sample
			population data in the form:

				sample	pop	super_pop	gender		
				HG00096	GBR	EUR	male
				HG00097	GBR	EUR	female
				NA19794	MXL	AMR	female
				NA19819	ASW	AFR	female
		
		save_dir: Directory where resulting score csv is saved.
		odds_ratio_score (bool, optional): Whether the 'OR or BETA' feild 
			in assoc_path data should be interpreted as an odds ratio 
			(when =True) or a beta value (when =False) when caluclating 
			scores. Odds ratios are typically used for case/control studies 
			and beta values for quantitative traits. Defaults to False.

	Returns:
		List of DataFrames containing same by sample score data that 
		would be saved as csv for each study. 
	"""
	# if just one string study id is input
	if isinstance(study_ids, str):
		study_id = study_ids
		print('\n' + study_id)

		# load all associations
		assoc_df = pd.read_csv(assoc_path, delimiter='\t')

		# get association data for sites significant in the study
		ss_df = assoc_df[assoc_df['STUDY ACCESSION'] == study_id]
		ss_df = ss_df[ss_df['OR or BETA'].notnull()]
		ss_df = ss_df[ss_df.CHR_POS.notnull()]

		ss_df['CHR_POS'] = ss_df.CHR_POS.astype(int)

		ss_df['OR or BETA'] = ss_df.groupby(['CHR_ID', 'CHR_POS'])['OR or BETA'].transform('mean')
		ss_df = ss_df.drop_duplicates(['CHR_ID', 'CHR_POS'], keep='first')

		# Load all data relevant to study
		site_data = None
		gts = []

		sites_w_weights = 0
		sites_used = 0

		print("CHR\tsites in study\tsites used")
		for this_chr in list(range(1, 23)) + ['X']:
			print(this_chr, end='\t')
			
			h5_path = 'data/preprocessed/chr{}.h5'.format(this_chr)
			callset = h5py.File(h5_path, mode='r')

			if not site_data:
				site_data = {}
				for k in callset['variants'].keys():
					site_data[k] = [] 

			# reduces to just sites on chromosome in callset data
			## for positions used twice, uses mean value
			## assumes alt allele is risk when ?
			this_chr_df = ss_df[ss_df.CHR_ID == str(this_chr)]
			this_chr_df = this_chr_df[this_chr_df.CHR_POS.isin(callset['variants']['POS'][:])]

			site_inds = []
			for p in this_chr_df.CHR_POS:
				site_inds.append(np.where(callset['variants']['POS'][:] == p)[0][0])

			# site_inds must be sorted to use with hdf5
			this_chr_df['site_inds'] = site_inds
			this_chr_df = this_chr_df.sort_values('site_inds')
			site_inds = this_chr_df.site_inds.values.astype(int)

			# this_chr_pos = ss_df[ss_df.CHR_ID == str(this_chr)].CHR_POS.astype(int).values
			# site_inds = np.where(np.isin(callset['variants']['POS'][:], this_chr_pos))[0]

			sites_w_weights += len(this_chr_df)
			sites_used += len(site_inds)

			print("{}\t\t{}".format(len(this_chr_df), len(site_inds)))

			# get genotype data
			gts.append(callset['calldata']['GT'][site_inds].sum(axis=2))

			# get site data
			for k in callset['variants'].keys():
				site_data[k].extend(callset['variants'][k][site_inds].tolist())

		print('Found {} of {} sites from study'.format(sites_w_weights, sites_used))

		df = pd.DataFrame(site_data)
		df['ALT'] = df.ALT.str.decode('utf-8')
		df['CHROM'] = df.CHROM.str.decode('utf-8')
		df['ID'] = df.ID.str.decode('utf-8')
		df['REF'] = df.REF.str.decode('utf-8')
		df['VT'] = df.VT.str.decode('utf-8')

		df['CHR_ID'] = df.CHROM.str.strip('chr')
		df['CHR_POS'] = df.POS

		cols_from_join = ['AC', 'AF', 'AFR_AF', 'ALT', 'AMR_AF', 'AN', 'CHROM', 'DP', 'EAS_AF',
		'EUR_AF', 'EX_TARGET', 'FILTER_PASS', 'ID', 'NS', 'POS', 'QUAL', 'REF',
		'SAS_AF', 'VT', 'altlen', 'is_snp', 'numalt', 'CHR_ID', 'CHR_POS',
		'REPORTED GENE(S)', 'MAPPED_GENE', 'UPSTREAM_GENE_ID',
		'DOWNSTREAM_GENE_ID', 'SNP_GENE_IDS', 'UPSTREAM_GENE_DISTANCE',
		'DOWNSTREAM_GENE_DISTANCE', 'STRONGEST SNP-RISK ALLELE', 'SNPS',
		'MERGED', 'SNP_ID_CURRENT', 'CONTEXT', 'INTERGENIC',
		'RISK ALLELE FREQUENCY', 'P-VALUE', 'OR or BETA']

		jdf = pd.merge(df, ss_df, how='left', on=['CHR_ID', 'CHR_POS'])#[cols_from_join]
		jdf['RISK ALLELE'] = np.array(jdf['STRONGEST SNP-RISK ALLELE'].str.split('-').tolist())[:,1]

		# get coeficients from genotypes, adjusting for which allele is risk allele
		coef_mat = np.vstack(gts)
		ref_risk = jdf['RISK ALLELE'] == jdf.REF
		coef_mat[ref_risk] = 2 - coef_mat[ref_risk]

		if odds_ratio_score:
			scores = np.exp(np.log(jdf['OR or BETA'].values) @ coef_mat)
		else:
			scores = jdf['OR or BETA'].values @ coef_mat

		# load and merge population data
		panel_df = pd.read_csv(panel_path, delim_whitespace=True)

		res_df = pd.DataFrame({'sample': callset['samples'][:].astype(str), 'scores': scores})
		res_df = pd.merge(res_df, panel_df, how='left', on='sample')
		res_df['STUDY ACCESSION'] = study_id
		res_df['n_loci_used'] = sites_used
		res_df['n_loci_skipped'] = sites_w_weights - sites_used

		# add by super pop risk allele freq (POP_RAF)
		jdf['AFR_RAF'] = jdf['AFR_AF']
		jdf.loc[ref_risk, 'AFR_RAF'] = 1 - jdf['AFR_RAF'][ref_risk]
		jdf['AMR_RAF'] = jdf['AMR_AF']
		jdf.loc[ref_risk, 'AMR_RAF'] = 1 - jdf['AMR_RAF'][ref_risk]
		jdf['EAS_RAF'] = jdf['EAS_AF']
		jdf.loc[ref_risk, 'EAS_RAF'] = 1 - jdf['EAS_RAF'][ref_risk]
		jdf['EUR_RAF'] = jdf['EUR_AF']
		jdf.loc[ref_risk, 'EUR_RAF'] = 1 - jdf['EUR_RAF'][ref_risk]
		jdf['SAS_RAF'] = jdf['SAS_AF']
		jdf.loc[ref_risk, 'SAS_RAF'] = 1 - jdf['SAS_RAF'][ref_risk]

		# add indicator of if beta vals or odds ratios used to calculate score
		if odds_ratio_score:
			res_df['weight_type'] = 'odds ratio'
		else:
			res_df['weight_type'] = 'beta'

		# Save
		if save_dir:
			res_df.to_csv(os.path.join(save_dir, '{}_scores.csv'.format(study_id)), index=False)
		
		return [res_df]
	else:
		return [
			get_gwas_scores(s_id, assoc_path, hdf5_dir, panel_path, save_dir, odds_ratio_score) for s_id in study_ids
		]	


def main():
	study_ids = 'GCST008904'
	odds_ratio_score = True

	assoc_path = 'data/GWAS_Catalog/gwas_catalog_v1.0.2-associations_e100_r2021-05-05.tsv'
	hdf5_dir = 'data/preprocessed/1000_genomes_GRCh38_20181129'
	panel_path = 'data/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel'
	save_dir = 'data/results'

	get_gwas_scores(study_ids, assoc_path, hdf5_dir, panel_path, save_dir,
					odds_ratio_score)


if __name__ == '__main__':
	main()