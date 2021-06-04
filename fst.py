import time
import os
import h5py
import json
from itertools import combinations_with_replacement

import numpy as np
import pandas as pd
import dask
from tqdm import tqdm

import allel


def recur_dictify(frame):
	"""Converts pandas DataFrame to a nested dict.

	The right most column's values are the innermost values of the nested
	dicts. The values of the other columns are used as keys in the nested
	dictionaries from left to right such that the left most column's values 
	are the keys of the outermost dictionary and the second most right 
	column's values are used as keys for the innermost dicts who's values 
	are values from the rightmost column.
	"""
	if len(frame.columns) == 1:
		if frame.values.size == 1: return frame.values[0][0]
		return frame.values.squeeze()
	grouped = frame.groupby(frame.columns[0])
	d = {k: recur_dictify(g.iloc[:,1:]) for k,g in grouped}
	return d


def pairwise_fst(hdf5_dir, panel_path, save_path,
					zero_with_self=True, method='hudson',
					block_length=10000):
	"""Calculates fixation index (Fst) between all population pairs.

	Results are saved as a csv that also includes the standard_error and
	as a JSON file used when analyzing gwas scores.

	Args:
		hdf5_dir: Directory containing preprocessed hdf5 data (see 
			preprocess.py).
		panel_path: Path to tab seperated value file containing sample
			population data in the form:

				sample	pop	super_pop	gender		
				HG00096	GBR	EUR	male
				HG00097	GBR	EUR	female
				NA19794	MXL	AMR	female
				NA19819	ASW	AFR	female

		save_path: Results will be saved as this with a ".csv" appened
			for the csv file and ".json" for the json.
		zero_with_self: If True, pairwise Fst of pop with self is 
			automatically set to 0 and the standard error to np.nan. If
			False, calculated like normal.
		method: 'hudson' or 'patterson'. Determines F-statistic calculation
			method. See scikit-allel for details.
		block_length: Block size in number of variants to use for 
			block-jackknife when calculating Fst.
	"""
	gts = []
	samps = []

	t0 = time.time()

	print("Loading chromosome data")
	for this_chr in list(range(1, 23)) + ['X']:
		print('\t' + str(this_chr))
		
		h5_path = os.path.join(hdf5_dir, 'chr{}.h5'.format(this_chr))
		callset = h5py.File(h5_path, mode='r')

		# get genotype data
		gts.append(
			allel.GenotypeDaskArray(callset['calldata']['GT'])
		)

		# get site data
		samps.append(callset['samples'][:])

	print('\tCombining')
	all_gts = gts[0].concatenate(gts[1:])

	print("\ttime: {}".format(time.strftime('%H:%M:%S', time.gmtime(time.time() - t0))))

	print("Combined genotype array shape: {}".format(all_gts.shape))

	print("All sampels in same order: {}".format(
		np.all([samps[i-1] == samps[1] for i in range(len(samps))])
	))
	samps = samps[0].astype(str)

	# load and merge population data
	panel_df = pd.read_csv(panel_path, delim_whitespace=True)

	samps_df = pd.DataFrame({'sample': samps})
	samps_df = pd.merge(samps_df, panel_df, how='left', on='sample')

	# get Fst between all pops
	pop_labels = np.append(samps_df['pop'].unique(), samps_df['super_pop'].unique())
	pop_labels = pop_labels[~pd.isna(pop_labels)]

	# dictionary mapping population names to sample indices
	subpops = {}
	for p in pop_labels:
		subpops[p] = samps_df[(samps_df['pop'] == p) | (samps_df['super_pop'] == p)].index

	# allele counts
	print("Getting allele counts")
	t0 = time.time()
	with dask.config.set(**{'array.slicing.split_large_chunks': False}):
		allele_counts = all_gts.count_alleles_subpops(subpops)
	print("\ttime: {}".format(time.strftime('%H:%M:%S', time.gmtime(time.time() - t0))))

	print("Getting allele counts by pop")
	pop_acs = dict()
	t0 = time.time()
	for p in pop_labels:
		print("\t" + p)
		t1 = time.time()
		pop_acs[p] = allel.AlleleCountsArray(allele_counts[p][:, :2])
		print("\t\ttime: {}".format(time.strftime('%H:%M:%S', time.gmtime(time.time() - t1))))
	print("\ttime: {}".format(time.strftime('%H:%M:%S', time.gmtime(time.time() - t0))))

	# for all pops, get pairwise Fsts
	all_dfs = []

	print("Method: {}".format(method))
	print("Block length: {}".format(block_length))

	label_pairs = list(combinations_with_replacement(pop_labels, 2))
	t = tqdm(label_pairs, desc='Pairwise Fst')
	for pop1, pop2 in t:
		t.set_description('Pairwise Fst pops: {} {} '.format(pop1, pop2))
		t.refresh()

		if pop1 == pop2 and zero_with_self:
			fst_df = pd.DataFrame(
				{
					'block_len': [block_length],
					'Fst': [0.],
					'standard_error': [np.nan]
				}
			)
		else:
			ac1 = pop_acs[pop1]
			ac2 = pop_acs[pop2]

			if method == 'hudson':
				fst, se, _, _ = allel.average_hudson_fst(ac1, ac2, block_length)
			elif method == 'patterson':
				fst, se, _, _ = allel.average_patterson_fst(ac1, ac2, block_length) # hudson seems fastest
			else:
				raise ValueError("method must be 'hudson' or 'patterson'")

			fst_df = pd.DataFrame(
				{
					'block_len': [block_length],
					'Fst': [fst],
					'standard_error': [se]
				}
			)

		if pop1 == pop2:
			fst_df['pop_1'] = pop1
			fst_df['pop_2'] = pop2
			all_dfs.append(fst_df)
		else:
			fst_df2 = fst_df.copy()
			fst_df['pop_1'] = pop1
			fst_df['pop_2'] = pop2
			fst_df2['pop_1'] = pop2
			fst_df2['pop_2'] = pop1
			all_dfs.append(fst_df.append(fst_df2))

	res_df = pd.concat(all_dfs)
	res_df.to_csv(save_path + '.csv', index=False)

	with open(save_path + '.json', 'w') as outfile:
		json.dump(recur_dictify(res_df[['pop_1', 'pop_2', 'Fst']]), outfile)


def main():
	"""Runs Fst calculations """

	hdf5_dir = 'data/preprocessed'
	save_path = 'data/preprocessed/Fst_bl10000'
	panel_path = 'data/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel'

	zero_with_self = True
	method='hudson'
	block_length=10000
	# block length info: https://digitalcommons.wayne.edu/cgi/viewcontent.cgi?article=1113&context=humbiol_preprints
	## The length of the blocks is the smallest value for
	## which increasing the length does not increase the standard error (the point at which
	## estimated standard errors converge; Schaefer et al., 2016)""

	pairwise_fst(hdf5_dir, panel_path, save_path,
					zero_with_self, method,
					block_length)


if __name__ == '__main__':
	main()