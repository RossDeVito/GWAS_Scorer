import os
import sys

import numpy as np
import pandas as pd
from tqdm import tqdm

import allel


def preprocess_to_hdf5(vcf_path, save_dir):
	"""Preprocesses vcf or vcf.gz files for each chromosome to HDF5 files.

	Args:
		vcf_path: String that when '{}' is replaced with {[1,22] and 'X'}
			becomes the path to each chromosome's vcf. 
				e.g.: 'data/chr{}.vcf'
		save_dir: Path to directory in which HFD5 files will be saved with
			the naming convention 'chr{}.h5'. If directory already exists
			will raise error. 

	Raises:
		ValueError when save_dir already exists.
	"""
	# Make save directory if doesn't exist
	try:
		os.mkdir(save_dir)
	except FileExistsError:
		raise ValueError("Save dir must not already exist")

	# Convert vcfs
	print('Saving vcfs as hdf5s')

	to_do = list(range(1, 23)) + ['X']
	to_do = tqdm(to_do)

	for chromosome in to_do:
		to_do.set_description('Chromosome {} '.format(chromosome))
		to_do.refresh()
		
		allel.vcf_to_hdf5(
			vcf_path.format(chromosome),
			os.path.join(save_dir, 'chr{}.h5'.format(chromosome)),
			fields='*',
			alt_number=1, 
			overwrite=True,
			log=sys.stdout
		)


def main():
	vcf_path = 'data/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr{}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz'
	save_dir = 'data/preprocessed/1000_genomes_GRCh38_20181129'

	preprocess_to_hdf5(vcf_path, save_dir)


if __name__ == '__main__':
	main()