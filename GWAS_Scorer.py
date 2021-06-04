import argparse
import sys

from preprocess_to_hdf5 import preprocess_to_hdf5
from score import get_gwas_scores
from fst import pairwise_fst


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

    
def run_preprocess(args):
    preprocess_to_hdf5(args.vcf_path, args.save_dir)


def run_score(args):
    get_gwas_scores(args.study_id, args.assoc_path, args.hdf5_dir, 
                    args.panel_path, args.save_dir, args.odds_ratio)


def run_fst(args):
    pairwise_fst(args.hdf5_dir, args.panel_path, args.save_path, 
                    zero_with_self=(not args.dont_zero_with_self), 
                    block_length=args.block_length
                )


parser = MyParser()
subparsers = parser.add_subparsers()

# Preprocess vcf to hdf5
preprocess_parser = subparsers.add_parser(
    'preprocess_vcf',
    description="Preprocess vcf or vcf.gz files for each chromosome to hdf5 files required by scorer.",
    epilog="Example usage: python GWAS_Scorer.py preprocess_vcf $VCFPATH $SAVEDIR"
)
preprocess_parser.add_argument(
	'vcf_path',
	help="String that when '{}' is replaced with {[1,22] and 'X'} becomes the path to each chromosome's vcf. \n\te.g.: 'data/chr{}.vcf'"
)
preprocess_parser.add_argument(
	'save_dir',
	help="Path to directory in which HFD5 files will be saved with the naming convention 'chr{}.h5'. If directory already exists will raise error."
)
preprocess_parser.set_defaults(func=run_preprocess)

# Score
score_parser = subparsers.add_parser(
    'score',
    description="For each study, generates score for all samples in hdf5_dir data and saves results as csv. Results are saved in save_dir with the format {study_id}_scores.csv.",
    epilog="Example usage for case/control studies: python GWAS_Scorer.py score -o $ASSOCPATH $HDF5DIR $PANELPATH $SAVEDIR GCST001757 GCST009336"
)
score_parser.add_argument(
	'assoc_path',
	help="Path to tsv containing studies' resulting set of sites and weights. Tab seperated value file should follow format of NHGRI-EBI GWAS Catalog 'all associations' data."
)
score_parser.add_argument(
    'hdf5_dir',
    help="Directory containing preprocessed hdf5 data (see preprocess_vcf)."
)
score_parser.add_argument(
    'panel_path',
    help="Path to tab seperated value file containing sample population data with the headings sample, pop, super_pop, and gender"
)
score_parser.add_argument(
    'save_dir',
    help="Directory where resulting score csv(s) is saved."
)
score_parser.add_argument(
    'study_id',
    nargs='+',
    help="String identifier for a study or a list of these strings. Each study ID should correspond to a value in the 'STUDY ACCESSION' column of assoc_path's resulting table."
)
score_parser.add_argument(
    '-o', '--odds_ratio', 
    default=False,
    help="Whether the 'OR or BETA' feild in assoc_path data should be interpreted as an odds ratio (when -o flag used) or a beta value (when flag not used) when caluclating scores. Odds ratios are typically used for case/control studies and beta values for quantitative traits. Defaults to beta values.",
    action="store_true"
)
score_parser.set_defaults(func=run_score)

# Fixation indices
fst_parser = subparsers.add_parser(
    'fst',
    description="Calculates fixation index (Fst) between all population pairs. Results are saved as a csv that also includes the standard_error and as a JSON file used when analyzing gwas scores.",
    epilog="Example usage: python GWAS_Scorer.py fst -b 1000 $HDF5DIR $PANELPATH $SAVEPATH"
)
fst_parser.add_argument(
    'hdf5_dir',
    help="Directory containing preprocessed hdf5 data (see preprocess_vcf)."
)
fst_parser.add_argument(
    'panel_path',
    help="Path to tab seperated value file containing sample population data with the headings sample, pop, super_pop, and gender"
)
fst_parser.add_argument(
    'save_path',
    help="Results will be saved as this with a '.csv' appened for the csv file and '.json' for the json."
)
fst_parser.add_argument(
    '-b', '--block_length',
    help="Block size in number of variants to use for block-jackknife when calculating Fst. Defaults to 10000.",
    type=int
)
fst_parser.add_argument(
    '-z', '--dont_zero_with_self',
    default=False,
    help="When flag not used, a population's fixation index with itself is set to zero. When used, Fst with self calculated, but may be negative.",
    action="store_true"
)
fst_parser.set_defaults(func=run_fst)


if __name__ == '__main__':
	'''
	example bash:

	FRAGPATH='data/fragments/chr20_1-1M/fragments.txt'
	VCFPATH='data/fragments/chr20_1-1M/2.0.realigned_genotypes.vcf'

	python GWAS_Scorer.py
	'''
	args = parser.parse_args()
	args.func(args)