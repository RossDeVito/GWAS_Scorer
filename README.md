# GWAS_Scorer
Python functions and command line tools for easily scoring samples using PRS from the NHGRI-EBI GWAS Catalog or other similarly formatted data. Additionallly provides a tool for calculating the fixation indices between all populations. Before samples can be scored or fixation indices calculated, the genotype data must be preprocessed from .vcf or .vcf.gz to HDF5.

Preprocessed data and 58 scores for the 1000 Genomes samples are available [here](https://drive.google.com/file/d/1U9N3YqD5n47g_3i8Gjwe0TKunjw1eFhi/view?usp=sharing). To use, unzip in data directory. The report using this data is [here]().

## 1. Preprocess Data
The samples' genotype data can originally be in .vcf or .vcf.gz format, with one chromosome per file. To improve runtime and use less memory when actually scoring samples or calculating fixation indices, these other functions use HDF5 files which are created and stored by this function.

### Command Line
```
usage: GWAS_Scorer.py preprocess_vcf [-h] vcf_path save_dir

Preprocess vcf or vcf.gz files for each chromosome to hdf5 files required by
scorer.

positional arguments:
  vcf_path    String that when '{}' is replaced with {[1,22] and 'X'} becomes
              the path to each chromosome's vcf. e.g.: 'data/chr{}.vcf'
  save_dir    Path to directory in which HFD5 files will be saved with the
              naming convention 'chr{}.h5'. If directory already exists will
              raise error.

optional arguments:
  -h, --help  show this help message and exit

Example usage: python GWAS_Scorer.py preprocess_vcf $VCFPATH $SAVEDIR
```
### Python
```
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
 ```
example use:
```
from preprocess_to_hdf5 import preprocess_to_hdf5
vcf_path = 'data/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr{}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz'
save_dir = 'data/preprocessed/1000_genomes_GRCh38_20181129'
preprocess_to_hdf5(vcf_path, save_dir)
```

## 2. Score

Generate scores for all samples based on given study or studies.

### Command Line
```
usage: GWAS_Scorer.py score [-h] [-o]
                            assoc_path hdf5_dir panel_path save_dir study_id
                            [study_id ...]

For each study, generates score for all samples in hdf5_dir data and saves
results as csv. Results are saved in save_dir with the format
{study_id}_scores.csv.

positional arguments:
  assoc_path        Path to tsv containing studies' resulting set of sites and
                    weights. Tab seperated value file should follow format of
                    NHGRI-EBI GWAS Catalog 'all associations' data.
  hdf5_dir          Directory containing preprocessed hdf5 data (see
                    preprocess_vcf).
  panel_path        Path to tab seperated value file containing sample
                    population data with the headings sample, pop, super_pop,
                    and gender
  save_dir          Directory where resulting score csv(s) is saved.
  study_id          String identifier for a study or a list of these strings.
                    Each study ID should correspond to a value in the 'STUDY
                    ACCESSION' column of assoc_path's resulting table.

optional arguments:
  -h, --help        show this help message and exit
  -o, --odds_ratio  Whether the 'OR or BETA' feild in assoc_path data should
                    be interpreted as an odds ratio (when -o flag used) or a
                    beta value (when flag not used) when caluclating scores.
                    Odds ratios are typically used for case/control studies
                    and beta values for quantitative traits. Defaults to beta
                    values.

Example usage for case/control studies: python GWAS_Scorer.py score -o
$ASSOCPATH $HDF5DIR $PANELPATH $SAVEDIR GCST001757 GCST009336
```

### Python
```
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
```
example use:
```
from score import get_gwas_scores

study_ids = 'GCST008904' # or something like ['GCST008904', 'GCST036147']
odds_ratio_score = True

assoc_path = 'data/GWAS_Catalog/gwas_catalog_v1.0.2-associations_e100_r2021-05-05.tsv'
hdf5_dir = 'data/preprocessed/1000_genomes_GRCh38_20181129'
panel_path = 'data/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel'
save_dir = 'data/results'

get_gwas_scores(study_ids, assoc_path, hdf5_dir, panel_path, save_dir,
				odds_ratio_score)
```

## 3. Fixation Indices
Calculates pairwise fixation indices between all populations and super populations. Results are saved as JSON and CSV files. For more information on block length selection see [this article](https://digitalcommons.wayne.edu/cgi/viewcontent.cgi?article=1113&context=humbiol_preprints). The general rule when selecting block length is that the length of the blocks should be the smallest value for which increasing the length does not increase the standard error.

### Command Line
```
usage: GWAS_Scorer.py fst [-h] [-b BLOCK_LENGTH] [-z] hdf5_dir panel_path save_path

Calculates fixation index (Fst) between all population pairs. Results are saved as a csv that
also includes the standard_error and as a JSON file used when analyzing gwas scores.

positional arguments:
  hdf5_dir              Directory containing preprocessed hdf5 data (see preprocess_vcf).
  panel_path            Path to tab seperated value file containing sample population data
                        with the headings sample, pop, super_pop, and gender
  save_path             Results will be saved as this with a '.csv' appened for the csv file
                        and '.json' for the json.

optional arguments:
  -h, --help            show this help message and exit
  -b BLOCK_LENGTH, --block_length BLOCK_LENGTH
                        Block size in number of variants to use for block-jackknife when
                        calculating Fst. Defaults to 10000.
  -z, --dont_zero_with_self
                        When flag not used, a population's fixation index with itself is set
                        to zero. When used, Fst with self calculated, but may be negative.

Example usage: python GWAS_Scorer.py fst -b 1000 $HDF5DIR $PANELPATH $SAVEPATH
```

### Python
```
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
```
example usage:
```
from fst import pairwise_fst
hdf5_dir = 'data/preprocessed'
save_path = 'data/preprocessed/Fst_bl10000'
panel_path = 'data/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel'

zero_with_self = True
method='hudson'
block_length=10000

pairwise_fst(hdf5_dir, panel_path, save_path,
				zero_with_self, method,
				block_length)
```
