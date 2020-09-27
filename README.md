# MHCshrubs
*MHC Superimposed Hierarchical Relations using Bayesian Statistics*

<b> /!\ This repository is under development /!\ </b>

This is a Python package for estimating HLA disease associations,
using the functional similarities between HLA molecules as prior information to aid the discovery of these associations.
We have described the algorithm and applications in the following preprint

> [Christiaan H. van Dorp and Can Kesmir](https://doi.org/10.1101/408302) Estimating HLA disease associations using similarity trees, bioRxiv preprint (2018)

## Installation

### Dependencies

* Clustal Omega (optional)
* JAGS
* NetMHCpan (optional)

Python packages

* numpy
* scipy
* matplotlib
* pystan (optional)
* tqdm
* networkx (version >= 2.0)
* ete3
* PyQt5 (dependency of ete3)

### Installation instructions

#### 0. Install JAGS.
On Ubuntu, one can install JAGS using `apt`:

```bash
$ sudo apt-get install jags
```

This will automatically make the `jags` command available system-wide.
For other OSs, follow the instructions on
[the JAGS website](http://mcmc-jags.sourceforge.net/).

#### 1. Clone the git repository
Either download the repository manually from github, or in the terminal type:

```bash
$ git clone git@github.com:chvandorp/MHCshrubs.git
```

#### 2. Install the python package (optional)
The package can be used either directly, or can be installed locally using `pip3`
If you don't install the package, the command line tool can be executed as follows:

```bash
$ cd /path/to/MHCshrubs
$ python3 -m mhcshrubs --help ## prints help message
```

If you choose to install, we advice installation in a `virtualenv` using

```bash
$ pip3 install virtualenv ## can be skipped if virtualenv is already installed
$ cd ~/path/to/work/dir ## choose a folder to work from
$ virtualenv shrubsenv ## create the virtualenv (choose any name you like)
$ source shrubsenv/bin/activate ## activates the shrubsenv
(shrubsenv) $ pip3 install /path/to/MHCshrubs ## installs the mhcshrubs package
(shrubsenv) $ mhcshrubs --help ## we can now use the command line tool
```

## Preparing data
The required data depends on the chosen options.
A table with MHC data and the trait of interest is always required.
In addition, the user can provide e.g. a list of MHC pseudo-sequences,
a `fasta` file with a proteome, a table with allele frequencies, etc.

### Meta-data using a `json` file
In order to run the program, the user must prepare a `json` file with some
basic information. An example is given by `example.json`, which contains the
following data:

```python
{
  "id" : "DurbanHIV1", ## used for output file names
  "outputFolder" : "data/durban", ## directory for writing output files
  "subjectFileName" : "data/durban/subjects_resolved_top0.99_all.tsv", ## file with subject data
  "alleleFreqFileName" : "data/hla/mhc-top0.99counts-ncbi-Sub-Saharan-Africa.tsv", ## file with allele counts
  "pSeqFileName" : "data/hla/MHC_pseudo.dat", ## file with MHC pseudo sequences
  "fastaFileName" : "data/hiv/proteome-clade-C/Ref.C.ZA.04.04ZASK146.AY772699.fasta", ## file with a pathogen's proteome
  "traitFieldName" : "VL", ## the field name of the trait of interest in the subject file
  "traitType" : "continuous", ## the type of trait: continuous or categorical
  "alleleFieldNames" : { ## the field names of the MHC loci in the subject file
    "A" : ["HLA_A1", "HLA_A2"],
    "B" : ["HLA_B1", "HLA_B2"],
    "C" : ["HLA_Cw1", "HLA_Cw2"]
  }
}
```

### MHC genotype and trait data from subjects

The most important file provided by the user contains the subject data. The file name must be given in the `json` meta-data file at the key `"subjectFileName"`. the subject file should be in the
`tsv` of `csv` format and each row should contain data for one
subject. Typically, the trait value of interest (e.g. viral load, or infection status), and HLA alleles at
a number of loci. The user also can provide optional censoring
information about the trait value. In the `json` file the column
name in the header of the subject file must be given at key
`"traitFieldName"`. If the trait field name is `X`,
then the censoring information should be given in a column with
name `X_censoring`. If no such column exists, all values are
assumed to be uncensored.

#### Censoring codes

The censoring codes for trait values are:
* `uncensored` for data that are not censored
* `left_censored` for left-censored data
* `right_censored` for right-censored data
* `missing` for when the observation is completely missing

For categorical data, only the codes `uncensored` and `missing` are allowed.

#### Admissible MHC alleles

The subject file should contain columns for HLA alleles. The user
can provide alleles at various loci (e.g. A, B, C), with two alleles per locus.
The names of the fields in the subject data file should be given at the `"alleleFieldNames"` key,
as indicated in the example above.
In order to only use a single locus (e.g. B) as predictor, leave out any other allele field names
in the `json` file. For instance:

```python
"alleleFieldNames" : { ## the field names of the MHC loci in the subject file
  "B" : ["HLA_B1", "HLA_B2"]
}
```

When the full 2 field type is not fully resolved (or an allele is completely missing), it is possible to enter multiple admissible alleles for one subject. These should be added in a single cell of the table, separated by semi-columns.

#### Covariates

If you want to include covariates in the regression, these have to
be present in the subject data file, and have to be specified in the
`.json` metadata file. For instance, if we have a continuous covariate `age` and a binary covariate `male_sex`, then we add a list to the metadata with key `covariateFieldNames`

```python
{
  "id" : "DurbanHIV1",
  ## other items as in example.json
  ## ...
  "covariateFieldNames" : ["age", "male_sex"] ## covariates included in the regression
}
```

Missing values are currently not supported, and covariates are only implemented for the categorical trait models. Only real values are allowed (but notice that this includes binary covariates, which must be encoded as 0 and 1). **/!\ TODO /!\**

#### Example

Here's a small example of a possible subject data table.

| `ID`    | `age` | `VL`       | `VL_censoring` | `HLA_A1`      | `HLA_A2` | ... |
|---------|-----|------------|----------------|---------------|----------|-----|
| `1003M` | `32`  | `286000.0` | `uncensored`   | `HLA-A*02:01` | `HLA-A*29:01;HLA-A*29:02;HLA-A*29:03` | ... |
| `1097M` | `53`  | `200`      | `left_censored`| `HLA-A*30:02` | `HLA-A*33:01;HLA-A*33:03` | ... |

## Command-line tool

Assuming the command line tool is installed with `pip3`, we can start a job
from any working folder. This working folder should contain the data files
as indicated in the `json` meta-data file.

### Options

Here we show some basic options. Further instructions can be found in the
manual (<b> /!\ TODO /!\ </b>)
Using the `--help` flag, we can get a list of command-line options

```bash
$ source /path/to/shrubsenv/bin/activate
(shrubsenv) $ mhcshrubs --help
```

To run the tree-model analysis described in the manuscript, we use the options
```bash
(shrubsenv) $ mhcshrubs -j metadata.json --chain-len 10000 --num-chains 4
```
This will run 4 MCMCs of length 10000. By default, the first half of the chains
is thrown away (burn-in). The default prior distribution of the branch weights
is the normal distribution. In order to use the Laplace (double exponential)
distribution instead, use the option `--prior dexp`.
The default MHC similarity measure is based on pseudo-sequences, using the
`PMBEC` amino-acid similarity matrix.
To see if the tree model is better than the model with independent MHC alleles,
we can use the trivial MHC similarity using the option `--mhc-similarity trivial`.
To find out if there is any significant MHC association with the trait of interest,
the user can run the null model by selecting the option `--model null`.
