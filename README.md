# MHCshrubs
*MHC Superimposed Hierarchical Relations using Bayesian Statistics*

**/!\ This repository is under development /!\**

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

0. Install JAGS. On Ubuntu, one can install JAGS using `apt`:

```bash
$ sudo apt-get install jags
```

This will automatically make the `jags` command available system-wide.
For other OSs, follow the instructions on
[the JAGS website](http://mcmc-jags.sourceforge.net/).

1. Clone the git repository

```bash
$ git clone git@github.com:chvandorp/MHCshrubs.git
```

2. The package can be used either directly, or can be installed locally using `pip3`
If you don't install the package, the command line tool can be executed as follows:

```bash
$ cd /path/to/MHCshrubs
$ python3 -m mhcshrubs --help ## prints help message
```

3. We advice installation in a `virtualenv` using

```bash
$ pip3 install virtualenv ## can be skipped if virtualenv is already installed
$ cd ~/path/to/work/dir ## choose a folder to work from
$ virtualenv shrubs ## create the virtualenv (choose any name you like)
$ source shrubs/bin/activate ## activates the shrubs env
(shrubs) $ pip3 install /path/to/MHCshrubs ## installs the mhcshrubs package
(shrubs) $ mhcshrubs --help ## we can now use the command line tool
```

## Preparing data

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
