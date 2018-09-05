# MHCshrubs
*MHC Superimposed Hierarchical Relations using Bayesian Statistics*

This is a Python package for estimating HLA disease associations, 
using the functional similarities between HLA molecules as prior information to aid the discovery of these associations.
We have described the algorithm and applications in the following preprint

> [Christiaan H. van Dorp and Can Kesmir](https://doi.org/10.1101/408302) Estimating HLA disease associations using similarity trees, bioRxiv preprint (2018)

## Installation

*instructions will follow soon*

### Dependencies

* Clustal Omega (optional)
* JAGS
* NetMHCpan
* pystan (optional)
* tqdm
* networkx (version >= 2.0)
* ete3

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
  "alleleFieldNames" : { ## the field names of the MHC loci in the subject file
    "A" : ["HLA_A1", "HLA_A2"],
    "B" : ["HLA_B1", "HLA_B2"],
    "C" : ["HLA_Cw1", "HLA_Cw2"]
  }
}
```
