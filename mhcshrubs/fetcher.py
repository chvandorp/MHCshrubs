"""@package fetcher
Functions for loading data

Import csv files with MHC alleles and disease traits
"""

from __future__ import (print_function, division)
from builtins import (map, str, zip, range)
import numpy as np
import csv
import itertools
import os
import warnings

from mhcshrubs import auxiliary as aux
from mhcshrubs import definitions as defn
from mhcshrubs import mhctools

def importAlleleFrequencyData(fileName): ## "mhc-top0.99counts-ncbi-Sub-Saharan-Africa.tsv"
    """
    input: name of the file with frequency data.
    output: a dictionary with HLA counts
    """
    with open(fileName, 'r') as fileHandle:
        table = [row.split('\t') for row in fileHandle.read().split('\n') if row != '']
    header = table[0]
    table = table[1:]
    countDict = dict((mhctools.MhcObject(row[0]), int(row[1])) for row in table)
    return countDict


def importSubjectData(fileName, traitFieldName, alleleFieldNames, subset=None,
                      verbose=False, traitTransform=(lambda x: x)):
    """
    Import MHC data and continuous disease traits.

    Args:
        fileName (str): name of the tsv/csv file with allele and trait data
            per subject
        TraitFieldName (str): the name of the trait of interest in the header
            of the tsv file
        alleleFieldNames (dict): a dictionary with the field names for the
            alleles.

    Kwargs:
        subset (int): take a sample from the complete data.
            If None (default), use all data.
        verbose (bool): print some basic statistics on the imported data.
        traitTransform (function): a funcion to transform the trait value.
            For instance numpy.log, or the identity (default)

    Returns:
        A dictionary containing
            - TraitValues: the trait values
            - CensCodes: censoring codes for the TraitValues
            - CensBounds: upper (or lower) bounds for left (or right) censored Trait
            - AlleleVecs: admissible alleles
            - Alleles: the alleles in the dataset, in the right order for the allele vectors

    Todo:
        - Intelligently handle censoring codes
        - Give a good discription of the expected file format
    """
    ## detect if the file is tsv or csv
    file_base, file_extension = os.path.splitext(fileName)
    if file_extension == ".csv":
        delim = ','
    elif file_extension == ".tsv":
        delim = '\t'
    else:
        raise Exception("invalid file format")
    ## read the contents of the file
    with open(fileName, 'r') as fileHandle:
        reader = csv.DictReader(fileHandle, delimiter=delim)
        data = [row for row in reader]
        ## FIXME: quick bodge for CMV: filter out negatives
        if traitFieldName == "CMV":
            data = [row for row in data if np.log10(float(row[traitFieldName])) > -0.3 ]
    ## take a sub sample
    if subset is not None:
        idxs = np.random.choice(len(data), size=subset, replace=False)
        data = [data[i] for i in idxs]
    ## find data in the table
    Ploidy = 2
    loci = sorted(alleleFieldNames.keys())
    Alleles = {}
    for locus in loci:
        alls = [row[alleleFieldNames[locus][i]].split(';') for i in range(Ploidy) for row in data]
        ## make MhcObject-s
        alls = [mhctools.MhcObject(a) for a in aux.unique(aux.flatten(alls))]
        ## remove non-expressed alleles and reduce to 2-field
        alls = [a.protein for a in alls] ## include Null alleles etc.
        ## find unique elements again
        Alleles[locus] = aux.unique(alls)

    ## for each subject, make a vector with zeros and ones
    ## indicating whether the allele is present
    def hasAllele(subject, allele, locus, i):
        fieldName = alleleFieldNames[locus][i]
        admissibles = [mhctools.MhcObject(a) for a in subject[fieldName].split(';')]
        in_adms = [adm_allele in allele for adm_allele in admissibles]
        return any(in_adms)

    AlleleVecs = {}
    for locus in sorted(alleleFieldNames.keys()):
        AlleleVecs[locus] = [[[hasAllele(subject, allele, locus, i)
                               for allele in Alleles[locus]]
                              for i in range(Ploidy)]
                             for subject in data]

    ## get the virus load and virus load censoring information
    traitCensFieldName = "{0}_censoring".format(traitFieldName)
    CensCodes = np.array([subject[traitCensFieldName]
                          if traitCensFieldName in subject.keys()
                          else defn.uncensored_code ## FIXME, make sure that this is consistant
                          for subject in data])

    TraitValues = np.array([traitTransform(float(subject[traitFieldName]))
                            if cens_code == defn.uncensored_code
                            else np.nan
                            for cens_code, subject in zip(CensCodes, data)])

    ## find the upper or lower bound for the left and (resp.) right censored VLs.
    ## FIXME For the non-interval-censored data, use an arbitrary (low) Trait value as lower bound
    CensBounds = np.array([traitTransform(float(subject[traitFieldName]))
                           if cens_code == defn.left_censored_code
                           or cens_code == defn.right_censored_code
                           else defn.auxiliaryLowerCensBound
                           for cens_code, subject in zip(CensCodes, data)])
    if verbose:
        ## number op subjects
        print("number of subjects: {0}".format(len(data)))
        ## Trait statistics
        mtrait = np.nanmedian(TraitValues)
        ltrait = np.nanpercentile(TraitValues, 2.5)
        htrait = np.nanpercentile(TraitValues, 97.5)
        print("median trait value: {0:0.2f}, 2.5 - 97.5 percentiles: {1:0.2f} - {2:0.2f}".format(mtrait, ltrait, htrait))
        ## allele statistics
        print("total number of alleles:", np.sum([len(Alleles[locus]) for locus in loci]))
        for locus in sorted(Alleles.keys()):
            print("number of {0} alleles: {1}".format(locus, len(Alleles[locus])))
        numCompleteHaplotypes = 0
        for idx in range(len(data)):
            compl = [sum(1 for x in alleleVec if x) == 1
                     for locus in sorted(Alleles.keys())
                     for alleleVec in AlleleVecs[locus][idx]]
            if all(compl):
                ## subject has with 2-field typed haplotype
                numCompleteHaplotypes += 1
        print("number of complete haplotypes: {0}".format(numCompleteHaplotypes))

    ## return relevant data in a dictionary
    dataDict = {
        "TraitValues" : TraitValues,
        "CensCodes" : CensCodes,
        "CensBounds" : CensBounds,
        "AlleleVecs" : AlleleVecs,
        "Alleles" : Alleles
    }
    return dataDict

def importBinarySubjectData(fileName): ## "hcv/subjects_resolved_nmp_all.tsv"
    """
    input: name of the tsv file with patient HLA and health status data
    output:
    - patientStatus -- the health status
    - StatusCensCodes -- censoring codes for health status
    - patientAlleleVecs -- admissible HLA alleles
    - hlaAlleles -- the HLA alleles in the dataset, in the right order

    FIXME: repair this function and the input file
    """
    data, header = aux.basicTsvImport(fileName, has_header=True)
    dataDicts = [dict(list(zip(header, row))) for row in data]

    ## get all HLA alleles in the data set.
    hlaKeys = dict((X, ["HLA_{}{}_resolved".format(X,i) for i in [1,2]]) for X in "ABC")

    hlaAlleles = {}
    for X in "ABC":
        alleles = [row[hlaKeys[X][i]].split(';') for i in range(2) for row in dataDicts]
        alleles = aux.unique(aux.flatten(alleles))
        hlaAlleles[X] = [mhctools.MhcObject(a) for a in alleles]

    ## for each patient, make a vector with zeroes and ones
    ## indicating whether the allele is present
    def hasAllele(p, hla, X, i):
        admissibles = [mhctools.MhcObject(a) for a in p[hlaKeys[X][i]].split(';')]
        return hla in admissibles ## "in" uses the __eq__ class method to determine equality

    patientAlleleVecs = {}
    for X in "ABC":
        patientAlleleVecs[X] = [[[hasAllele(p, h, X, i)
            for h in hlaAlleles[X]] for i in range(2)] for p in dataDicts]

    ## get the health status
    patientHealthStatuss = np.array([int(row["HCV_clr"]) for row in dataDicts]) ## TODO/FIXME: make this generic

    ## return relevant data
    return (patientHealthStatuss, patientAlleleVecs, hlaAlleles)
