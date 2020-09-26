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
    output: a dictionary with HLA counts, of the form
        MhcObject -> int
    """
    with open(fileName, 'r') as fileHandle:
        table = [row.split('\t') for row in fileHandle.read().split('\n') if row != '']
    header = table[0]
    table = table[1:]
    countDict = dict((mhctools.MhcObject(row[0]), int(row[1])) for row in table)
    return countDict

def fetchContinuousValue(fieldName, censFieldName, data):
    """
    Extract values and censoring codes for field with name 'fieldName'
    from data.
    """
    CensCodes = np.array([subject[censFieldName] if censFieldName in subject.keys()
                          else defn.uncensored_code for subject in data])
    CensBounds = np.array([float(subject[fieldName])
                           if cens_code == defn.left_censored_code
                           or cens_code == defn.right_censored_code
                           else defn.auxiliaryLowerCensBound
                           for cens_code, subject in zip(CensCodes, data)])
    Values = np.array([float(subject[fieldName]) if cens_code == defn.uncensored_code
                       else np.nan for cens_code, subject in zip(CensCodes, data)])
    return Values, CensCodes, CensBounds



def importSubjectData(fileName, traitFieldName, alleleFieldNames,
                      covariateFieldNames=None,
                      traitCensFieldName=None, subset=None,
                      verbose=False, traitTransform=(lambda x: x), categorical=False):
    """
    Import MHC data and continuous disease traits.

    Args:
        fileName (str): name of the tsv/csv file with allele and trait data
            per subject
        traitFieldName (str): the name of the trait of interest in the header
            of the tsv file
        alleleFieldNames (dict): a dictionary with the field names for the
            alleles.

    Kwargs:
        covariateFieldNames (list of str): list with covariate names
            to be imported from the data file.
        traitCensFieldName (str or bool): title of the column with censoring values.
            If False, all data is assumed to be uncensored.
            If True (default) the key is assumed to be X_censoring,
            where X is the traitFieldName.
        subset (int): take a sample from the complete data.
            If None (default), use all data.
        verbose (bool): print some basic statistics on the imported data.
        traitTransform (function): a funcion to transform the trait value.
            For instance numpy.log, or the identity (default)
        categorical (bool): True if the trait is categorical, False otherwise (default).
            Use traitTransform to map string to required values (e.g. 0/1)

    Returns:
        A dictionary containing
            - TraitValues: the trait values
            - CensCodes: censoring codes for the TraitValues
            - CensBounds: upper (or lower) bounds for left (or right) censored Trait
                Only included if not categorical
            - Categories: a list of possible traitValues
                Only included if categorical
            - AlleleVecs: admissible alleles
            - Alleles: the alleles in the dataset, in the right order for the allele vectors

    Todo:
        - Intelligently handle censoring codes
        - Give a good discription of the expected file format
        - Allow for covariates
    """
    ## detect if the file is tsv or csv
    file_base, file_extension = os.path.splitext(fileName)
    if file_extension == ".csv":
        delim = ','
    elif file_extension == ".tsv":
        delim = '\t'
    elif file_extension == ".xlsx":
        ## @todo: implement excel files
        raise Exception("xlsx file format not yet implemented")
    else:
        raise Exception("invalid file format (csv or tsv expected)")
    ## read the contents of the file
    with open(fileName, 'r') as fileHandle:
        reader = csv.DictReader(fileHandle, delimiter=delim)
        data = [row for row in reader]
    ## take a sub-sample
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

    ## get the trait values and censoring information
    if type(traitCensFieldName) is bool and traitCensFieldName == True:
        ## use the default key, based on the traitFieldName
        traitCensFieldName = "{0}_censoring".format(traitFieldName)
    ## now, get censoring from data set if traitCensFieldName is str
    if type(traitCensFieldName) is str:
        CensCodes = np.array([subject[traitCensFieldName]
                              if traitCensFieldName in subject.keys()
                              else defn.uncensored_code ## FIXME, make sure that this is consistant
                              for subject in data])
    else: ## assume everything is uncensored
        CensCodes = np.array([defn.uncensored_code for _ in data])

    ## get trait values
    if not categorical:
        TraitValues = np.array([traitTransform(float(subject[traitFieldName]))
                                if cens_code == defn.uncensored_code
                                else np.nan
                                for cens_code, subject in zip(CensCodes, data)])
        ## find the upper or lower bound for the left and (resp.) right censored values
        ## FIXME For the non-interval-censored data, use an arbitrary (low) Trait value as lower bound
        CensBounds = np.array([traitTransform(float(subject[traitFieldName]))
                               if cens_code == defn.left_censored_code
                               or cens_code == defn.right_censored_code
                               else defn.auxiliaryLowerCensBound
                               for cens_code, subject in zip(CensCodes, data)])
    else: ## the trait is categorical: Censoring can only be missing or uncensored
        TraitValues = np.array([traitTransform(subject[traitFieldName])
                                if cens_code == defn.uncensored_code
                                else np.nan
                                for cens_code, subject in zip(CensCodes, data)])
    ## make a list of possible trait values
    Categories = aux.unique([x for x, c in zip(TraitValues, CensCodes)
                             if c == defn.uncensored_code])

    ## import covariates
    if covariateFieldNames is None:
        covariateFieldNames = [] ## empty list
    Covariates = {} ## dictionary indexed by covariate name, empty if no covariates
    for cfn in covariateFieldNames:
        ccfn = "{0}_censoring".format(cfn)
        covVals, covCCs, covCBs = fetchContinuousValue(cfn, ccfn, data)
        Covariates[cfn] = {
            "Values" : covVals,
            "CensCodes" : covCCs,
            "CensBounds" : covCBs
        }

    ## do some printing...
    if verbose:
        ## number op subjects
        print("number of subjects: {0}".format(len(data)))
        ## Trait statistics
        if not categorical:
            mtrait = np.nanmedian(TraitValues)
            ltrait = np.nanpercentile(TraitValues, 2.5)
            htrait = np.nanpercentile(TraitValues, 97.5)
            print("median trait value: {0:0.2f}, 2.5 - 97.5 percentiles: {1:0.2f} - {2:0.2f}".format(mtrait, ltrait, htrait))
        else:
            catcounts = {cat : len([x for x in TraitValues if x == cat]) for cat in Categories}
            print("category histogram:")
            for cat in Categories:
                print(f"\t{cat}:\t{catcounts[cat]}")
        ## allele statistics
        print("total number of alleles:", np.sum([len(Alleles[locus]) for locus in loci]))
        for locus in sorted(Alleles.keys()):
            print(f"number of {locus} alleles: {len(Alleles[locus])}")
        numCompleteHaplotypes = 0
        for idx in range(len(data)):
            compl = [sum(1 for x in alleleVec if x) == 1
                     for locus in sorted(Alleles.keys())
                     for alleleVec in AlleleVecs[locus][idx]]
            if all(compl):
                ## subject has with 2-field typed haplotype
                numCompleteHaplotypes += 1
        print("number of complete haplotypes: {0}".format(numCompleteHaplotypes))
        ## Covariates
        if len(covariateFieldNames) > 0:
            print("covariate statistics:")
        for cfn in covariateFieldNames:
            mcov = np.nanmean(Covariates[cfn]["Values"])
            lcov = np.nanpercentile(Covariates[cfn]["Values"], 2.5)
            hcov = np.nanpercentile(Covariates[cfn]["Values"], 97.5)
            print(f"\t'{cfn}' mean value: {mcov:0.2f}, 2.5 - 97.5 percentiles: {lcov:0.2f} - {hcov:0.2f}")

    ## return relevant data in a dictionary
    dataDict = {
        "TraitValues" : TraitValues,
        "CensCodes" : CensCodes,
        "AlleleVecs" : AlleleVecs,
        "Alleles" : Alleles,
        "Covariates" : Covariates
    }
    ## include censoring bounds if the trait value is a real number
    if not categorical:
        dataDict.update({"CensBounds" : CensBounds})
    else:
        dataDict.update({"Categories" : Categories})
    return dataDict
