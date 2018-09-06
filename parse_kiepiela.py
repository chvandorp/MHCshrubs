# Import HIV VL and HLA data
# 1. take care of censoring
# 2. find admissible HLA alleles

import csv
import numpy as np
from mhcshrubs import auxiliary as aux
from mhcshrubs import mhctools

data_folder = "data"
example_folder = os.path.join(data_folder, "example")

linesep = "-----------\n"

## import Kiepiela data
with open(example_folder + "/subjects_kiepiela.csv", 'r') as f:
    reader = csv.DictReader(f, delimiter=',', quotechar='"')
    records = [r for r in reader]
    variables = sorted(records[0].keys())

## Print keys
print("keys: \n" + linesep +"\n".join(variables))

subjectIDKey = "Patient_ID"
subjectIDs = aux.unique(r[subjectIDKey] for r in records)
print(linesep + "number of patients:{}".format(len(subjectIDs)))
print(linesep + "number of data points: {}".format(len(records)))

## Sort by patient
def collapseDicts(dicts):
    keys = aux.unique(aux.flatten([d.keys() for d in dicts]))
    cdict = dict((key, [d[key] if key in d.keys() else None for d in dicts]) for key in keys)
    cdict["Number_Observations"] = len(dicts)
    return cdict

recordsByID = {
    ID : collapseDicts(filter(lambda subjd: subjd[subjectIDKey] == ID, records)))
    for ID in subjectIDs
}

## Print an example
exampleID = np.random.choice(subjectIDs)
exampleRecord = recordsByID[exampleID]
print(linesep)
for key in variables:
    print(key + "\t" + "\t".join(exampleRecord[key]))


def mustGetUnique(values):
    """
    returns a unique value from values
    /!\ Raise an exception when different samples yield different values.
    """
    unique_vals = aux.unique(values)
    if len(unique_vals) != 1:
        raise Exception("different samples yield different HLA alleles.")
    return values[0]

## resolve VL sensoring and multiple values

missing_code = "missing"
left_censored_code = "left_censored"
right_censored_code = "right_censored"
uncensored_code = "uncensored"

def getVLcensor(vlstring):
    """
    given a string representing VL, a censoring code and value is returned
    """
    if vlstring == '-':
        return (np.nan, missing_code) ## TODO: parameter
    elif vlstring[0] == '>': ## the value is a lower bound: right censoring
        vlbound = float(vlstring[1:])
        return (vlbound, right_censored_code)
    elif vlstring[0] == '<': ## the value is an upper bound: left censoring
        vlbound = float(vlstring[1:])
        return (vlbound, left_censored_code)
    else: ## plain VL value
        vl = float(vlstring)
        return (vl, uncensored_code)


def getVL(vlstrings):
    """
    get a representative VL from a list of (possibly) different VLs.
    There may be missing values and censoring.
    Return a censoring code.
    The VL is computed as follows:
    1. determine for each value the type of censoring
    2. get all non-censored values and return the average
    3. if there are no non-censored values, take the average of the right or left censored values
    /!\ if both left and right censored values exist, raise an exception
    4. if there are only missing values, return a NaN
    /!\ if there are no valid censor codes, ranse an exception
    """
    censors = [getVLcensor(vlstr) for vlstr in vlstrings]
    uncensored_vals = [x for x, c in censors if c==uncensored_code]
    left_censored_vals = [x for x, c in censors if c==left_censored_code]
    right_censored_vals = [x for x, c in censors if c==right_censored_code]
    missing_vals = [x for x, c in censors if c==missing_code]
    if len(uncensored_vals) > 0: ## simplest case: return the mean of all uncensored values
        return (np.mean(uncensored_vals), uncensored_code)
    elif len(right_censored_vals) > 0: ## more difficult: no uncensored values exist: get left or right censored values
        if len(left_censored_vals) > 0:
            raise Exception("can not handle both left and right censored values. FIXME.")
        return (np.mean(right_censored_vals), right_censored_code)
    elif len(left_censored_vals) > 0:
        return (np.mean(left_censored_vals), left_censored_code)
    elif len(missing_vals) > 0:
        return (np.nan, missing_code)
    else:
        raise Exception("no valid censoring codes.")

## resolve HLA censoring and multiple values

missing_code = "missing"
one_field_code = "one_field"
uncensored_code = "uncensored"

manualCorrectionDict = {
    "C06var" : "C06",
    "B8108" : "B81",
    "A69" : "A6901"
}

## NOTES: the only member of the serotype A69 is A6901. unfortunately, A6901 is not common in the SSA population sample

def getHLAcensor(hlastring):
    """
    Determine the kind of censoring of the HLA alleles.
    perform some small manual corrections that would otherwise be hard to automize
    """
    ## first pass: manual correction
    if hlastring in manualCorrectionDict.keys():
        hlastring = manualCorrectionDict[hlastring]
    ## check censoring
    if hlastring == '-': ## missing
        return (hlastring, missing_code)
    if len(hlastring) == 3: ## one field (or serotype)
        return (hlastring, one_field_code)
    elif len(hlastring) == 5: ## uncensored
        return (hlastring, uncensored_code)
    else:
        raise Exception("invalid HLA string.")

def parseHLAstring(hlastring, censor_code, locus=None):
    ## parse the HLA based on the censor code
    if censor_code == missing_code:
        return mhctools.MhcObject(locus=locus)
    loc = hlastring[0]
    if locus is not None and loc != locus:
        raise Exception("the passed locus does not correspond to locus inferred from the HLA string")
    elif censor_code == one_field_code:
        field1 = hlastring[1:3]
        return mhctools.MhcObject(locus=loc, field1=field1)
    elif censor_code == uncensored_code:
        field1 = hlastring[1:3]
        field2 = hlastring[3:5]
        return mhctools.MhcObject(locus=loc, field1=field1, field2=field2)
    else:
        raise Exception("invalid censor code.")


### find admissible HLA alleles for censored HLA values

## make a list of all 2-field alleles in the dataset

hlaLocusKeyDict = {
    'A' : ["HLA_A1", "HLA_A2"],
    'B' : ["HLA_B1", "HLA_B2"],
    'C' : ["HLA_Cw1", "HLA_Cw2"]
}

cohortHlaAlleles = {}
for X in "ABC":
    uniques = [mustGetUnique(recordsByID[ID][key]) for key in hlaLocusKeyDict[X] for ID in subjectIDs]
    censors = [getHLAcensor(hlastring) for hlastring in uniques]
    ## only use uncensored values
    alleles = [parseHLAstring(x, c, X) for x, c in censors if c == uncensored_code]
    cohortHlaAlleles[X] = aux.unique(alleles)

## import HLA alleles from another cohort

with open(rootFolder + "hla/mhc-top0.99counts-ncbi-Sub-Saharan-Africa.tsv", 'r') as f:
    hla_table = [row.split('\t') for row in f.read().split('\n') if row != '']
    hla_header = hla_table[0]
    hla_table = hla_table[1:]
hlaCountDict = dict((mhctools.MhcObject(row[0]), int(row[1])) for row in hla_table)

popHlaAlleles = {
    X : sorted(filter(lambda x: x.locus == X, hlaCountDict.keys()),
               key=lambda x: hlaCountDict[x], reverse=True)
    for X in "ABC"
}

## make a combined list of admissible alleles
admHlaAlleles = {
    X : aux.unique(cohortHlaAlleles[X] + popHlaAlleles[X])
    for X in "ABC"
}

## import alleles belonging to serotypes
with open(os.path.join(data_folder, "HLA_serotypes.tsv"), 'r') as f:
    serotype_table = [row.split('\t') for row in f.read().split('\n') if row !='']
serotype_header = serotype_table[0]
serotype_table = serotype_table[1:]
serotypes = aux.unique([row[1] for row in serotype_table])
serotypeDict = {
    ser : [mhctools.MhcObject(row[0]) for row in serotype_table if row[1] == ser]
    for ser in serotypes
}

def getSerotypeAdmissibleAlleles(hlastring, locus):
    hlastring, censor_code = getHLAcensor(hlastring)
    ## assume serotyping for censored HLAs
    members = []
    if censor_code == one_field_code and hlastring in serotypeDict.keys():
        members = serotypeDict[hlastring]
    ## filter members that are NOT in admAlleles
    admMembers = filter(lambda x: x in admHlaAlleles[locus], members)
    return (admMembers, censor_code)

## FIXME: if a serotype is present in the cohort, but not in the population sample, this will not work.
## see e.g. A69, which was added to the 'manual corrections' for this reason


def getGroupAdmissibleAlleles(hlastring, locus):
    hlastring, censor_code = getHLAcensor(hlastring)
    ## assume single field typing for censored HLAs
    hla = parseHLAstring(hlastring, censor_code, locus)
    admHlas = [x for x in admHlaAlleles[hla.locus] if x in hla]
    return (admHlas, censor_code)

def getAllAdmissibleAlleles(hlastring, locus):
    grpAdm, censor_code1 = getGroupAdmissibleAlleles(hlastring, locus)
    serAdm, censor_code2 = getSerotypeAdmissibleAlleles(hlastring, locus)
    union = aux.unique(grpAdm + serAdm)
    return (union, censor_code1)


# resolve HLA alleles in the subject disctionaries and make a new tsv file.
# apply the following filters:
# 1. remove patients with missing VL data
# 2. remove patients with ONLY missing HLA data

## make resolved dictionaries
resRecordsByID = {}
for ID in subjectIDs:
    record = recordsByID[ID]
    ## make a new record
    record_resolved = {}
    ## resolve missing values
    VL_val, VL_censcode = getVL(record["VL"])
    sex = mustGetUnique(record["Sex"])
    for X in "ABC":
        for key in hlaLocusKeyDict[X]:
            hlastring = mustGetUnique(record[key])
            admAlleles, hla_censor_code = getAllAdmissibleAlleles(hlastring, X)
            hlastring_resolved = ";".join(hla.subtype_str() for hla in admAlleles)
            record_resolved[key] = hlastring_resolved
            record_resolved[key + "_censoring"] = hla_censor_code
    record_resolved["Patient_ID"] = ID
    record_resolved["VL"] = str(VL_val)
    record_resolved["VL_censoring"] = VL_censcode
    record_resolved["Sex"] = sex
    ## add to resolved records
    resRecordsByID[ID] = record_resolved

print("number of subjects:", len(subjectIDs))

## filter out missing VL and completely missing HLA
filteredSubjectIDs = []
for ID in subjectIDs:
    record = resRecordsByID[ID]
    has_VL = record["VL_censoring"] != missing_code
    has_HLA = False
    for X in "ABC":
        for key in hlaLocusKeyDict[X]:
            has_HLA = has_HLA or record[key + "_censoring"] != missing_code
    if has_VL and has_HLA:
        filteredSubjectIDs += [ID]

print("number of filtered subjects:", len(filteredSubjectIDs))


fileName = rootFolder + "durban/subjects_resolved_top0.99_all.tsv"
fieldnames = sorted(resRecordsByID[filteredSubjectIDs[0]])
with open(fileName, 'w') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    for ID in filteredSubjectIDs:
        writer.writerow(resRecordsByID[ID])
