#!/usr/bin/env python

from __future__ import (print_function, division)
from builtins import (map, range)
import requests
import warnings
import urllib.request
import urllib.parse
import lxml.html
import csv
import re
import os

def require_folder(folder_name):
    if not os.path.isdir(folder_name):
        print("new folder ({}) created".format(folder_name))
        try:
            os.makedirs(folder_name)
        except OSError:
            if not os.path.isdir(folder_name): raise


data_folder = "data"
example_folder = os.path.join(data_folder, "example")
HIV_refseq_folder = os.path.join(example_folder, "HIV1_subtype_reference")


kiepiela_url = "https://www.hiv.lanl.gov/content/immunology/hlatem/study3/subjects.csv"
kiepiela_filename = "subjects_kiepiela.csv"

response = urllib.request.urlopen(kiepiela_url)
page = response.read()

require_folder(example_folder)
with open(os.path.join(example_folder, kiepiela_filename), 'w') as f:
    f.write(page.decode("utf-8"))

with open(os.path.join(example_folder, kiepiela_filename), 'r') as f:
    reader = csv.DictReader(f, delimiter=',', quotechar='"')
    table = [row for row in reader]

## TODO: resolve data

## alignments

#alignments_url = "https://www.hiv.lanl.gov/content/sequence/NEWALIGN/align.html"
alignments_url = "https://www.hiv.lanl.gov/cgi-bin/NEWALIGN/align.cgi"
regions = ["ENV", "GAG", "NEF", "POL", "REV", "TAT", "VIF", "VPR", "VPU"]
sequence_dicts = {}
for region in regions:
    payload = {
        "ORGANISM" : "HIV",
        "ALIGN_TYPE" : "REF",
        "SUBORGANISM" : "HIV1",
        "PRE_USER" : "predefined",
        "REGION" : region,
        "BASETYPE" : "PRO",
        "YEAR" : "2010",
        "FORMAT" : "fasta",
        "GENO_SUB" : "A-K",
        "submit" : "Get Alignment"
    }
    response = requests.post(alignments_url, data=payload)
    tree = tree = lxml.html.fromstring(response.text)
    all_pre_content = [elt.text for elt in tree.findall(".//pre")]
    all_fasta_content = [x for x in all_pre_content if len(x) > 0 and x[0] == ">"]
    if len(all_fasta_content) == 0:
        print(response.text)
        raise Exception("unable to find any fasta content on requested LANL page")
    elif len(all_fasta_content) > 1:
        warnings.warn("multiple fasta elements found. using the first one")
    fasta_content = all_fasta_content[0]
    def remove_naa(s): return "".join(ch for ch in s if ch not in "-*")
    sequences = {
        seq.split("\n")[0] : remove_naa("".join(seq.split("\n")[1:]))
        for seq in fasta_content.split(">")[1:]
    }
    sequence_dicts[region] = sequences

refseq_ids = [k for sequences in sequence_dicts.values() for k in sequences.keys()]
refseq_ids = sorted(list(set(refseq_ids)))

fasta_file_dicts = {
    refseq_id : {
        region : sequence_dicts[region][refseq_id]
        for region in regions
    }
    for refseq_id in refseq_ids
}

require_folder(HIV_refseq_folder)
for refseq_id in refseq_ids:
    with open(os.path.join(HIV_refseq_folder, refseq_id + ".fasta"), 'w') as f:
        for region in regions:
            f.write(">" + region + "\n" + fasta_file_dicts[refseq_id][region] + "\n")
