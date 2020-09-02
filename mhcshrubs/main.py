"""@package main
This module contains the main function

The main function reads command line arguments and choses what model to fit
"""
from __future__ import (print_function, division)
from builtins import (map, str, zip, range)
import sys
import os
import json
import argparse
import warnings
import numpy as np ## log transform
from mhcshrubs import (crossval, fittrees, mhcclus, sequences, compare, mhctools, fetcher)
from mhcshrubs import auxiliary as aux
from mhcshrubs import definitions as defn


def main():
    parser = argparse.ArgumentParser(prog="mhcshrubs",
                description="Find HLA associations with disease traits using HLA trees.",
                epilog="""The program makes use of the ETE tree rendering methods,
                and therefore requires an available X-server to run correctly.
                If the program is run on a server or e.g. in a screen environment,
                then one can make use of xvfb:
                `xvfb-run mhcshrubs [options]`,
                or, alternatively:
                `Xvfb :1 &; export DISPLAY=:1; mhcshrubs [options]`
                where 1 can be replaced with other numbers""")
    parser.add_argument("-j", dest="json_file", type=str, default="example.json",
                        help="""the name of a json file containing names of files
                        containing data""")
    parser.add_argument("-y", dest="no_confirm", action="store_true", default=False,
                        help="""don't wait for user confirmation, but automatically choose 'yes' (y)""")
    parser.add_argument("--dry-run", dest="dry_run", action="store_true", default=False,
                        help="""do not run lengthy computations (e.g. JAGS and NetMHCpan)""")
    parser.add_argument("--use-cache", dest="use_cache", action="store_true", default=False,
                        help="""try to use cached output (e.g. from JAGS and NetMHCpan)""")
    parser.add_argument("--model", dest="model", choices=["tree", "pmm", "null"], default="tree",
                        help="""select a model...
                        tree: MHC tree model (default),
                        pmm: phylogenetic mixed model,
                        null: null model""")
    parser.add_argument("--mhc-similarity", dest="mhc_similarity", choices=["pseudoseq", "kir", "binding", "group", "trivial"],
                        default="pseudoseq", help="""compute MHC similarity based on...
                        pseudoseq: MHC pseudo protein sequence (default),
                        kir: amino-acids in the MHC that interact with KIR,
                        binding: NetMHCpan binding predictions,
                        group: cluster at one-field typing,
                        trivial: independent model""")
    parser.add_argument("--aa-similarity", dest="aa_similarity", choices=["pmbec", "blosum50", "blosum62", "manhattan"],
                        default="pmbec", help="""compute distance between amino-acid sequences using...
                        pmbec: the PMBEC matrix (default),
                        blosum50, blosum62: a BLOSUM matrix,
                        manhattan: the Manhattan distance""") ## TODO
    parser.add_argument("--prior", dest="prior", choices=["norm", "dexp"], default="norm",
                        help="""prior distribution used for branch weights...
                        norm: the Normal distribution (default),
                        dexp: double exponential distribution""")
    parser.add_argument("--multiplicity", dest="multiplicity", choices=["single", "double"], default="double",
                        help="""multiplicity of homozygous branch weights...
                        single: the 'generalized homozygosity' model,
                        double: count double weights twice (default)""")
    parser.add_argument("--hetr-adv", dest="hetr_adv", type=str, default=None, choices=["tree", "discrete"],
                        help="""add a parameter representing the heterozygote advantage.
                        discrete: The level of heterozygosity is calculated by allele identity,
                        tree: use the MHC tree to compute 'fractional' heterozygosity.""")
    parser.add_argument("--chain-len", dest="chain_len", type=int, default=10000,
                        help="""the length of the MCMC before thinning (default: 10000).
                        warmup is of equal length.""")
    parser.add_argument("--chain-thin", dest="chain_thin", type=int, default=10,
                        help="""amount of thinning of the MCMC. decreases auto-correlation (default: 10)""")
    parser.add_argument("--num-chains", dest="num_chains", type=int, default=4,
                        help="""number of independent MCMC runs (needs to be > 1 for valid Rhat)""")
    parser.add_argument("--name-base", dest="name_base", type=str, default=None,
                        help="""optional string used to form output file names.
                        if not specified, a name base is generated automatically.""")
    parser.add_argument("--sampler", dest="sampler", choices=["jags", "stan"], default="jags",
                        help="""MCMC sampler...
                        jags: use JAGS (Gibbs sampling),
                        stan: use Stan (Hamiltonian Monte-Carlo)""")
    parser.add_argument("--cross-val", dest="cross_val", choices=["alleles"], default=None,
                        help="""specify leave-one-out cross-validation (optinal)...
                        alleles: leave-out alleles by masking the response variable of their owners""")
    parser.add_argument("--no-concurrency", dest="no_concurrency", action="store_true", default=False,
                        help="""disable all parallel computations (e.g. for debugging)""")
    parser.add_argument("--compare", nargs=2, dest="models", default=None, help="""compare two models""")
    args = parser.parse_args()

    if args.models is not None:
        compare.compareModelsSEM(*args.models)
        return 0

    ## show the chosen options
    print("\nAbout to run MHCshrubs with the following options:\n")
    max_key_len = max(map(len, args.__dict__.keys()))
    for key in sorted(args.__dict__.keys()):
        val = args.__dict__[key]
        if val is not None:
            print("{0:{1}s}:".format(key, max_key_len+1), val)

    with open(args.json_file) as fileHandle:
        parDict = json.loads(fileHandle.read())

    ## make a name-base for the output
    if args.name_base is None:
        name_base = ""
        if "id" in parDict.keys():
            name_base += parDict["id"]
        name_base += "_" + args.model
        if args.model != "null":
            if args.model == "tree":
                name_base += "_" + args.prior + "_" + args.multiplicity
            name_base += "_" + args.mhc_similarity
            if args.mhc_similarity in ["pseudoseq", "kir"]:
                name_base += "_" + args.aa_similarity
            if args.hetr_adv is not None:
                name_base += "_" + args.hetr_adv + "adv"
        cl = str(args.chain_len)
        cl = "{0}e{1}".format(cl[0], len(cl)-1)
        name_base += "_" + cl
    else:
        name_base = args.name_base
    print("using model name:", name_base)

    ## check that an X-server is available for ETE rendering functions
    if not aux.isXServerAvailable():
        warnings.warn("no X-server avaliable. Run the program with argument '--help' for advice.")

    ## the user can check the options, and decide if (s)he wants to proceed
    while not args.no_confirm:
        answer = input("\nDo you want to continue? [Y/n]")
        if answer in ['y', 'Y', '']:
            break ## continue by breaking the while loop
        elif answer in ['n', 'N']:
            return 0 ## abort

    ## input files
    subjectFileName = parDict["subjectFileName"]
    hlaFileName = parDict["alleleFreqFileName"] if "alleleFreqFileName" in parDict.keys() else None
    pSeqFileName = parDict["pSeqFileName"] if "pSeqFileName" in parDict.keys() else None
    fastaFileName = parDict["fastaFileName"] if "fastaFileName" in parDict.keys() else None

    ## select a file with an amino-acid similarity matrix
    if args.aa_similarity == "pmbec":
        aaCovFileName = os.path.join("resources", "PMBEC.MAT")
    elif args.aa_similarity == "blosum62":
        aaCovFileName = os.path.join("resources", "BLOSUM62")
    elif args.aa_similarity == "blosum50":
        aaCovFileName = os.path.join("resources", "BLOSUM50")
    elif args.aa_similarity == "manhattan":
        raise Exception("Manhattan distance is not implemented yet") ## TODO: use a different system
    else:
        raise Exception("invalid aa-similarity option")

    ## output files
    if args.cross_val is None:
        summaryFileName = "data/summary.{}.tsv".format(name_base)
    else:
        summaryFileName = "data/cv_{0}_summary.{1}.tsv".format(args.cross_val, name_base)

    ## add the right directory to the file names
    work_folder = os.getcwd()

    ## check that the output directory exists, and make if not dry-run
    if "outputFolder" in parDict.keys():
        outputFolder = os.path.join(work_folder, parDict["outputFolder"])
    else:
        outputFolder = os.path.join(work_folder, "data")
    figureFolder = os.path.join(outputFolder, "figures")
    if not os.path.isdir(figureFolder):
        if args.dry_run:
            print("folder for writing data and/or figures does not exist")
        else:
            print("new folder for writing data and/or figures created")
    try:
        os.makedirs(figureFolder)
    except OSError:
        if not os.path.isdir(figureFolder): raise

    print("\nworking from folder:", work_folder)

    subjectFileName = os.path.join(work_folder, subjectFileName)
    if hlaFileName is not None:
        hlaFileName = os.path.join(work_folder, hlaFileName)
    if pSeqFileName is not None:
        pSeqFileName = os.path.join(work_folder, pSeqFileName)
    if fastaFileName is not None:
        fastaFileName = os.path.join(work_folder, fastaFileName)
    aaCovFileName = os.path.join(defn.ROOT_DIR, aaCovFileName)
    summaryFileName = os.path.join(work_folder, summaryFileName)

    ## by default: pass parallel=True to functions
    parallel = False if args.no_concurrency else True

    if args.model == "null":
        result = funFitNull(subjectFileName, summaryFileName, parDict,
                            dry_run=args.dry_run, chain_len=args.chain_len,
                            chain_thin=args.chain_thin, name_base=name_base)
    elif args.model == "tree":
        kwargs = {
            "dry_run"       : args.dry_run,
            "use_cache"     : args.use_cache,
            "chain_len"     : args.chain_len,
            "chain_thin"    : args.chain_thin,
            "num_chains"    : args.num_chains,
            "name_base"     : name_base,
            "prior"         : args.prior,
            "multiplicity"  : args.multiplicity,
            "hetr_adv"      : args.hetr_adv,
            "sampler"       : args.sampler,
            "cross_val"     : args.cross_val,
            "parallel"      : parallel
        }
        if args.mhc_similarity == "trivial":
            result = funFitNotree(subjectFileName, hlaFileName, summaryFileName,
                parDict, **kwargs)
        elif args.mhc_similarity == "group":
            result = funFitGroup(subjectFileName, hlaFileName, summaryFileName,
                parDict, **kwargs)
        elif args.mhc_similarity == "pseudoseq":
            if pSeqFileName is None:
                raise Exception("no field 'pSeqFileName' found in json file, which is required for option 'pseudoseq'")
            result = funFitTree(subjectFileName, hlaFileName, pSeqFileName,
                       aaCovFileName, summaryFileName, parDict, **kwargs)
        elif args.mhc_similarity == "kir":
            result = funFitKirTree(subjectFileName, hlaFileName, aaCovFileName,
                summaryFileName, parDict, **kwargs)
        elif args.mhc_similarity == "binding":
            if fastaFileName is None:
                raise Exception("no field 'fastaFileName' found in json file, which is required for option 'binding'")
            result = funFitBindingTree(subjectFileName, hlaFileName, fastaFileName,
                summaryFileName, parDict, **kwargs)
        else:
            raise Exception("MHC similarity '{}' not implemented".format(args.mhc_similarity))
    elif args.model == "pmm":
        result = funFitPMM(subjectFileName, hlaFileName, pSeqFileName,
                           aaCovFileName, summaryFileName, parDict,
                           dry_run=args.dry_run, use_cache=args.use_cache,
                           chain_len=args.chain_len, chain_thin=args.chain_thin,
                           name_base=name_base, sampler=args.sampler,
                           cross_val=args.cross_val, parallel=parallel,
                           num_chains=args.num_chains)
    else:
        raise Exception("model '{}' not implemented".format(args.model))
    return 0


def funFitTree(subjectFileName, hlaFileName, pSeqFileName, aaCovFileName, summaryFileName,
               parDict, chain_len=1000, chain_thin=10, num_chains=4, dry_run=False,
               use_cache=False, name_base="anon", prior="norm", multiplicity="double", hetr_adv=None,
               sampler="jags", cross_val=None, parallel=True):
    ## get data
    traitFieldName = parDict["traitFieldName"]
    ## determine the type of trait (continuous or categorical)
    traitType = parDict["traitType"]
    if traitType == "categorical":
        categorical = True
    elif traitType == "continuous":
        categorical = False
    else:
        raise Exception(f"invalid traitType '{traitType}' in json file")
    traitTransform = (lambda x: x) if categorical else np.log10 ## FIXME!!
    alleleFieldNames = parDict["alleleFieldNames"]
    dataDict = fetcher.importSubjectData(subjectFileName, traitFieldName, alleleFieldNames,
                                         verbose=True, traitTransform=traitTransform,
                                         categorical=categorical)
    hlaAlleles = dataDict["Alleles"]
    ## define models by their HLA tree
    pSeqDict = mhcclus.mkPseqDict(pSeqFileName)
    ## restrict the pSeqDict to alleles in the dataset
    loci = sorted(hlaAlleles.keys())
    pSeqDictRestr = dict((hla, pSeqDict[hla]) for X in loci for hla in hlaAlleles[X] if hla.suffix != 'N')
    ## choose distances between amino-acids
    aaDistDict = sequences.mkAaDistDict(aaCovFileName, rescale=True, ignore_aas=sequences.ambiguous_aas)
    ## use the data to make a tree (in Newick format)
    ## using representatives for NetMHCpan equivalence classes
    newick_str_pseq, represDict = mhcclus.mkNewickTree(pSeqDictRestr, aaDistDict=aaDistDict, collapse=True)

    ## run the model
    if cross_val is None: ## simple case: fit one model
        if sampler == "jags":
            if categorical:
                result = fittrees.fitTreeWeightsCat(dataDict, parDict, hlaFileName, newick_str_pseq,
                                                    represDict=represDict, chain_len=chain_len,
                                                    chain_thin=chain_thin, num_chains=num_chains, modelName=name_base,
                                                    verbose=True, dry_run=dry_run, use_cache=use_cache,
                                                    prior=prior, multiplicity=multiplicity, hetr_adv=hetr_adv, parallel=parallel)
            else:
                result = fittrees.fitTreeWeights(dataDict, parDict, hlaFileName, newick_str_pseq,
                                                 represDict=represDict, chain_len=chain_len,
                                                 chain_thin=chain_thin, num_chains=num_chains, modelName=name_base,
                                                 verbose=True, dry_run=dry_run, use_cache=use_cache,
                                                 prior=prior, multiplicity=multiplicity, hetr_adv=hetr_adv, parallel=parallel)
        elif sampler == "stan":
            if categorical:
                raise Exception("categorical Stan sampler not implemented")
            else:
                result = fittrees.fitTreeWeightsStan(dataDict, parDict, hlaFileName, newick_str_pseq,
                                                     represDict=represDict, chain_len=chain_len,
                                                     chain_thin=chain_thin, num_chains=num_chains, modelName=name_base,
                                                     verbose=True, dry_run=dry_run, use_cache=use_cache,
                                                     prior=prior, multiplicity=multiplicity, hetr_adv=hetr_adv, wbic_sampling=False)
        else:
            raise Exception(f"invalid sampler '{sampler}' given")
    elif cross_val == "alleles": ## use loo-allele cross-validation
        ## TODO: categorical data
        result = crossval.crossValidateAlleles(dataDict, parDict, hlaFileName, newick_str_pseq,
                                               summaryFileName, represDict=represDict,
                                               chain_len=chain_len, chain_thin=chain_thin,
                                               prior=prior, modelName=name_base,
                                               multiplicity=multiplicity, hetr_adv=hetr_adv,
                                               dry_run=dry_run, use_cache=use_cache,
                                               parallel=parallel, verbose=True)
    else:
        raise Exception("invalid cross-validation mode '{}' given".format(cross_val))
    return result


def funFitKirTree(subjectFileName, hlaFileName, aaCovFileName, summaryFileName, parDict,
                  chain_len=1000, chain_thin=10, num_chains=4, dry_run=False, use_cache=False,
                  name_base="anon", prior="norm", multiplicity="double", hetr_adv=None,
                  sampler="jags", cross_val=None, parallel=True):
    traitFieldName = parDict["traitFieldName"]
    alleleFieldNames = parDict["alleleFieldNames"]
    dataDict = fetcher.importSubjectData(subjectFileName, traitFieldName, alleleFieldNames,
                                         verbose=True, traitTransform=np.log10)
    hlaAlleles = dataDict["Alleles"]

    kir_positions = [77, 80, 81, 82, 83] ## TODO load file (e.g. Martin et al, Nature (2002))
    start = 23 ## same as NetMHCpan psuedo sequences. TODO: verify: check the motif

    ## choose distances between amino-acids
    aaDistDict = sequences.mkAaDistDict(aaCovFileName)

    newick_str_kir, represDict = mhcclus.mkSubSeqNewickTree(hlaAlleles, kir_positions, start=start,
                                                            aaDistDict=aaDistDict, verbose=True)

    ## run the model
    result = fittrees.fitTreeWeights(dataDict, parDict, hlaFileName, newick_str_kir,
                                     represDict=represDict, chain_len=chain_len, chain_thin=chain_thin,
                                     num_chains=num_chains, modelName=name_base,
                                     verbose=True, dry_run=dry_run, use_cache=use_cache,
                                     prior=prior, multiplicity=multiplicity, hetr_adv=hetr_adv,
                                     parallel=parallel)
    return result


def funFitBindingTree(subjectFileName, hlaFileName, fastaFileName, summaryFileName, parDict,
                      chain_len=1000, chain_thin=10, num_chains=4, dry_run=False, use_cache=False,
                      name_base="anon", prior="norm", multiplicity="double", hetr_adv=None,
                      sampler="jags", cross_val=None, parallel=True):
    traitFieldName = parDict["traitFieldName"]
    alleleFieldNames = parDict["alleleFieldNames"]
    dataDict = fetcher.importSubjectData(subjectFileName, traitFieldName, alleleFieldNames,
                                         verbose=True, traitTransform=np.log10)
    hlaAlleles = dataDict["Alleles"]
    ## choose a working folder for netMHC
    work_folder = os.getcwd()
    cache_dir = os.path.join(work_folder, "netmhc-cache")
    ## define models by their HLA tree
    newick_str_bind, represDict = \
        mhcclus.mkPeptideBindingNewickTree(fastaFileName, hlaAlleles, collapse=True,
                                           verbose=True, dry_run=dry_run, parallel=parallel,
                                           use_cache=use_cache, cache_dir=cache_dir)
    ## run the model
    result = fittrees.fitTreeWeights(dataDict, parDict, hlaFileName, newick_str_bind,
                                     represDict=represDict, chain_len=chain_len,
                                     chain_thin=chain_thin, num_chains=num_chains,
                                     modelName=name_base, verbose=True, dry_run=dry_run,
                                     use_cache=use_cache, prior=prior,
                                     multiplicity=multiplicity, hetr_adv=hetr_adv,
                                     parallel=parallel)
    return result


def funFitNotree(subjectFileName, hlaFileName, summaryFileName, parDict,
                 chain_len=1000, chain_thin=10, num_chains=4, dry_run=False,
                 use_cache=False, name_base="anon", prior="norm",
                 multiplicity="double", hetr_adv=None, sampler="jags",
                 cross_val=None, parallel=True):
    traitFieldName = parDict["traitFieldName"]
    alleleFieldNames = parDict["alleleFieldNames"]
    dataDict = fetcher.importSubjectData(subjectFileName, traitFieldName, alleleFieldNames,
                                         verbose=True, traitTransform=np.log10)
    hlaAlleles = dataDict["Alleles"]
    ## make a trivial newick string
    newick_str_notree = mhcclus.mkTrivialNewickTree(hlaAlleles)

    ## run the model
    if cross_val is None:
        if sampler == "jags":
            result = fittrees.fitTreeWeights(dataDict, parDict, hlaFileName,
                                             newick_str_notree, chain_len=chain_len,
                                             chain_thin=chain_thin, num_chains=num_chains,
                                             modelName=name_base, verbose=True,
                                             dry_run=dry_run, use_cache=use_cache,
                                             prior=prior,
                                             multiplicity=multiplicity, hetr_adv=hetr_adv,
                                             parallel=parallel)
        elif sampler == "stan":
            result = fittrees.fitTreeWeightsStan(dataDict, parDict, hlaFileName,
                                                 newick_str_notree, chain_len=chain_len,
                                                 chain_thin=chain_thin, num_chains=num_chains,
                                                 modelName=name_base, verbose=True,
                                                 dry_run=dry_run, prior=prior,
                                                 multiplicity=multiplicity, hetr_adv=hetr_adv,
                                                 wbic_sampling=False)
        else:
            raise Exception("invalid sampler '{}' given".format(sampler))
    elif cross_val == "alleles":
        result = crossval.crossValidateAlleles(dataDict, parDict, hlaFileName,
                                               newick_str_notree, summaryFileName,
                                               chain_len=chain_len, chain_thin=chain_thin,
                                               prior=prior, modelName=name_base,
                                               multiplicity=multiplicity, hetr_adv=hetr_adv,
                                               dry_run=dry_run, use_cache=use_cache,
                                               parallel=parallel, verbose=True)
    else:
        raise Exception("invalid cross-validation mode '{}' given".format(cross_val))
    return result


def funFitGroup(subjectFileName, hlaFileName, summaryFileName, parDict,
                chain_len=1000, chain_thin=10, num_chains=4, dry_run=False,
                use_cache=False, name_base="anon", prior="norm",
                multiplicity="double", hetr_adv=None, sampler="jags",
                cross_val=None, parallel=True):
    traitFieldName = parDict["traitFieldName"]
    alleleFieldNames = parDict["alleleFieldNames"]
    dataDict = fetcher.importSubjectData(subjectFileName, traitFieldName, alleleFieldNames,
                                         verbose=True, traitTransform=np.log10)
    hlaAlleles = dataDict["Alleles"]
    ## make a trivial newick string
    newick_str_group = mhcclus.mkGroupNewickTree(hlaAlleles)
    ## run the model
    result = fittrees.fitTreeWeights(dataDict, parDict, hlaFileName, newick_str_group,
                                     chain_len=chain_len, chain_thin=chain_thin, num_chains=num_chains,
                                     modelName=name_base, verbose=True,
                                     dry_run=dry_run, use_cache=use_cache, prior=prior,
                                     multiplicity=multiplicity, hetr_adv=hetr_adv,
                                     parallel=parallel)
    return result


def funFitPMM(subjectFileName, hlaFileName, pSeqFileName,
              aaCovFileName, summaryFileName, parDict,
              chain_len=1000, chain_thin=10, num_chains=4,
              dry_run=False, use_cache=False,
              name_base="anon", sampler="jags", cross_val=None, parallel=True):
    """Todo: implement cross_val, PMM could also be used for many other situations..."""
    traitFieldName = parDict["traitFieldName"]
    alleleFieldNames = parDict["alleleFieldNames"]
    dataDict = fetcher.importSubjectData(subjectFileName, traitFieldName, alleleFieldNames,
                                         verbose=True, traitTransform=np.log10)
    hlaAlleles = dataDict["Alleles"]
    ## define models by their HLA tree
    pSeqDict = mhcclus.mkPseqDict(pSeqFileName)
    ## restrict the pSeqDict to alleles in the dataset
    loci = sorted(hlaAlleles.keys())
    pSeqDictRestr = dict((hla, pSeqDict[hla]) for X in loci for hla in hlaAlleles[X])
    ## choose distances between amino-acids
    aaDistDict = sequences.mkAaDistDict(aaCovFileName, rescale=True, ignore_aas=sequences.ambiguous_aas)
    ## use the data to make a tree (in Newick format)
    ## using representatives for NetMHCpan equivalence classes
    newick_str_pseq, represDict = mhcclus.mkNewickTree(pSeqDictRestr, aaDistDict=aaDistDict, collapse=True)
    hlaCovMat, alleleList = mhcclus.treeToCovmat(newick_str_pseq)
    hlaCovMatHeader = [mhctools.MhcObject(allele, fmt="Newick") for allele in alleleList]
    ## run the model
    result = fittrees.fitPMMweights(dataDict, parDict, hlaFileName, hlaCovMat, hlaCovMatHeader,
                                    represDict=represDict, chain_len=chain_len,
                                    chain_thin=chain_thin, num_chains=num_chains,
                                    modelName=name_base, verbose=True, dry_run=dry_run,
                                    use_cache=use_cache, parallel=parallel)
    return result


def funFitNull(subjectFileName, summaryFileName, parDict,
               chain_len=1000, chain_thin=10, num_chains=4,
               dry_run=False, use_cache=False, name_base="anon", sampler="jags",
               parallel=True):
    traitFieldName = parDict["traitFieldName"]
    alleleFieldNames = parDict["alleleFieldNames"]
    dataDict = fetcher.importSubjectData(subjectFileName, traitFieldName, alleleFieldNames,
                                         verbose=True, traitTransform=np.log10)
    ## the null model does not require any HLA data
    result = fittrees.fitNullModel(dataDict, parDict, chain_len=chain_len, chain_thin=chain_thin,
                                   num_chains=num_chains, modelName=name_base,
                                   verbose=True, dry_run=dry_run,
                                   use_cache=use_cache, parallel=parallel)
    return result
