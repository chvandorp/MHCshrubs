#!/usr/bin/env python3
"""
Validate the models by simulating data.
If the model is run as a script, the function
validateModel is called.
"""
from __future__ import (print_function, division)
from builtins import (map, str, zip, range)
import numpy as np
import scipy.stats as sts
import sys
import matplotlib.pyplot as plt
from mhcshrubs import (stantools, statistics)
from mhcshrubs import auxiliary as aux


def simulateData(Ploidy=2, NumLoci=3, MeanNumAlleles=10, NumSubjects=100,
                 AddlSampleSize=1000, MaskingProb=0.1, CensorProb=0.1,
                 sigmaTraitValue=1.0, WBIC=0, NodeWeightPrior=stantools.normal_prior_code):
  NumAlleles = [aux.posPoisson(MeanNumAlleles) for _ in range(NumLoci)]
  AlleleFreqs = [sts.dirichlet.rvs(np.ones(n))[0] for n in NumAlleles]
  AddlAlleles = [np.random.choice(range(1, n+1), p=p, size=AddlSampleSize)
                 for n, p in zip(NumAlleles, AlleleFreqs)]
  AddlAlleleData = [[len([b for b in sample if b==a]) for a in range(1, n+1)]
			   for n, sample in zip(NumAlleles, AddlAlleles)]
  AlleleGroups = [aux.randomGroupMasking(n) for n in NumAlleles]
  ## sample subject alleles
  Alleles = [[[np.random.choice(range(1, n+1), p=p) for n, p in zip(NumAlleles, AlleleFreqs)]
              for _ in range(Ploidy)] for s in range(NumSubjects)]
  ## mask alleles using AlleleGroups ## TODO: completely missing alleles
  AdmAlleles = [[[x[a] if sts.bernoulli.rvs(MaskingProb)==1 else [a]
                  for a, x in zip(Alleles[s][p], AlleleGroups)]
                 for p in range(Ploidy)] for s in range(NumSubjects)]
  NumAdmAlleles = [[[len(x) for x in xs] for xs in xss] for xss in AdmAlleles]
  ## simulate weights
  sigmaAlleleWeight = 1.0 ## TODO: parameter
  AlleleWeights = [sts.norm.rvs(scale=sigmaAlleleWeight, size=n) for n in NumAlleles]
  ## create trait values
  TraitValue = [np.sum([np.sum([w[a-1] for a, w in zip(Alleles[s][p], AlleleWeights)])
                         for p in range(Ploidy)]) + sts.norm.rvs(scale=sigmaTraitValue)
                 for s in range(NumSubjects)]
  ## compute mean and sd of TraitValue in order to sample censoring limits
  meanTrait = np.mean(TraitValue)
  sdTrait = np.std(TraitValue)
  ## censor some TraitValues
  cp = [1.0-CensorProb, CensorProb/3, CensorProb/3, CensorProb/3]
  cc = [0, 1, 2, 3] ## uncensored, left, right, missing
  TraitCensorType = [np.random.choice(cc, p=cp) for _ in TraitValue]
  ## sample censoring limits
  for i, c in enumerate(TraitCensorType):
    if c == 1:
      TraitValue[i] = aux.truncNorm(meanTrait, sdTrait, left=TraitValue[i])
    elif c == 2:
      TraitValue[i] = aux.truncNorm(meanTrait, sdTrait, right=TraitValue[i])
  ## create a data dictionary that can be passed to Stan.
  data = {
    "NumLoci" : NumLoci,
    "NumAlleles" : NumAlleles,
    "AddlAlleleData" : aux.flatten_recursive(AddlAlleleData),
    "Ploidy" : Ploidy,
    "NumSubjects" : NumSubjects,
    "TraitValue" : TraitValue,
    "TraitCensorType" : TraitCensorType,
    "NumAdmAlleles" : NumAdmAlleles,
    "AdmAlleles" : aux.flatten_recursive(AdmAlleles),
    "NumNodes" : np.sum(NumAlleles),
    "LengthEdgeToParent" : np.ones(np.sum(NumAlleles)),
    "TreeMatrix" : np.eye(np.sum(NumAlleles)),
    "WBIC" : WBIC,
    "NodeWeightPrior" : NodeWeightPrior
  }
  ## create a parameter dictionary for diagnostic purposes
  parameters = {
    "AlleleFreqs" : AlleleFreqs,
    "AlleleWeights" : AlleleWeights,
    "Alleles" : Alleles,
    "sigmaTraitValue" : sigmaTraitValue,
    "sigmaAlleleWeight" : sigmaAlleleWeight
  }
  return (data, parameters)

def validateModelFit(fit, parameters, data):
  """
  compare a model fit with known (simulated) parameters
  input:
  - fit -- the result of pystan.StanModel.sampling
  - parameters -- the second return value of simulateData
  output: None
  """
  ## compare allele frequencies
  parameter_samples = fit.extract()
  allele_freq_samples = parameter_samples["alleleFreqs"]
  allele_freqs = parameters["AlleleFreqs"]
  Ploidy = data["Ploidy"]
  NumLoci = data["NumLoci"]
  NumAlleles = data["NumAlleles"]
  NumSubjects = data["NumSubjects"]
  NumAdmAlleles = data["NumAdmAlleles"]
  ## compare allele weights
  allele_weight_samples = parameter_samples["rescaledAlleleWeights"]
  allele_weights = parameters["AlleleWeights"]
  ## compare conditional allele frequencies for missing alleles
  adm_allele_prob_samples = parameter_samples["admAlleleProbs"]
  adm_allele_prob_means = aux.un_flatten(NumAdmAlleles, np.mean(adm_allele_prob_samples, axis=0))
  adm_alleles = aux.un_flatten(NumAdmAlleles, data["AdmAlleles"])
  nontriv_adm_allele_prob_means = [d for xss in adm_allele_prob_means
                                   for xs in xss for d in xs if len(d) > 1]
  ## find the 'real' conditional admissible allele probabilities
  adm_allele_probs = []
  for s in range(NumSubjects):
    xss = []
    for p in range(Ploidy):
      xs = []
      for ell in range(NumLoci):
        aap = np.array([allele_freqs[ell][a-1] for a in adm_alleles[s][p][ell]])
        aap /= np.sum(aap) ## normalize to get conditional probabilities
        xs.append(aap)
      xss.append(xs)
    adm_allele_probs.append(xss)
  nontriv_adm_allele_probs = [d for xss in adm_allele_probs
                              for xs in xss for d in xs if len(d) > 1]
  ## create a figure object
  fig, axs = plt.subplots(3, 1, figsize=(10,10))
  ## plot (estimated) allele frequencies
  pos = range(1, sum(NumAlleles)+1)
  axs[0].violinplot(allele_freq_samples)
  axs[0].scatter(pos, aux.flatten(allele_freqs), color='k', s=10)
  axs[0].set_ylabel("allele frequencies")
  axs[0].set_xlim(0, sum(NumAlleles)+1)
  ## plot (estimated) allele weights
  axs[1].violinplot(allele_weight_samples, positions=pos)
  axs[1].scatter(pos, aux.flatten(allele_weights), color='k', s=10)
  axs[1].set_ylabel("allele weights")
  axs[1].set_xlim(0, sum(NumAlleles)+1)
  axs[1].axhline(y=0, color='red', linestyle='--')
  ## separate individual alleles
  for ax in axs[:2]: ## first two panels require locus separation
    for ell in range(1, NumLoci): ## when there are 3 loci, plot only 2 lines
      ax.axvline(x=np.sum(NumAlleles[:ell])+0.5, color='k', linestyle='--')
  ## plot conditional allele probabilities
  for s in range(len(nontriv_adm_allele_probs)):
    ps = nontriv_adm_allele_probs[s]
    qs = nontriv_adm_allele_prob_means[s]
    for i in range(len(ps)):
      axs[2].plot([s, s+1], [np.sum(ps[:i]), np.sum(qs[:i])], color='k')
      axs[2].axvline(x=s, color='k')
  axs[2].set_ylabel("probability")
  axs[2].set_ylim(0, 1)
  axs[2].set_xlim(0, len(nontriv_adm_allele_probs))
  ## save the resulting figure...
  fig.savefig("../stan_model_validation.png", dpi=200, bbox_inches='tight')


def validateModel(): ## TODO: iter, chains, ...
  ## compile the Stan model
  with open("stan/hla-tree-model.stan") as f: model_code = f.read() ## FIXME: pass filename
  sm = stantools.CachedStanModel(model_code=model_code, path="stan/cache")
  ## create some random data
  data, parameters = simulateData(WBIC=0, MeanNumAlleles=5, NumSubjects=100)
  ## determine which parameters to follow
  monitor = [
    "rescaledNodeWeights",
    "rescaledAlleleWeights",
    "sigmaTraitValue",
    "sigmaNodeWeight",
    "alleleFreqs",
    "admAlleleProbs",
    "traitValueLoglikes",## for WAIC
    "sumTraitValueLoglikes" ## for WBIC
  ]
  fit = sm.sampling(data=data, pars=monitor, chains=4, iter=2000)
  ## print fit statistics and plot a validation figure
  print("Stan fit statistics:", fit)
  chain = fit.extract()
  statistics.calcWAIC(aux.transpose(chain["traitValueLoglikes"]), verbose=True)
  validateModelFit(fit, parameters, data)
  ## now sample to compute WBIC
  data["WBIC"] = 1
  monitor_wbic = ["sumTraitValueLoglikes"]
  fit_wbic = sm.sampling(data=data, pars=monitor_wbic, chains=4, iter=2000)
  print("Stan WBIC fit statistics:", fit_wbic)
  chain_wbic = fit_wbic.extract()
  statistics.calcWBIC(chain_wbic["sumTraitValueLoglikes"], verbose=True)



if __name__ == '__main__':
  validateModel()
  sys.exit(0)
