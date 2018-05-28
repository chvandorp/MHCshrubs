/* fit the HLA tree model with Stan, taking missing HLAs into account.
 * ------------------------------------------------------------------
 * common abbreviations:
 * ---------------------
 * Num: number
 * Adm: admissible
 * Comb: combination
 * Idx: index
 * Bnd: boundary/bound
 * Addl: Additional
 * ------------------------------------------------------------------
 * codes for options:
 * ---------------------
 * regular sampling: set WBIC = 0
 * sampling to compute WBIC: set WBIC = 1
 * normal prior for node weights: set NodeWeightPrior = 0
 * laplace prior for node weights: set NodeWeightPrior = 1
 * ------------------------------------------------------------------
 * TODO:
 * add a vector of loglikes corresponding to observations of allele frequencies
 * This vector must then be used for the WAIC and WBIC computation
 * The WatanabeBeta-scaled sum must be added to target.
 */

functions {
  vector get_simplex(vector u) {
    /* turns a positive vector u of length k-1 into a k-simplex.
     * the simplex is represented as a vector
     */
    vector[num_elements(u)] v; // v is an ordered vector
    // make a ordered vector v
    v = cumulative_sum(u);
    // transform to (0, 1)
    v = v ./ (1+v); // note the space between "v" and "./"
    // take differences to get a simplex
    return append_row(v, rep_vector(1.0, 1)) - append_row(rep_vector(0.0, 1), v);
  }
  real get_simplex_logdetjac(vector u) {
    /* compute the log determinant of the jacobian
     * for the above transformation get_simplex
     *
     * TODO tripple check that this is correct
     */
    return -2.0 * sum(log(1+cumulative_sum(u)));
  }
  int count_combs(int[,,] nss) {
    int combs[size(nss)];

   for ( i in 1:size(nss) ) {
      combs[i] = prod(to_array_1d(nss[i,:,:]));
    }
    return sum(combs);
  }
}

data {
  int<lower=0> Ploidy; // usually equal to 2
  int<lower=0> NumLoci; // e.g. equals 3 for HLA A, B, C
  int<lower=0> NumSubjects; // number of subjects
  real TraitValue[NumSubjects]; // disease trait
  int<lower=0,upper=3> TraitCensorType[NumSubjects]; // 0 = uncensored, 1 = left censored, 2 = right censored, 3 = missing
  int<lower=0> NumAlleles[NumLoci]; // number of alleles
  int<lower=1> NumAdmAlleles[NumSubjects, Ploidy, NumLoci]; // each subject must at least have one admissible allele
  int<lower=1> AdmAlleles[sum(to_array_1d(NumAdmAlleles))]; // all possible alleles concatenated into one big array
  int<lower=0> AddlAlleleData[sum(NumAlleles)]; // additional allele counts
  int<lower=sum(NumAlleles)> NumNodes; // nodes in the allele tree
  matrix<lower=0, upper=1>[sum(NumAlleles), NumNodes] TreeMatrix; // binary matrix encoding the tree
  vector<lower=0>[NumNodes] LengthEdgeToParent; // correct weights with distance
  // some options...
  int<lower=0, upper=1> WBIC; // 0 -> regular sampling, 1 -> WBIC sampling
  int<lower=0, upper=1> NodeWeightPrior; // 0 -> normal, 1 -> laplace
}

transformed data {
  int<lower=0> SumNumAlleles; // the total number of alleles
  int<lower=1> LeftAlleleBnds[NumLoci]; // alleleFreqs and alleleWeights
  int<lower=1> RightAlleleBnds[NumLoci];
  int<lower=1> LeftAlleleFreqParamBnds[NumLoci]; // alleleFreqParams
  int<lower=1> RightAlleleFreqParamBnds[NumLoci];
  int<lower=0> SumNumAdmAlleles;
  int<lower=1> LeftAdmAlleleBnds[NumSubjects, Ploidy, NumLoci]; // admAlleleProbs, AdmAlleles, AdmAlleleIdxs
  int<lower=1> RightAdmAlleleBnds[NumSubjects, Ploidy, NumLoci];
  int<lower=1> LeftAdmAlleleProbParamBnds[NumSubjects, Ploidy, NumLoci]; // admAlleleProbParams
  int<lower=0> RightAdmAlleleProbParamBnds[NumSubjects, Ploidy, NumLoci];
  int<lower=1> AdmAlleleIdxs[sum(to_array_1d(NumAdmAlleles))];
  // a matrix listing all possible allele combinations
  int<lower=1> LeftAdmAlleleCombBnds[NumSubjects]; // AdmAlleleIdxCombs, AdmAlleleProbIdxCombs
  int<lower=1> RightAdmAlleleCombBnds[NumSubjects];
  int<lower=1> AdmAlleleIdxCombs[count_combs(NumAdmAlleles), Ploidy, NumLoci]; // alleleFreqs, alleleWeights
  int<lower=1> AdmAlleleProbIdxCombs[count_combs(NumAdmAlleles), Ploidy, NumLoci]; // admAlleleProbs
  // centralize and scale the TraitValue
  real LocTraitValue;
  real ScaleTraitValue;
  real NormTraitValue[NumSubjects];
  int NumUncensoredTraitValues;
  real<lower=0> WatanabeBeta; // determines sampling temperature: 1/log(N), where N is the number of observations

  // define integer constants
  SumNumAlleles = sum(NumAlleles);
  SumNumAdmAlleles = sum(to_array_1d(NumAdmAlleles));
  print("INFO: ", SumNumAdmAlleles, " adm. alleles / ", NumSubjects*Ploidy*NumLoci, " alleles");

  // define indices
  for ( ell in 1:NumLoci ) {
    LeftAlleleBnds[ell] = sum(NumAlleles[:ell-1]) + 1; // indexing starts at 1
    RightAlleleBnds[ell] = sum(NumAlleles[:ell]);

    LeftAlleleFreqParamBnds[ell] = sum(NumAlleles[:ell-1]) - (ell - 1) + 1;
    RightAlleleFreqParamBnds[ell] = sum(NumAlleles[:ell]) - ell;
  }

  // define indices for admissible alleles
  for ( s in 1:NumSubjects ) {
    for ( p in 1:Ploidy ) {
      for ( ell in 1:NumLoci ) {
        int i;
        int al; int ar;

        i = (s-1)*Ploidy*NumLoci + (p-1)*NumLoci + ell;
        al = sum(to_array_1d(NumAdmAlleles)[:i-1]) + 1;
        ar = sum(to_array_1d(NumAdmAlleles)[:i]);

        LeftAdmAlleleBnds[s, p, ell] = al;
        RightAdmAlleleBnds[s, p, ell] = ar;

        LeftAdmAlleleProbParamBnds[s, p, ell] = al - i + 1;
        RightAdmAlleleProbParamBnds[s, p, ell] = ar - i;

        for ( a in al:ar ) { // find index in alleleFreqs
          AdmAlleleIdxs[a] = AdmAlleles[a] + sum(NumAlleles[:ell-1]);
        }
      }
    }

    // construct boundaries of all possible allele combinations
    LeftAdmAlleleCombBnds[s] = count_combs(NumAdmAlleles[:s-1,:,:]) + 1;
    RightAdmAlleleCombBnds[s] = count_combs(NumAdmAlleles[:s,:,:]);
  }

  // construct all possible allele combinations
  for ( s in 1:NumSubjects ) {
    int sl; int sr;

    sl = LeftAdmAlleleCombBnds[s];
    sr = RightAdmAlleleCombBnds[s];
    for ( p in 1:Ploidy ) {
      for ( ell in 1:NumLoci ) {
        int al; int ar;

        al = LeftAdmAlleleBnds[s, p, ell]; ar = RightAdmAlleleBnds[s, p, ell];
        for ( a in al:ar ) {
          for ( i in 1:(sr-sl+1)/(ar-al+1) ) {
            AdmAlleleIdxCombs[sl + (i-1)*(ar-al+1) + (a-al), p, ell] = AdmAlleleIdxs[a];
            AdmAlleleProbIdxCombs[sl + (i-1)*(ar-al+1) + (a-al), p, ell] = a;
          }
        } // loop over alleles
      } // loop over loci
    } // loop over ploidy
  } // loop over subjects

  // compute normalized TraitValue
  LocTraitValue = 0.0;
  ScaleTraitValue = 0.0;
  NumUncensoredTraitValues = 0;

  for ( s in 1:NumSubjects ) {
    if ( TraitCensorType[s] == 0 ) {
      LocTraitValue += TraitValue[s];
      NumUncensoredTraitValues += 1;
    }
  }
  if ( NumUncensoredTraitValues > 0 ) {
    LocTraitValue /= NumUncensoredTraitValues;
  }

  for ( s in 1:NumSubjects ) {
    if ( TraitCensorType[s] == 0 ) {
      ScaleTraitValue += square(TraitValue[s] - LocTraitValue);
    }
  }
  if ( NumUncensoredTraitValues > 1 ) {
    ScaleTraitValue = sqrt(ScaleTraitValue / (NumUncensoredTraitValues-1));
  } else {
   ScaleTraitValue = 1.0;
  }

  for ( s in 1:NumSubjects ) {
    NormTraitValue[s] = (TraitValue[s] - LocTraitValue) / ScaleTraitValue;
  }

  // determine the sampling temperature
  if ( WBIC == 1 ) {
    print("USING SAMPLING TEMPERATURE FOR WBIC");
    WatanabeBeta = 1.0 / log(NumUncensoredTraitValues);
  } else {
    print("USING REGULAR SAMPLING TEMPERATURE");
    WatanabeBeta = 1.0;
  }
}

parameters {
  vector<lower=0>[SumNumAlleles-NumLoci] alleleFreqParams; // must be transformed to get allele frequencies
  vector<lower=0>[SumNumAdmAlleles-NumSubjects*Ploidy*NumLoci] admAlleleProbParams; // must be transformed to get likelihood
  real<lower=0> sigmaNodeWeight; // hypo-parameter
  vector[NumNodes] nodeWeights; // must be transformed to make allele weights
  real<lower=0> sigmaTraitValue; // hypo-parameter
}

transformed parameters {
  vector<lower=0, upper=1>[SumNumAlleles] alleleFreqs; // concatenated simplices
  vector<lower=0, upper=1>[SumNumAdmAlleles] admAlleleProbs; // likelihood of admissible allele for each subject, locus, zygote
  vector[SumNumAlleles] alleleWeights; // effects of alleles on the TraitValue
  vector[NumSubjects] traitValueLoglikes; // used for WAIC and should be added to target in the model block

  // compute allele frequencies
  for ( ell in 1:NumLoci ) {
    int al; int ar;
    int bl; int br;

    al = LeftAlleleBnds[ell]; ar = RightAlleleBnds[ell]; // alias
    bl = LeftAlleleFreqParamBnds[ell]; br = RightAlleleFreqParamBnds[ell];
    alleleFreqs[al:ar] = get_simplex(alleleFreqParams[bl:br]);
  }

  // compute probabilities of admissible alleles
  for ( s in 1:NumSubjects ) {
    for ( p in 1:Ploidy ) {
      for ( ell in 1:NumLoci ) {
        int al; int ar;
        int bl; int br;

        al = LeftAdmAlleleBnds[s, p, ell]; ar = RightAdmAlleleBnds[s, p, ell];
        bl = LeftAdmAlleleProbParamBnds[s, p, ell]; br = RightAdmAlleleProbParamBnds[s, p, ell];
        admAlleleProbs[al:ar] = get_simplex(admAlleleProbParams[bl:br]);
      }
    }
  }

  // compute the allele weights using the tree
  alleleWeights = TreeMatrix * (LengthEdgeToParent .* nodeWeights);

  // model the TraitValue
  for ( s in 1:NumSubjects ) {
    int sl; int sr;
    // likelihoods of the TraitValue, given combination of alleles
    vector[RightAdmAlleleCombBnds[s] - LeftAdmAlleleCombBnds[s] + 1] loglikes;
    // probability for the allele combination
    vector[RightAdmAlleleCombBnds[s] - LeftAdmAlleleCombBnds[s] + 1] logprobs;

    sl = LeftAdmAlleleCombBnds[s]; sr = RightAdmAlleleCombBnds[s];
    for ( sc in sl:sr ) {
      real sum_allele_weights;

      logprobs[sc-sl+1] = log(prod(admAlleleProbs[to_array_1d(AdmAlleleProbIdxCombs[sc,:,:])]));
      sum_allele_weights = sum(alleleWeights[to_array_1d(AdmAlleleIdxCombs[sc,:,:])]);

      if ( TraitCensorType[s] == 0 ) { // uncensored
	       loglikes[sc-sl+1] = normal_lpdf(NormTraitValue[s] | sum_allele_weights, sigmaTraitValue);
      } else if ( TraitCensorType[s] == 1 ) { // left censored
	       loglikes[sc-sl+1] = normal_lcdf(NormTraitValue[s] | sum_allele_weights, sigmaTraitValue);
      } else if ( TraitCensorType[s] == 2 ) { // right censored
	       loglikes[sc-sl+1] = normal_lccdf(NormTraitValue[s] | sum_allele_weights, sigmaTraitValue);
      } else { // missing value
	       loglikes[sc-sl+1] = 0.0;
      }
    }
    traitValueLoglikes[s] = log_sum_exp(loglikes + logprobs);
  }
}

model {
  for ( ell in 1:NumLoci ) {
    int al; int ar;

    al = LeftAlleleBnds[ell]; ar = RightAlleleBnds[ell]; // alias
    // Dirichlet priors for the allele frequencies
    alleleFreqs[al:ar] ~ dirichlet(rep_vector(0.5, NumAlleles[ell])); // Jeffreys prior ?
    // correct mor manual simplex transformations
    target += get_simplex_logdetjac(alleleFreqParams[al-ell+1:ar-ell]);
    // Additional Allele frequency data
    AddlAlleleData[al:ar] ~ multinomial(alleleFreqs[al:ar]);
  }

  // handle missing alleles and inform allele frequencies by subjects
  for ( s in 1:NumSubjects ) {
    for ( p in 1:Ploidy ) {
      for ( ell in 1:NumLoci ) {
        int al; int ar;
        int bl; int br;

        al = LeftAdmAlleleBnds[s, p, ell]; ar = RightAdmAlleleBnds[s, p, ell];
        bl = LeftAdmAlleleProbParamBnds[s, p, ell]; br = RightAdmAlleleProbParamBnds[s, p, ell];
        // likelihood of admissible alleles, given conditional allele frequencies
        admAlleleProbs[al:ar] ~ dirichlet(alleleFreqs[AdmAlleleIdxs[al:ar]] / sum(alleleFreqs[AdmAlleleIdxs[al:ar]]));
        // count frequencies
        target += log(dot_product(admAlleleProbs[al:ar], alleleFreqs[AdmAlleleIdxs[al:ar]]));
        // correct for manual simplex transformation
        target += get_simplex_logdetjac(admAlleleProbParams[bl:br]);
      }
    }
  }

  // prior for the nodeWeight
  sigmaNodeWeight ~ normal(0.0, 10.0); // sigma is restricted to be positive
  if ( NodeWeightPrior == 0 ) { // 0 -> normal
    nodeWeights ~ normal(0.0, sigmaNodeWeight);
  } else if ( NodeWeightPrior == 1 ) { // 1 -> laplace
    nodeWeights ~ double_exponential(0.0, sigmaNodeWeight);
  } else {
    print("ERROR: invalid NoneWeight prior");
  }

  // prior for mean and sd TraitValue
  sigmaTraitValue ~ normal(0.0, 10.0); // sigma is restricted to be positive

  // add loglikes for the TraitValue to the target
  target += WatanabeBeta * sum(traitValueLoglikes);
}

generated quantities {
  vector[NumNodes] rescaledNodeWeights; // transform node weights back to TraitValue scale
  vector[SumNumAlleles] rescaledAlleleWeights; // transform node weights back to TraitValue scale
  real sumTraitValueLoglikes; // for WBIC computation

  rescaledNodeWeights = ScaleTraitValue * nodeWeights;
  rescaledAlleleWeights = ScaleTraitValue * alleleWeights;
  sumTraitValueLoglikes = sum(traitValueLoglikes);
}
