# new_asymptotic_formulae
This is to document some codes to apply the new asymptotic formulae proposed in the paper [2101.06944](https://arxiv.org/abs/2101.06944) in a search for a bottom quark supersymmetric partner by the ATLAS collaboration [JHEP12(2019)060](https://arxiv.org/abs/1908.03122). The public likelihoods can be found from [HepData](https://www.hepdata.net/record/ins1748602).

## Step1: do fits and produce toys
The toys are produced using the [PyHF](https://joss.theoj.org/papers/10.21105/joss.02823) framework. We also need to do some fits to get the number of signal and background events to build the binned model proposed in [2101.06944](https://arxiv.org/abs/2101.06944) to estimate the effect of the low statistics. 

## Step2: Analysis of the toys
In the directory "bottom-squarks_application", we have some python scripts for analyzing the toys and compare with the classic and new asymptotic formulae.
