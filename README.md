# Higgs Boson Rediscovery: Bayesian Inference and Discovery Sensitivity
## A statistical analysis for the $H \rightarrow ZZ^* \rightarrow 4 \ell$ and $H \rightarrow \gamma \gamma$ Decay Channel

## 1. Overview
This project investigates the rediscovery of the Higgs boson using 13 TeV proton-proton collision data from the LHC, available through the CERN Open Data Portal. Using Run 2 data collected from 2015–2018, we analyze two Higgs decay channels: the four-lepton ($H\rightarrow ZZ^*\rightarrow4\ell$) “golden channel” and the diphoton ($H\rightarrow\gamma\gamma$) channel.

We impose physics-motivated and detector-level criteria for selection to isolate events that are consistent with each decay channel and reconstruct their corresponding invariant-mass distributions. The cut data are modeled with a signal-plus-background Poisson likelihood, which is combined with a prior on the signal strength parameter to perform Bayesian inference using BAT.jl.

The project then utilizes toy Monte Carlo (pseudo-)experiments to evaluate discovery sensitivity by comparing the background-only (no Higgs) hypothesis with pseudo-data containing an injected Higgs signal. The goal is to determine whether the analysis can provide sufficient evidence to establish the presence of the Higgs boson with $5\sigma$ significance against the background-only hypothesis. Future progress will extend this statistical framework to perform limit setting on the signal strength.

# 2. Data
## 2.1 LHC Run 2 Data
## 2.2 Four-Lepton Dataset ($H \rightarrow ZZ^* \rightarrow 4 \ell$)
## 2.3 Diphoton Dataset ($H \rightarrow \gamma\gamma$)

# 3. Event Selection
## 3.1 Physics Motivation
## 3.2 Four-Lepton Selection Criteria / Cuts
## 3.3 Diphoton Selection Criteria / Cuts

# 4. Event Reconstruction
## 4.1 Four Lepton ("Golden") Channel
## 4.2 Diphoton Channel

# 5. Statistical Model
## 5.1 Common Structure of the Model
## 5.2 Four-Lepton Statistical Model
## 5.3 Diphoton Statistical Model

# 6. Bayesian Inference
## 6.1 Prior on Signal Strength
## 6.2 Posterior Distribution
## 6.3 BAT.jl Implementation
## 6.4 Four-Lepton Results
## 6.5 Diphoton Results

# 7. Discovery Sensitivity
## 7.1 Toy Monte Carlo Experiments
## 7.2 Background-Only Hypthesis
## 7.3 Signal-Plus-Background Hypothesis
## 7.4 Discovery Test Statistic
## 7.5 $5\sigma$ Sensitivity
## 7.6 Four-Lepton Results
## 7.7 Diphoton Results

## References
ATLAS Collaboration (2020). ATLAS 13 TeV samples collection at least four leptons (electron or muon), for 2020 Open Data release. CERN Open Data Portal. DOI:10.7483/OPENDATA.ATLAS.2Y1T.TLGL

ATLAS Collaboration (2020). ATLAS 13 TeV samples collection Gamma-Gamma, for 2020 Open Data release. CERN Open Data Portal. DOI:10.7483/OPENDATA.ATLAS.B5BJ.3SGS
