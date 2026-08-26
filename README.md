# Higgs Boson Rediscovery: Bayesian Inference and Discovery Sensitivity
## A statistical analysis for the $H \rightarrow ZZ^* \rightarrow 4 \ell$ and $H \rightarrow \gamma \gamma$ Decay Channel

## 1. Overview
This project investigates the rediscovery of the Higgs boson using 13 TeV proton-proton collision data from the LHC, available through the CERN Open Data Portal. Using Run 2 data collected from 2015–2018, we analyze two Higgs decay channels: the four-lepton ($H\rightarrow ZZ^*\rightarrow4\ell$) “golden channel” and the diphoton ($H\rightarrow\gamma\gamma$) channel.

We impose physics-motivated and detector-level criteria for selection to isolate events that are consistent with each decay channel and reconstruct their corresponding invariant-mass distributions. The cut data are modeled with a signal-plus-background Poisson likelihood, which is combined with a prior on the signal strength parameter to perform Bayesian inference using BAT.jl.

The project then utilizes toy Monte Carlo (pseudo-)experiments to evaluate discovery sensitivity by comparing the background-only (no Higgs) hypothesis with pseudo-data containing an injected Higgs signal. The goal is to determine whether the analysis can provide sufficient evidence to establish the presence of the Higgs boson with $5\sigma$ significance against the background-only hypothesis. Future progress will extend this statistical framework to perform limit setting on the signal strength.

# 2. Data
## 2.1 LHC Run 2 Data
The analysis uses publicly available proton-proton collision data from the ATLAS Collaboration through the CERN Open Data Portal. Both datasets were collected during the 2016 LHC data-taking period at a center-of-mass energy of 13 TeV and were released as part of the ATLAS 2020 Open Data release for educational use. The datasets contain both real collision data and corresponding simulated samples of Standard Model processes and selected Beyond the Standard Model signals.

Two complementary datasets are used in this analysis:

Four-lepton channel: The ATLAS 13 TeV samples collection at least four leptons (electron or muon) contains events preselected to include at least four electrons or muons. The collection consists of 111 files totaling approximately 930.5 MiB.
Diphoton channel: The ATLAS 13 TeV samples collection Gamma-Gamma contains events preselected to include at least two photons. The collection consists of 10 files totaling approximately 3.5 GiB.

Both collections apply a loose preselection at the object and event levels to reduce the number of events requiring further analysis. The selections described in Section 3 are subsequently applied to these preselected samples to isolate events consistent with the respective Higgs decay channels. Documentation released by CERN for the properties of each event and its identification in the raw code is provided in the references for both the four-lepton and diphoton data files. 

## 2.2 Four-Lepton Dataset ($H \rightarrow ZZ^* \rightarrow 4 \ell$)
The four-lepton analysis uses the ATLAS 13 TeV samples collection at least four leptons (electron or muon) from the 2020 Open Data release. The dataset contains both real collision data and simulated Monte Carlo samples consistent with Standard Model processes and Higgs signal events.

For the $H\rightarrow ZZ^*\rightarrow4\ell$ analysis, we use events containing combinations of electrons and muons and their associated kinematic and identification information. The primary quantities used in the analysis include each lepton's transverse momentum ($p_T$), pseudorapidity ($\eta$), azimuthal angle ($\phi$), energy ($E$), charge, lepton type, and identification and isolation variables.

The real collision data provide the observed event sample, while simulated signal and background samples are used to model the expected contributions to the four-lepton invariant-mass distribution. The dataset has already undergone a loose ATLAS preselection requiring at least four leptons; additional event-selection criteria are applied in this analysis to isolate events consistent with the $H\rightarrow ZZ^*\rightarrow4\ell$ decay.

## 2.3 Diphoton Dataset ($H \rightarrow \gamma\gamma$)
The diphoton analysis uses the ATLAS 13 TeV samples collection Gamma-Gamma from the 2020 Open Data release. The dataset also contains real collision data together with Monte Carlo simulated samples of Standard Model processes and selected signal processes.

For the $H\rightarrow\gamma\gamma$ analysis, we use events containing at least two photons and their associated kinematic and identification information. The primary quantities used in the analysis include photon transverse momentum ($p_T$), pseudorapidity ($\eta$), azimuthal angle ($\phi$), energy ($E$), and photon identification and isolation variables.

The real collision data provide the observed event sample, while the simulated samples are used to characterize the expected signal and background contributions where applicable. The dataset has already undergone a loose ATLAS preselection requiring at least two photons; additional selection criteria are applied in this analysis to isolate events consistent with the $H\rightarrow\gamma\gamma$ decay.

# 3. Event Selection
## 3.1 Physics Motivation
## 3.2 Four-Lepton Selection Criteria / Cuts
## 3.3 Diphoton Selection Criteria / Cuts

# 4. Event Reconstruction
## 4.1 Four-Lepton ("Golden") Channel
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

