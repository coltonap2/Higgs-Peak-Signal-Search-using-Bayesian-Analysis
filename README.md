# Higgs Boson Rediscovery: Bayesian Inference and Discovery Sensitivity
## A statistical analysis for the $H \rightarrow ZZ^* \rightarrow 4 \ell$ and $H \rightarrow \gamma \gamma$ Decay Channel

## 1. Overview
This project develops a reproducible statistical analysis workflow for particle-physics discovery searches, demonstrated through the rediscovery of the Higgs boson using publicly available 13 TeV proton-proton collision data from the CERN Open Data Portal. The analysis uses Run 2 data collected from 2015–2018 and considers two Higgs decay channels: the four-lepton ($H \rightarrow ZZ^* \rightarrow 4\ell$) “golden channel” and the diphoton ($H \rightarrow \gamma\gamma$) channel.

The workflow proceeds from publicly available collision data through event selection, invariant-mass reconstruction, statistical modeling, and hypothesis testing. Physics-motivated and detector-level selection criteria are applied to isolate events consistent with each decay channel, after which the resulting invariant-mass distributions are used to construct signal and background models.

At the core of the analysis is a signal-plus-background statistical model in which the expected number of events in each mass bin is parameterized by the signal strength $\mu$ and the corresponding signal and background contributions:

$$
\lambda_i(\mu) = \mu s_i + b_i,
$$

where $\lambda_i(\mu)$ is the expected number of events in bin $i$, $s_i$ is the expected signal contribution, and $b_i$ is the expected background contribution. The resulting Poisson likelihood is combined with a prior on the signal strength parameter to perform Bayesian inference using BAT.jl.

The same statistical structure is applied across the two decay channels, while the channel-specific event selections and signal and background models are defined from the corresponding data. This separation between the common statistical model and channel-specific physics inputs allows the analysis workflow to be reproduced and adapted to other datasets or search channels.

To evaluate discovery sensitivity, the framework uses toy Monte Carlo (pseudo-)experiments to generate datasets under both the background-only ($\mu=0$) and signal-plus-background ($\mu>0$) hypotheses. These pseudo-experiments are then analyzed using the same inference procedure applied to the observed data, allowing the expected sensitivity of the analysis to be studied and compared with the evidence required for a $5\sigma$ discovery.

The repository is organized as a reproducible workflow rather than solely as a record of the final result. Each stage of the analysis—from data preparation and event selection to statistical modeling and pseudo-experiment generation—is documented so that the analysis can be reproduced, modified, and extended. Although the Higgs boson serves as the primary physics case study, the underlying statistical workflow is intended to provide a foundation that can be adapted to other counting-based particle-physics searches and datasets.

Future development will extend the framework to signal-strength limit setting and further generalize the statistical and modeling components for use beyond the Higgs rediscovery analysis.

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
This crucial step reduces noise from irrelevant detections by eliminating events that fail to meet the desired criteria.  The ultimate goal is to produce an invariant mass distribution, which can then be fit and analyzed for the desired Higgs spike.  The two detection channels each require a distinct set of cuts, as detailed below:

## 3.2 Four-Lepton Selection Criteria / Cuts
The four-lepton channel is treated with a series of six cuts.  These cuts are motivated by the selection criteria laid out in the original Higgs Boson discovery paper, and are summarized as follows:

i. Transverse momenta

Electrons are identified as having transverse momenta $p_t > 7$ GeV and pseudorapidity $|\eta|< 2.47$ , while muons require a minimum transverse momentum of 5 GeV and pseudorapidity $|\eta| < 2.7$.  Leptons must meet additional thresholds based on the order of detection.  Specifically, the first must have a transverse momentum of $p_t > 20$ GeV, the second > 15 GeV, and the third > 10 GeV.

ii. Charge/flavor

Charge conservation requires that each event must contain two pairs of leptons of opposite sign (charge) and same flavor (OSSF).

iii. Mass pairings

The OSSF pair whose invariant mass is closest to that of the $Z$ boson (91.18 GeV) is chosen as the leading pair $Z_1$, with the other denoted as the subleading pair $Z_2$

iv. Mass bounds

The leading pair must have invariant mass between 50 and 106 GeV, the subleading pair between 12 and 115 GeV

v. Resonance veto

To eliminate similar events caused by $J/\psi$ and $\Upsilon$ meson decays, we require that all OSSFs have invariant mass $m_ll > 5$ GeV

vi. Lepton separation

Finally, to remove any overlapping detections from closely-spaced leptons, we require a minimum separation of 0.1 for OSSF pairs and 0.2 for differently-flavored pairs.

## 3.3 Diphoton Selection Criteria / Cuts
The noisier diphoton channel requires eight cuts, some of which are already flagged in the CERN dataset.  These cuts are made in accordance with the selection criteria from a newer 2022 paper for $\sqrt{s} = 13$ TeV:

i. Photon trigger

Selects all events already identified by CERN's photon trigger

ii. Pseudorapidity

Requires a pseudorapidity of $|\eta| < 2.37$ for each photon

iii. Crack veto

Excludes all photons whose with pseudorapidity within “crack” caused by detector’s geometry, namely the range $1.37 < |\eta| < 1.52$

iv. Transverse momentum

Requires minimum transverse momentum of $p_t > 40$ for leading photon, 30 GeV for subleading photon

v. Tight photon identification

Excludes all candidates that conflict with CERN’s photon identification criteria

vi. Photon isolation

Requires maximum track isolation of 5%, calorimetric isolation less than 6.5%

vii. Transverse energy to invariant mass ratio

Requires a minimum $E_T/m_\gamma\gamma$ ratio of 0.35 for leading, 0.25 for subleading photon

viii. Diphoton mass

Limits diphoton candidates to mass range between 105 and 160 GeV

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

ATLAS Collaboration (2025). ATLAS 13 TeV Open Data: 2025 Beta Release. CERN Open Data Portal. https://opendata.atlas.cern/docs/category/13-tev-2025-beta-release
