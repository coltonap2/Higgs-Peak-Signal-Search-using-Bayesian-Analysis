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
The analysis begins with publicly available proton-proton collision data from the ATLAS Collaboration through the CERN Open Data Portal. The datasets used in this analysis were collected during the 2016 LHC data-taking period at a center-of-mass energy of $13$ TeV and were released as part of the ATLAS 2020 Open Data release for educational use.

The Open Data release provides both real collision data and simulated Monte Carlo samples. The collision data contain events recorded by the ATLAS detector, while the simulated samples contain modeled events representing Standard Model processes and selected signal processes. These samples are distinguished in the released datasets through the dataset variable and are used differently in later stages of the analysis.

This analysis uses two complementary ATLAS data collections:

Four-lepton channel: The [ATLAS 13 TeV samples collection containing at least four leptons (electron or muon)] provides events preselected to contain at least four electrons or muons. The collection consists of 111 files, totaling approximately 930.5 MiB.
Diphoton channel: The [ATLAS 13 TeV $\gamma\gamma$ samples collection] provides events preselected to contain at least two photons. The collection consists of 10 files, totaling approximately 3.5 GiB.

The released collections therefore serve as the starting point for the analysis rather than the final datasets used for inference. The ATLAS Open Data samples have already undergone a loose preselection designed to identify events containing the relevant reconstructed objects. The analysis then reads these released event-level variables and applies the additional processing described in Sections 3 and 4.

The overall relationship between the data and subsequent analysis stages is:

$$ \text{ATLAS Open Data} \rightarrow \text{Data Preparation} \rightarrow \text{Event Selection} \rightarrow \text{Event Reconstruction} \rightarrow \text{Statistical Analysis}. $$

The remainder of this section describes how the two released datasets are organized and identifies the event-level variables used by the analysis. The specific selection criteria and reconstruction procedures are described separately in Sections 3 and 4.

Official ATLAS documentation for the datasets and their event-level properties is provided in the References. These resources should be consulted when reproducing the analysis or adapting the workflow to another ATLAS Open Data sample.

## 2.2 Four-Lepton Dataset ($H \rightarrow ZZ^* \rightarrow 4 \ell$)
The four-lepton analysis uses the ATLAS 13 TeV samples collection containing at least four leptons (electron or muon) from the 2020 Open Data release. The collection contains both real collision data and simulated Monte Carlo samples representing relevant background and signal processes.

Understanding the Event Data

Each event in the released four-lepton collection contains information about the reconstructed leptons in that event. The dataset provides separate variables for up to four leptons, with the suffix _1 through _4 identifying the corresponding lepton.

For example, lep_pt_1 refers to the transverse momentum of the first lepton, while lep_pt_4 refers to that of the fourth. This naming convention allows the same quantity to be accessed for each reconstructed lepton.

The primary lepton-level variables used by this analysis are:

| Variable                        | Meaning                     | Purpose                                                           |
| :------------------------------ | :-------------------------- | :---------------------------------------------------------------- |
| `lep_pt_1` – `lep_pt_4`         | Transverse momentum ($p_T$) | Provides the transverse momentum of each reconstructed lepton     |
| `lep_eta_1` – `lep_eta_4`       | Pseudorapidity ($\eta$)     | Provides the direction of each lepton relative to the beam axis   |
| `lep_phi_1` – `lep_phi_4`       | Azimuthal angle ($\phi$)    | Provides the angular position of each lepton around the beam axis |
| `lep_E_1` – `lep_E_4`           | Energy ($E$)                | Provides the energy of each reconstructed lepton                  |
| `lep_charge_1` – `lep_charge_4` | Electric charge             | Identifies the charge associated with each lepton                 |
| `lep_type_1` – `lep_type_4`     | Lepton type                 | Distinguishes electrons from muons                                |

These variables provide the basic kinematic and identification information required by the later stages of the analysis. The individual quantities are not themselves the final analysis observable; rather, they provide the information from which the event-selection and reconstruction procedures in Sections 3 and 4 operate.

Dataset Organization

The released files contain both observed and simulated events. For the observed-data analysis, events are identified using the dataset classification provided in the file. Simulated samples are retained separately so that they can be used to characterize the expected signal and background contributions in the statistical analysis.

The analysis code reads the released CSV data into a Julia DataFrame. The event-level columns can then be accessed directly by their dataset identifiers. For example: df.lep_pt_1, df.lep_eta_1, and df.lep_charge_1.

This provides a direct connection between the variables documented by ATLAS and the quantities used by the analysis code.

The units of the released variables should also be checked before analysis. In the workflow used here, lepton transverse momenta and energies are converted from MeV to GeV during data preparation where necessary. The corresponding conversion is performed in the analysis code before the variables are used downstream.

The specific criteria applied to these variables are described in Section 3.2, while the construction of the four-lepton observable is described in Section 4.1.

## 2.3 Diphoton Dataset ($H \rightarrow \gamma\gamma$)
The diphoton analysis uses the ATLAS 13 TeV $\gamma\gamma$ samples collection from the 2020 Open Data release. Like the four-lepton collection, it contains real collision data together with Monte Carlo simulated samples representing relevant background and signal processes.

Understanding the Event Data

Each event in the released diphoton collection contains information about the reconstructed photons in that event. Variables associated with the two photons are identified using the suffixes _1 and _2. For example, photon_pt_1 and photon_pt_2 contain the transverse momenta of the two reconstructed photons.

The primary variables used by this analysis are:

| Variable                                   | Meaning                              | Purpose                                                                      |
| :----------------------------------------- | :----------------------------------- | :--------------------------------------------------------------------------- |
| `dataset`                                  | Dataset classification               | Identifies the type of event sample, including real collision data (`data`)  |
| `trigP_1`, `trigP_2`                       | Photon trigger flags                 | Provide the photon-trigger information recorded in the dataset               |
| `photon_pt_1`, `photon_pt_2`               | Transverse momentum ($p_T$)          | Provides the transverse momentum of each reconstructed photon                |
| `photon_eta_1`, `photon_eta_2`             | Pseudorapidity ($\eta$)              | Provides the direction of each photon relative to the beam axis              |
| `photon_phi_1`, `photon_phi_2`             | Azimuthal angle ($\phi$)             | Provides the angular position of each photon around the beam axis            |
| `photon_E_1`, `photon_E_2`                 | Energy ($E$)                         | Provides the energy of each reconstructed photon                             |
| `photon_isTightID_1`, `photon_isTightID_2` | Tight photon identification flags    | Provide the ATLAS photon-identification information                          |
| `photon_ptcone30_1`, `photon_ptcone30_2`   | Track-based isolation variable       | Provides the track activity measured around each photon                      |
| `photon_etcone20_1`, `photon_etcone20_2`   | Calorimeter-based isolation variable | Provides the calorimeter activity measured around each photon                |
| `event_weight`                             | Event weight                         | Provides the weight associated with an event for weighted-yield calculations |
| `m_gg`                                     | Diphoton invariant mass              | Dataset variable corresponding to the diphoton invariant mass                |

These variables provide the event-level information used by the subsequent selection and reconstruction stages. Their specific application is intentionally described in later sections rather than duplicated here.

Dataset Organization

As with the four-lepton collection, the diphoton files contain both observed and simulated events. The dataset variable distinguishes the different samples, allowing the observed collision data to be separated from simulated signal and background samples.

The analysis code reads the released CSV file into a Julia DataFrame, allowing individual event-level variables to be accessed directly. For example: df.photon_pt_1, df.photon_eta_1, df.photon_isTightID_1, and df.photon_ptcone30_1.

The photon transverse momenta and energies are converted from MeV to GeV during data preparation before they are used in the subsequent analysis.

Several variables in the released dataset describe detector-level identification and isolation information. In particular, photon_isTightID, photon_ptcone30, and photon_etcone20 provide information that is later used when applying the diphoton event-selection criteria.

The specific requirements applied to these variables are described in Section 3.3, while the construction of the diphoton observable is described in Section 4.2.

Working with the Released Data

For reproducibility, the important distinction is between variables provided by ATLAS and quantities derived by the analysis. The released dataset supplies the underlying event-level measurements and flags, while later analysis code combines these quantities to produce the selected event sample and reconstructed observables.

This separation allows the same workflow to be adapted to another compatible ATLAS dataset: the dataset-specific variables can be mapped to the corresponding analysis inputs, while the selection, reconstruction, and statistical procedures can be applied independently.

Dataset Documentation

The official ATLAS documentation should be used to verify the meaning, units, and representation of the variables in the released files. The documentation for both the four-lepton and diphoton collections is provided in the References.

This documentation is particularly important when reproducing or modifying the analysis because variable names, units, available samples, and event-level definitions may differ between ATLAS Open Data releases.


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

# 4. Statistical Model
## 4.1 Common Structure of the Model
The statistical model provides the common framework used to analyze the selected invariant-mass distributions in both Higgs decay channels. While the four-lepton and diphoton channels have different event selections and signal and background models, both are reduced to the same statistical problem: comparing the observed number of events in each invariant-mass bin with the number expected from signal and background contributions.

This section describes the statistical structure shared by both channels. The channel-specific construction of the signal and background models is described separately in Sections 5.2 and 5.3.

Binned Signal-Plus-Background Model

After event reconstruction, the invariant-mass distribution is divided into a set of bins. For each bin $i$, we define:

$n_i$ as the observed number of events,
$s_i$ as the expected signal contribution,
$b_i$ as the expected background contribution.

The expected number of events in each bin is modeled as

$$ \lambda_i(\mu) = \mu s_i + b_i, $$

where $\mu$ is the signal strength parameter.

The signal strength scales the expected signal contribution while leaving the background contribution unchanged. This provides a common parameterization for testing different signal hypotheses:

$$ \mu = 0 $$

represents the background-only hypothesis, while

$$ \mu = 1 $$

represents the nominal signal-plus-background hypothesis.

The same parameterization is used for both decay channels. The values of $s_i$ and $b_i$, however, are determined independently for the four-lepton and diphoton analyses and are described in Sections 5.2 and 5.3.

Poisson Likelihood

Because $n_i$ represents a count of observed events in each bin, the number of events in a bin is modeled using a Poisson distribution. For a single bin,

```math
P(n_i|\mu) = \frac{\lambda_i(\mu)^{n_i}e^{-\lambda_i(\mu)}}{n_i!}.
```

Assuming the bins are independent, the likelihood for the complete invariant-mass distribution is

```math
\mathcal{L}(\mu) = \prod_i P(n_i|\mu) = \prod_i \frac{\lambda_i(\mu)^{n_i}e^{-\lambda_i(\mu)}}{n_i!}.
```

Substituting the signal-plus-background model gives

```math
\mathcal{L}(\mu) = \prod_i \frac{(\mu s_i+b_i)^{n_i} e^{-(\mu s_i+b_i)}}{n_i!}.
```

This likelihood is the central statistical object passed to the Bayesian inference procedure in Section 6.

Separating the Statistical Model from the Physics Inputs

An important feature of the framework is the separation between the common statistical structure and the channel-specific physics inputs.

The statistical model requires:

| Quantity | Meaning                               | Source                            |
| :------- | :------------------------------------ | :-------------------------------- |
| $n_i$    | Observed number of events in bin $i$  | Selected collision data           |
| $s_i$    | Expected signal events in bin $i$     | Channel-specific signal model     |
| $b_i$    | Expected background events in bin $i$ | Channel-specific background model |
| $\mu$    | Signal strength                       | Parameter being inferred          |

The procedure for constructing $n_i$, $s_i$, and $b_i$ can differ between analyses without changing the underlying likelihood. This means that the statistical framework can be reused for another decay channel or dataset by replacing the corresponding physics inputs while retaining the same signal-plus-background structure.

For replication, the key interface between the physics analysis and statistical model is therefore the set of binned quantities $(n_i,s_i,b_i)$. Once these quantities have been constructed for a given analysis, the common likelihood can be applied without modifying its fundamental structure.

Model Workflow

The common statistical workflow can be summarized as:

$$ \text{Invariant-Mass Distribution} \rightarrow \text{Binning} \rightarrow (n_i,s_i,b_i) \rightarrow \lambda_i(\mu)=\mu s_i+b_i \rightarrow \mathcal{L}(\mu). $$

The four-lepton and diphoton analyses follow this same statistical structure. Their differences arise in how the signal and background contributions are obtained from their respective datasets, which is described in the following sections.

The Bayesian prior, posterior distribution, and implementation of this likelihood in BAT.jl are intentionally treated separately in Section 6.

## 4.2 Four-Lepton Statistical Model
## 4.3 Diphoton Statistical Model

# 5. Bayesian Inference
## 5.1 Prior on Signal Strength
## 5.2 Posterior Distribution
## 5.3 BAT.jl Implementation
## 5.4 Four-Lepton Results
## 5.5 Diphoton Results

# 6. Discovery Sensitivity
## 6.1 Toy Monte Carlo Experiments
## 6.2 Background-Only Hypthesis
## 6.3 Signal-Plus-Background Hypothesis
## 6.4 Discovery Test Statistic
## 6.5 $5\sigma$ Sensitivity
## 6.6 Four-Lepton Results
## 6.7 Diphoton Results

## References
ATLAS Collaboration (2020). ATLAS 13 TeV samples collection at least four leptons (electron or muon), for 2020 Open Data release. CERN Open Data Portal. DOI:10.7483/OPENDATA.ATLAS.2Y1T.TLGL

ATLAS Collaboration (2020). ATLAS 13 TeV samples collection Gamma-Gamma, for 2020 Open Data release. CERN Open Data Portal. DOI:10.7483/OPENDATA.ATLAS.B5BJ.3SGS

ATLAS Collaboration (2025). ATLAS 13 TeV Open Data: 2025 Beta Release. CERN Open Data Portal. https://opendata.atlas.cern/docs/category/13-tev-2025-beta-release
