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
The analysis begins with publicly available proton-proton collision data from the ATLAS Collaboration through the CERN Open Data Portal. The datasets used here were collected during the 2016 LHC data-taking period at a center-of-mass energy of $13$ TeV and were released as part of the ATLAS 2020 Open Data release for educational use.

The Open Data release provides both real collision data and corresponding simulated Monte Carlo samples. The collision data represent the events observed by the ATLAS detector, while the simulated samples provide modeled examples of Standard Model background processes and selected signal processes. Keeping these two types of samples separate is important later in the analysis, where the observed data are compared with expected signal and background contributions.

This analysis uses two complementary ATLAS data collections:

Four-lepton channel: The ATLAS 13 TeV samples collection containing at least four leptons (electron or muon) provides events preselected to contain at least four electrons or muons. The collection consists of 111 files, totaling approximately 930.5 MiB.
Diphoton channel: The ATLAS 13 TeV $\gamma\gamma$ samples collection provides events preselected to contain at least two photons. The collection consists of 10 files, totaling approximately 3.5 GiB.

These collections are not yet the final analysis samples. Instead, they provide the starting point for the analysis workflow. Each collection has already undergone a loose ATLAS preselection at the object and event levels, reducing the number of events that need to be considered while retaining events potentially relevant to the corresponding Higgs decay channel.

The workflow used in this repository can therefore be viewed as:

$$
\text{ATLAS Open Data}
\rightarrow
\text{Dataset Preparation}
\rightarrow
\text{Analysis Selection}
\rightarrow
\text{Reconstructed Observable}
\rightarrow
\text{Statistical Analysis}.
$$

The following sections describe the two datasets in more detail, including the structure of their event-level variables and how those variables are used to construct the analysis samples.

Official ATLAS documentation for the datasets and their event-level properties is provided in the References. These resources are particularly useful when reproducing the analysis or adapting the workflow to a different dataset. 

## 2.2 Four-Lepton Dataset ($H \rightarrow ZZ^* \rightarrow 4 \ell$)
The four-lepton analysis uses the ATLAS 13 TeV samples collection containing at least four leptons (electron or muon) from the 2020 Open Data release. The collection contains both real collision data and simulated Monte Carlo samples representing relevant Standard Model backgrounds and Higgs signal events.

Understanding the Event Data

For the $H \rightarrow ZZ^* \rightarrow 4\ell$ analysis, we are interested in events containing four reconstructed leptons. These leptons can be electrons or muons, giving final states such as

$$
4e,\qquad 4\mu,\qquad 2e2\mu.
$$

The released event data contain information about the individual leptons that can be used to determine whether an event satisfies the additional analysis selections and to reconstruct the four-lepton invariant mass.

The primary lepton-level quantities used in this analysis include:

| Variable        | Meaning                              | Role in the analysis                                       |
| :-------------- | :----------------------------------- | :--------------------------------------------------------- |
| `lep_pt`        | Transverse momentum $p_T$            | Used for lepton momentum requirements                      |
| `lep_eta`       | Pseudorapidity $\eta$                | Describes the lepton's direction relative to the beam axis |
| `lep_phi`       | Azimuthal angle $\phi$               | Used to determine angular separations between leptons      |
| `lep_E`         | Energy $E$                           | Used when reconstructing the invariant mass                |
| `lep_charge`    | Electric charge                      | Used to identify opposite-sign lepton pairs                |
| `lep_type`      | Lepton type                          | Distinguishes electrons from muons                         |
| `lep_isTightID` | Identification requirement           | Used to select well-identified leptons                     |
| `lep_ptcone30`  | Track-based isolation variable       | Used to characterize nearby activity around a lepton       |
| `lep_etcone20`  | Calorimeter-based isolation variable | Used to characterize nearby energy deposition              |


The notation used in the released data is also reflected in the code. For example, lep_pt_1 through lep_pt_4 refer to the transverse momenta of the four leptons assigned to an event, while the corresponding _eta, _phi, and _E variables provide their other kinematic properties.

The transverse momentum is the component of a particle's momentum perpendicular to the proton-proton beam axis:

$$
p_T = \sqrt{p_x^2+p_y^2}.
$$

The pseudorapidity is defined as

$$
\eta = -\ln\left(\tan\frac{\theta}{2}\right),
$$

where $\theta$ is the polar angle measured relative to the beam axis. The azimuthal angle $\phi$ describes the direction of the particle in the plane transverse to the beam.

These quantities are sufficient to construct several useful observables. In particular, the angular separation between two leptons is calculated using

$$
\Delta R =
\sqrt{(\Delta\eta)^2+(\Delta\phi)^2}.
$$

This quantity is later used as part of the event-selection procedure.

From the Released Sample to the Analysis Sample

The ATLAS collection has already been preselected to contain at least four leptons. We then apply additional analysis-level requirements to the events. These requirements use the kinematic, identification, isolation, and angular information provided in the dataset.

The resulting workflow is:

$$
\text{ATLAS 4$\ell$ sample}
\rightarrow
\text{Lepton selection}
\rightarrow
\text{Event selection}
\rightarrow
\text{Four-lepton system}
\rightarrow
m_{4\ell}.
$$

After the selection requirements are applied, the four selected leptons are combined to reconstruct the invariant mass of the four-lepton system. This distribution is the primary observable used to search for the Higgs signal.

The real collision data provide the observed distribution, while the simulated signal and background samples provide the components needed to construct the statistical model.

The specific selection criteria and their implementation are described in Section 3: Event Selection.

## 2.3 Diphoton Dataset ($H \rightarrow \gamma\gamma$)
The diphoton analysis uses the ATLAS 13 TeV $\gamma\gamma$ samples collection from the 2020 Open Data release. Like the four-lepton collection, it contains real collision data together with Monte Carlo simulated samples representing relevant background and signal processes.

Understanding the Event Data

For the $H \rightarrow \gamma\gamma$ analysis, we select events containing at least two reconstructed photons. The relevant photon-level information includes their kinematic properties together with identification and isolation variables.

The primary quantities used in this analysis include:
| Variable                        | Meaning                                          | Role in the analysis                                   |
| :------------------------------ | :----------------------------------------------- | :----------------------------------------------------- |
| `photon_pt`                     | Transverse momentum (`pT`)                       | Used for photon momentum requirements                  |
| `photon_eta`                    | Pseudorapidity (`η`)                             | Describes the photon's direction                       |
| `photon_phi`                    | Azimuthal angle (`ϕ`)                            | Used to determine angular relationships                |
| `photon_E`                      | Energy (`E`)                                     | Used when reconstructing the invariant mass            |
| Photon identification variables | Photon reconstruction/identification information | Used to select suitable photon candidates              |
| Photon isolation variables      | Activity surrounding the photon                  | Used to reduce contamination from non-isolated objects |


As with the four-lepton dataset, the released variables provide the information required to move from the preselected events to the final analysis sample.

The photon transverse momentum and angular variables are used to impose additional selection requirements, while the photon energies and directions are used to reconstruct the invariant mass of the diphoton system.

For two photons with four-momenta $p_1$ and $p_2$, the diphoton invariant mass is obtained from

$$
m_{\gamma\gamma}^2 = (p_1+p_2)^2.
$$

The resulting $m_{\gamma\gamma}$ distribution provides the observable in which a Higgs boson contribution can appear as an excess above the smoothly varying background.

From the Released Sample to the Analysis Sample

The ATLAS $\gamma\gamma$ collection has already undergone a loose preselection requiring at least two photons. We then apply additional analysis-level requirements to isolate events consistent with the $H \rightarrow \gamma\gamma$ decay.

The workflow is therefore:

$$
\text{ATLAS }\gamma\gamma\text{ sample}
\rightarrow
\text{Photon selection}
\rightarrow
\text{Event selection}
\rightarrow
\text{Diphoton system}
\rightarrow
m_{\gamma\gamma}.
$$

The selected collision data are then used to construct the observed diphoton invariant-mass distribution. Simulated samples provide information about the expected signal and background contributions where applicable.

As with the four-lepton channel, the goal is not simply to obtain a final histogram, but to establish a sequence of reproducible transformations from the public ATLAS dataset to the observable used in the statistical analysis.

The specific selection criteria and their implementation are described in Section 3: Event Selection.

Dataset Documentation

When reproducing or modifying this analysis, the official ATLAS documentation should be used to verify the meaning, units, and representation of the released event variables. The documentation for both the four-lepton and diphoton collections is provided in the References.

This is particularly important when adapting the workflow to another ATLAS Open Data release or another analysis channel, since the available variables and their definitions may differ between datasets.


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
