import uproot
import awkward as ak
import vector
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Register vector behavior
vector.register_awkward()

# CONFIGURATION
LUMI = 10.0 # fb^-1
from pathlib import Path

base_path = Path(__file__).resolve().parent.parent / "RawData"

# Cross Section Overrides (pb)
XSEC_OVERRIDES = {
    'ZH': 0.00012,
    'WH': 0.00015
}

samples = {
    'data': {
        'files': [
            str(base_path / 'data_A.4lep.root'),
            str(base_path / 'data_B.4lep.root'),
            str(base_path / 'data_C.4lep.root'),
            str(base_path / 'data_D.4lep.root')
        ],
        'type': 'data', 'label': 'Data'
    },
    'background': {
        'files': [str(base_path / 'mc_363490.llll.4lep.root')],
        'type': 'mc', 'label': r'ZZ $\rightarrow 4\ell$', 'color': 'royalblue'
    },
    # Signals
    'signal_ggF': {
        'files': [str(base_path / 'mc_345060.ggH125_ZZ4lep.4lep.root')],
        'type': 'mc', 'label': 'ggF H(125)', 'color': 'red'
    },
    'signal_VBF': {
        'files': [str(base_path / 'mc_344235.VBFH125_ZZ4lep.4lep.root')],
        'type': 'mc', 'label': 'VBF H(125)', 'color': 'salmon' 
    },
    'signal_ZH': {
        'files': [str(base_path / 'mc_341947.ZH125_ZZ4lep.4lep.root')],
        'type': 'mc', 'label': 'ZH(125)', 'color': 'orange'
    },
    'signal_WH': {
        'files': [str(base_path / 'mc_341964.WH125_ZZ4lep.4lep.root')],
        'type': 'mc', 'label': 'WH(125)', 'color': 'coral' 
    }
}

branches = [
    'lep_pt', 'lep_eta', 'lep_phi', 'lep_E', 'lep_charge', 'lep_type', 
    'lep_isTightID', 'lep_ptcone30', 'lep_etcone20', 'lep_n',
    'lep_tracksigd0pvunbiased', 'lep_z0', 
    'trigE', 'trigM', 'scaleFactor_PILEUP', 'scaleFactor_ELE', 
    'scaleFactor_MUON', 'scaleFactor_LepTRIGGER', 
    'mcWeight', 'XSection', 'SumWeights', 'runNumber', 'eventNumber'
]

print("Configuration updated.")


def get_data(files, branches):
    """Lazily reads data from a list of files."""
    data_arrays = []
    for file in files:
        try:
            with uproot.open(f"{file}:mini") as tree:
                data_arrays.append(tree.arrays(branches, library="ak"))
        except FileNotFoundError:
            print(f"Warning: File {file} not found. Skipping.")
    return ak.concatenate(data_arrays) if data_arrays else None

def calculate_weight(events, luminosity, sample_name=""):
    """
    Calculates event weights using standard formula + XSection overrides.
    """
    # 1. Determine Cross Section
    xsec = events.XSection
    
    # Check for overrides (e.g. ZH or WH) based on sample name
    for key, val in XSEC_OVERRIDES.items():
        if key in sample_name: # e.g. if 'ZH' in 'signal_ZH'
            xsec = val
            
    # 2. Calculate Base Weight (convert xsec pb -> fb using *1000)
    weight = (luminosity * xsec * 1000 * events.mcWeight) / events.SumWeights
    
    # 3. Apply Scale Factors
    weight = weight * events.scaleFactor_PILEUP * \
             events.scaleFactor_ELE * \
             events.scaleFactor_MUON * \
             events.scaleFactor_LepTRIGGER
             
    return weight

print("Helper functions defined.")

results = {}
cutflow = {}

# Constants
M_Z = 91.2  # Z boson mass in GeV

print("Starting analysis loop")

for sample_name, config in samples.items():
    print(f"Processing {config['label']}...")
    
    # 1. LOAD DATA
    events = get_data(config['files'], branches)
    if events is None: continue
    
    cutflow[sample_name] = {'Initial': len(events)}
    
    # 2. TRIGGER
    trigger_mask = events.trigE | events.trigM
    events = events[trigger_mask]
    cutflow[sample_name]['Trigger'] = len(events)

    # 3. OBJECT SELECTION
    leptons = ak.zip({
        "pt": events.lep_pt,
        "eta": events.lep_eta,
        "phi": events.lep_phi,
        "E": events.lep_E,
        "charge": events.lep_charge,
        "d0sig": events.lep_tracksigd0pvunbiased,
        "z0": events.lep_z0,
        "ptcone30": events.lep_ptcone30,
        "etcone20": events.lep_etcone20
    }, with_name="Momentum4D")

    lep_type = events.lep_type
    
    # Calculate z0 * sin(theta)
    theta = 2 * np.arctan(np.exp(-leptons.eta))
    z0_sin_theta = np.abs(leptons.z0 * np.sin(theta))

    # --- LEPTON CUTS (vectorized) ---
    is_electron = np.abs(lep_type) == 11
    ele_crack_veto = (np.abs(leptons.eta) > 1.37) & (np.abs(leptons.eta) < 1.52)
    ele_cut = is_electron & (leptons.pt > 7000) & (np.abs(leptons.eta) < 2.47) & \
              (np.abs(leptons.d0sig) < 6.5) & (~ele_crack_veto)

    is_muon = np.abs(lep_type) == 13
    muon_cut = is_muon & (leptons.pt > 6000) & (np.abs(leptons.eta) < 2.7) & \
               (np.abs(leptons.d0sig) < 3.5)

    iso_cut = ((leptons.ptcone30 / leptons.pt) < 0.15) & \
              ((leptons.etcone20 / leptons.pt) < 0.20)
    common_cut = (z0_sin_theta < 0.5) & iso_cut

    good_lep_mask = (ele_cut | muon_cut) & common_cut
    events["leptons"] = leptons[good_lep_mask]
    events["lep_type_good"] = lep_type[good_lep_mask]
    
    cutflow[sample_name]['Lepton Object'] = len(events)

    # 4. EVENT SELECTION
    
    # Cut: EXACTLY 4 good leptons
    has_four = ak.num(events.leptons) == 4
    events = events[has_four]
    cutflow[sample_name]['== 4 Leptons'] = len(events)
    
    # Cut: Total Charge == 0
    total_charge = ak.sum(events.leptons.charge, axis=1)
    events = events[total_charge == 0]
    cutflow[sample_name]['Total Charge 0'] = len(events)
    
    # VECTORIZED OSSF Z CANDIDATE FINDING
    # Get the 4 leptons as separate arrays for vectorized ops
    l0 = events.leptons[:, 0]
    l1 = events.leptons[:, 1]
    l2 = events.leptons[:, 2]
    l3 = events.leptons[:, 3]
    
    t0 = events.lep_type_good[:, 0]
    t1 = events.lep_type_good[:, 1]
    t2 = events.lep_type_good[:, 2]
    t3 = events.lep_type_good[:, 3]
    
    c0 = l0.charge
    c1 = l1.charge
    c2 = l2.charge
    c3 = l3.charge
    
    # Calculate all 6 pair masses (in GeV)
    m01 = (l0 + l1).mass / 1000.0
    m02 = (l0 + l2).mass / 1000.0
    m03 = (l0 + l3).mass / 1000.0
    m12 = (l1 + l2).mass / 1000.0
    m13 = (l1 + l3).mass / 1000.0
    m23 = (l2 + l3).mass / 1000.0
    
    # Check OSSF for each pair (same flavor AND opposite sign)
    ossf_01 = (t0 == t1) & (c0 * c1 < 0)
    ossf_02 = (t0 == t2) & (c0 * c2 < 0)
    ossf_03 = (t0 == t3) & (c0 * c3 < 0)
    ossf_12 = (t1 == t2) & (c1 * c2 < 0)
    ossf_13 = (t1 == t3) & (c1 * c3 < 0)
    ossf_23 = (t2 == t3) & (c2 * c3 < 0)
    
    # Three possible non-overlapping pairings: (01,23), (02,13), (03,12)
    # For each pairing, check if BOTH pairs are OSSF
    valid_A = ossf_01 & ossf_23  # pairing (0,1) + (2,3)
    valid_B = ossf_02 & ossf_13  # pairing (0,2) + (1,3)
    valid_C = ossf_03 & ossf_12  # pairing (0,3) + (1,2)
    
    # At least one valid pairing must exist
    has_ossf = valid_A | valid_B | valid_C
    
    # For valid pairings, determine Z1 (closer to M_Z) and Z2
    # Initialize with NaN
    n_evt = len(events)
    mZ1 = np.full(n_evt, np.nan)
    mZ2 = np.full(n_evt, np.nan)
    
    # Convert to numpy for conditional assignment
    m01_np = ak.to_numpy(m01)
    m02_np = ak.to_numpy(m02)
    m03_np = ak.to_numpy(m03)
    m12_np = ak.to_numpy(m12)
    m13_np = ak.to_numpy(m13)
    m23_np = ak.to_numpy(m23)
    
    valid_A_np = ak.to_numpy(valid_A)
    valid_B_np = ak.to_numpy(valid_B)
    valid_C_np = ak.to_numpy(valid_C)
    
    # For each event, find the best pairing (minimize |mZ1 - M_Z|)
    for i in range(n_evt):
        candidates = []
        if valid_A_np[i]:
            ma, mb = m01_np[i], m23_np[i]
            if abs(ma - M_Z) < abs(mb - M_Z):
                candidates.append((ma, mb, abs(ma - M_Z)))
            else:
                candidates.append((mb, ma, abs(mb - M_Z)))
        if valid_B_np[i]:
            ma, mb = m02_np[i], m13_np[i]
            if abs(ma - M_Z) < abs(mb - M_Z):
                candidates.append((ma, mb, abs(ma - M_Z)))
            else:
                candidates.append((mb, ma, abs(mb - M_Z)))
        if valid_C_np[i]:
            ma, mb = m03_np[i], m12_np[i]
            if abs(ma - M_Z) < abs(mb - M_Z):
                candidates.append((ma, mb, abs(ma - M_Z)))
            else:
                candidates.append((mb, ma, abs(mb - M_Z)))
        
        if candidates:
            # Pick the one with Z1 closest to M_Z
            best = min(candidates, key=lambda x: x[2])
            mZ1[i] = best[0]
            mZ2[i] = best[1]
    
    # Apply OSSF cut
    events = events[has_ossf]
    mZ1 = mZ1[ak.to_numpy(has_ossf)]
    mZ2 = mZ2[ak.to_numpy(has_ossf)]
    cutflow[sample_name]['2 OSSF Pairs'] = len(events)
    
    # --- VECTORIZED Delta R Separation ---
    l0 = events.leptons[:, 0]
    l1 = events.leptons[:, 1]
    l2 = events.leptons[:, 2]
    l3 = events.leptons[:, 3]
    
    t0 = events.lep_type_good[:, 0]
    t1 = events.lep_type_good[:, 1]
    t2 = events.lep_type_good[:, 2]
    t3 = events.lep_type_good[:, 3]
    
    # Calculate all deltaR values
    dr01 = l0.deltaR(l1)
    dr02 = l0.deltaR(l2)
    dr03 = l0.deltaR(l3)
    dr12 = l1.deltaR(l2)
    dr13 = l1.deltaR(l3)
    dr23 = l2.deltaR(l3)
    
    # Same flavor: ΔR > 0.1, Different flavor: ΔR > 0.2
    pass_dr01 = ak.where(t0 == t1, dr01 > 0.1, dr01 > 0.2)
    pass_dr02 = ak.where(t0 == t2, dr02 > 0.1, dr02 > 0.2)
    pass_dr03 = ak.where(t0 == t3, dr03 > 0.1, dr03 > 0.2)
    pass_dr12 = ak.where(t1 == t2, dr12 > 0.1, dr12 > 0.2)
    pass_dr13 = ak.where(t1 == t3, dr13 > 0.1, dr13 > 0.2)
    pass_dr23 = ak.where(t2 == t3, dr23 > 0.1, dr23 > 0.2)
    
    pass_dr = pass_dr01 & pass_dr02 & pass_dr03 & pass_dr12 & pass_dr13 & pass_dr23
    
    events = events[pass_dr]
    mZ1 = mZ1[ak.to_numpy(pass_dr)]
    mZ2 = mZ2[ak.to_numpy(pass_dr)]
    cutflow[sample_name]['Lepton Separation'] = len(events)
    
    # --- VECTORIZED J/Psi Veto ---
    l0 = events.leptons[:, 0]
    l1 = events.leptons[:, 1]
    l2 = events.leptons[:, 2]
    l3 = events.leptons[:, 3]
    
    t0 = events.lep_type_good[:, 0]
    t1 = events.lep_type_good[:, 1]
    t2 = events.lep_type_good[:, 2]
    t3 = events.lep_type_good[:, 3]
    
    c0 = l0.charge
    c1 = l1.charge
    c2 = l2.charge
    c3 = l3.charge
    
    # Recalculate masses for remaining events
    m01 = (l0 + l1).mass / 1000.0
    m02 = (l0 + l2).mass / 1000.0
    m03 = (l0 + l3).mass / 1000.0
    m12 = (l1 + l2).mass / 1000.0
    m13 = (l1 + l3).mass / 1000.0
    m23 = (l2 + l3).mass / 1000.0
    
    # OSSF flags
    ossf_01 = (t0 == t1) & (c0 * c1 < 0)
    ossf_02 = (t0 == t2) & (c0 * c2 < 0)
    ossf_03 = (t0 == t3) & (c0 * c3 < 0)
    ossf_12 = (t1 == t2) & (c1 * c2 < 0)
    ossf_13 = (t1 == t3) & (c1 * c3 < 0)
    ossf_23 = (t2 == t3) & (c2 * c3 < 0)
    
    # J/psi veto: no OSSF pair with mass <= 5 GeV
    jpsi_fail_01 = ossf_01 & (m01 <= 5.0)
    jpsi_fail_02 = ossf_02 & (m02 <= 5.0)
    jpsi_fail_03 = ossf_03 & (m03 <= 5.0)
    jpsi_fail_12 = ossf_12 & (m12 <= 5.0)
    jpsi_fail_13 = ossf_13 & (m13 <= 5.0)
    jpsi_fail_23 = ossf_23 & (m23 <= 5.0)
    
    pass_jpsi = ~(jpsi_fail_01 | jpsi_fail_02 | jpsi_fail_03 | 
                  jpsi_fail_12 | jpsi_fail_13 | jpsi_fail_23)
    
    events = events[pass_jpsi]
    mZ1 = mZ1[ak.to_numpy(pass_jpsi)]
    mZ2 = mZ2[ak.to_numpy(pass_jpsi)]
    cutflow[sample_name]['J/Psi Veto'] = len(events)
    
    # --- Ranked pT Thresholds ---
    pt = events.leptons.pt
    pt_sorted = ak.sort(pt, axis=1, ascending=False)
    pass_pt = (pt_sorted[:,0] > 20000) & (pt_sorted[:,1] > 15000) & (pt_sorted[:,2] > 10000)
    events = events[pass_pt]
    mZ1 = mZ1[ak.to_numpy(pass_pt)]
    mZ2 = mZ2[ak.to_numpy(pass_pt)]
    cutflow[sample_name]['Ranked pT'] = len(events)
    
    # --- Z1 Mass Window ---
    pass_z1 = (mZ1 > 50) & (mZ1 < 106)
    events = events[pass_z1]
    mZ1 = mZ1[pass_z1]
    mZ2 = mZ2[pass_z1]
    cutflow[sample_name]['Z1 Mass'] = len(events)
    
    # --- Z2 Mass Window ---
    pass_z2 = (mZ2 > 12) & (mZ2 < 115)
    events = events[pass_z2]
    mZ1 = mZ1[pass_z2]
    mZ2 = mZ2[pass_z2]
    cutflow[sample_name]['Z2 Mass'] = len(events)
    
    # 5. FINAL MASS & WEIGHT
    p4_sum = events.leptons[:,0] + events.leptons[:,1] + events.leptons[:,2] + events.leptons[:,3]
    events["m4l"] = p4_sum.mass / 1000.0
    
    if config['type'] == 'mc':
        events["final_weight"] = calculate_weight(events, LUMI, sample_name)
    else:
        events["final_weight"] = np.ones(len(events))
        
    results[sample_name] = {
        'mass': events.m4l,
        'weight': events.final_weight,
        'color': config.get('color', 'black'),
        'label': config.get('label', 'Data'),
        'mZ1': mZ1,
        'mZ2': mZ2
    }
    
    print(f"  -> {len(events)} events pass all cuts")

print("\nAnalysis loop complete!")

# 5. CUTFLOW TABLE (Raw Counts vs Weighted Yields)

if 'data' in cutflow:
    steps = list(cutflow['data'].keys())
    signal_keys = [k for k in cutflow.keys() if k.startswith('signal_')]
    
    # --- Calculate average weights for each sample (from final events) ---
    # These will be used to estimate weighted yields at each cut step
    avg_weights = {}
    
    # Data: weight = 1
    avg_weights['data'] = 1.0
    
    # Background
    if 'background' in results:
        b_raw_final = cutflow['background']['Z2 Mass']
        b_weighted_final = float(np.sum(results['background']['weight']))
        avg_weights['background'] = b_weighted_final / b_raw_final if b_raw_final > 0 else 0
    
    # Signals (combined average)
    total_sig_raw = sum(cutflow[k]['Z2 Mass'] for k in signal_keys if k in cutflow)
    total_sig_weighted = sum(float(np.sum(results[k]['weight'])) for k in signal_keys if k in results)
    avg_weights['signal'] = total_sig_weighted / total_sig_raw if total_sig_raw > 0 else 0
    
    # --- MAIN CUTFLOW TABLE ---
    print(f"CUTFLOW TABLE: Raw Counts vs Weighted Yields (L = {LUMI} fb⁻¹)")

    print(f"{'Cut':<20} | {'Data':<10} | {'Signal':^25} | {'Background':^25}")
    print(f"{'':<20} | {'(wt=1)':<10} | {'Raw':<10} {'Weighted':<15} | {'Raw':<10} {'Weighted':<15}")
    print("-" * 110)
    
    for step in steps:
        d_raw = cutflow['data'].get(step, 0)
        s_raw = sum(cutflow[k].get(step, 0) for k in signal_keys)
        b_raw = cutflow.get('background', {}).get(step, 0)
        
        # Weighted yields (estimated using final average weights)
        d_wt = d_raw * avg_weights['data']
        s_wt = s_raw * avg_weights['signal']
        b_wt = b_raw * avg_weights['background']
        
        print(f"{step:<20} | {d_raw:<10} | {s_raw:<10} {s_wt:<15.2f} | {b_raw:<10} {b_wt:<15.2f}")
    
    print("-" * 110)
    
    # --- FINAL YIELDS (Exact values using actual weights) ---
    print("\n")

    print("FINAL YIELDS AFTER ALL CUTS (Exact weighted values)")
    print(f"{'Sample':<25} | {'Raw Events':<15} | {'Weighted Yield':<15}")

    # Data
    d_final_raw = cutflow['data']['Z2 Mass']
    print(f"{'Data':<25} | {d_final_raw:<15} | {d_final_raw:<15.2f}")
    
    # Background
    b_final_raw = cutflow['background']['Z2 Mass']
    b_final_wt = float(np.sum(results['background']['weight']))
    print(f"{'Background (ZZ→4ℓ)':<25} | {b_final_raw:<15} | {b_final_wt:<15.2f}")
    
    # Signal breakdown
    print("-" * 80)
    total_sig_raw = 0
    total_sig_wt = 0.0
    for k in signal_keys:
        if k in results:
            s_raw = cutflow[k]['Z2 Mass']
            s_wt = float(np.sum(results[k]['weight']))
            total_sig_raw += s_raw
            total_sig_wt += s_wt
            print(f"  {results[k]['label']:<23} | {s_raw:<15} | {s_wt:<15.4f}")
    
    print("-" * 80)
    print(f"{'TOTAL SIGNAL':<25} | {total_sig_raw:<15} | {total_sig_wt:<15.4f}")
    print(f"{'TOTAL MC (Sig+Bkg)':<25} | {total_sig_raw + b_final_raw:<15} | {total_sig_wt + b_final_wt:<15.2f}")

# Combine all signal masses and weights
signal_masses = np.concatenate([results[k]['mass'] for k in results if k.startswith('signal_')])
signal_weights = np.concatenate([results[k]['weight'] for k in results if k.startswith('signal_')])

# --- 1. FULL SPECTRUM HISTOGRAM DATA (80 to 250, bin width 5) ---
full_bins = np.arange(80, 255, 5)
full_bin_min = full_bins[:-1]
full_bin_max = full_bins[1:]
full_bin_center = 0.5 * (full_bin_min + full_bin_max)
full_bin_width = np.full_like(full_bin_min, 5.0)

full_data_obs, _ = np.histogram(results['data']['mass'], bins=full_bins)
full_sig_exp, _ = np.histogram(signal_masses, bins=full_bins, weights=signal_weights)
full_bkg_exp, _ = np.histogram(results['background']['mass'], bins=full_bins, weights=results['background']['weight'])
full_data_err = np.sqrt(full_data_obs)

df_full = pd.DataFrame({
    'bin_min': full_bin_min,
    'bin_max': full_bin_max,
    'bin_center': full_bin_center,
    'bin_width': full_bin_width,
    'data_observed': full_data_obs,
    'data_stat_error': full_data_err,
    'signal_expected': full_sig_exp,
    'background_expected': full_bkg_exp
})
df_full.to_csv('full_spectrum_histogram_data.csv', index=False)

# --- 2. HIGGS FIT HISTOGRAM DATA (110 to 140, bin width derived automatically) ---
# To create integer-friendly bins from 110 to 140, width = 1 GeV is used (30 bins)
fit_bins = np.arange(110, 141, 1)
fit_bin_min = fit_bins[:-1]
fit_bin_max = fit_bins[1:]
fit_bin_center = 0.5 * (fit_bin_min + fit_bin_max)
fit_bin_width = fit_bin_max - fit_bin_min

fit_data_obs, _ = np.histogram(results['data']['mass'], bins=fit_bins)
fit_sig_exp, _ = np.histogram(signal_masses, bins=fit_bins, weights=signal_weights)
fit_bkg_exp, _ = np.histogram(results['background']['mass'], bins=fit_bins, weights=results['background']['weight'])
fit_data_err = np.sqrt(fit_data_obs)

df_fit = pd.DataFrame({
    'bin_min': fit_bin_min,
    'bin_max': fit_bin_max,
    'bin_center': fit_bin_center,
    'bin_width': fit_bin_width,
    'data_observed': fit_data_obs,
    'data_stat_error': fit_data_err,
    'signal_expected': fit_sig_exp,
    'background_expected': fit_bkg_exp
})
df_fit.to_csv('higgs_fit_histogram_data.csv', index=False)

# --- 3. BACKGROUND TEMPLATE (111 to 139) ---
# Bin centers starting at 111 up to 139 (step of 1)
bkg_centers = np.arange(111, 140, 1)
# Bins corresponding to integer centers 111-139 -> edges 110.5 to 139.5
bkg_template_bins = np.arange(110.5, 140.5, 1)
bkg_counts, _ = np.histogram(
    results['background']['mass'], 
    bins=bkg_template_bins, 
    weights=results['background']['weight']
)

df_bkg_template = pd.DataFrame({
    'bin_center': bkg_centers,
    'background': bkg_counts
})
df_bkg_template.to_csv('background_template.csv', index=False)

# --- 4. PYTHON BEST FIT / SIGNIFICANCE CALCULATIONS ---
# Simple Poisson profile likelihood calculation on fit region (110-140 GeV)
n_obs = np.sum(fit_data_obs)
b_exp = np.sum(fit_bkg_exp)
s_exp = np.sum(fit_sig_exp)

# Best fit signal strength mu = (N_obs - B) / S
best_fit_mu = (n_obs - b_exp) / s_exp if s_exp > 0 else 0.0
mu_error = np.sqrt(n_obs) / s_exp if s_exp > 0 else 0.0

# Asymptotic Significance Z = sqrt(2 * ((S+B)*ln(1 + S/B) - S))
if b_exp > 0 and s_exp > 0:
    significance = np.sqrt(2 * ((s_exp + b_exp) * np.log(1 + s_exp / b_exp) - s_exp))
else:
    significance = 0.0

# --- 5. METADATA CSV ---
df_metadata = pd.DataFrame([{
    'fit_min': 110,
    'fit_max': 140,
    'n_bins_fit': len(fit_bin_center),
    'bin_width_fit': fit_bin_width[0],
    'full_min': 80,
    'full_max': 250,
    'n_bins_full': len(full_bin_center),
    'bin_width_full': 5,
    'luminosity_fb': LUMI,
    'signal_total_fit': float(np.sum(fit_sig_exp)),
    'background_total_fit': float(np.sum(fit_bkg_exp)),
    'data_total_fit': int(np.sum(fit_data_obs)),
    'signal_total_full': float(np.sum(full_sig_exp)),
    'background_total_full': float(np.sum(full_bkg_exp)),
    'data_total_full': int(np.sum(full_data_obs)),
    'python_best_fit_mu': float(best_fit_mu),
    'python_mu_error': float(mu_error),
    'python_significance': float(significance)
}])

df_metadata.to_csv('metadata.csv', index=False)


