"""
Generate Table 1: parameter uncertainty propagated in run_scenarios.py.

Loads the top-N parsets from the shrunk calibration
(results/gabon_calib.obj) and compares them to the calibration search
bounds defined in run_sims.make_calib_pars. Writes:
  - results/table1_parameter_uncertainty.csv
  - TABLE1.md
"""
import sciris as sc
import pandas as pd
import numpy as np

N_PARSETS = 10

# Priors from run_sims.make_calib_pars (v3), format [best_or_init, low, high].
# beta is fixed at 0.120 in make_sim (not calibrated on v3).
PRIORS = {
    'm_cross_layer':                   [0.15, 0.1,  0.7],
    'f_cross_layer':                   [0.1,  0.05, 0.5],
    'network_m_partners_casual':       [0.2,  0.1,  0.6],
    'network_f_partners_casual':       [0.2,  0.1,  0.6],
    'cross_immunity_rel_sev_loc':      [1.0,  0.5,  1.5],
    'hi5_cancer_fn_transform_prob':    [1.5e-3, 0.5e-3, 2.5e-3],
    'hi5_cin_fn_k':                    [0.15, 0.1,  0.25],
    'hi5_dur_cin_mean':                [4.5,  3.5,  5.5],
    'hi5_dur_cin_std':                 [20.0, 16.0, 24.0],
    'ohr_cancer_fn_transform_prob':    [1.5e-3, 0.5e-3, 2.5e-3],
    'ohr_cin_fn_k':                    [0.15, 0.1,  0.25],
    'ohr_dur_cin_mean':                [4.5,  3.5,  5.5],
    'ohr_dur_cin_std':                 [20.0, 16.0, 24.0],
}

# Flat key in the parset dict returned by rs.load_top_parsets (dotted).
KEY = {
    'm_cross_layer':                   'm_cross_layer',
    'f_cross_layer':                   'f_cross_layer',
    'network_m_partners_casual':       'network.m_partners_casual',
    'network_f_partners_casual':       'network.f_partners_casual',
    'cross_immunity_rel_sev_loc':      'cross_immunity.rel_sev.loc',
    'hi5_cancer_fn_transform_prob':    'hi5.cancer_fn.transform_prob',
    'hi5_cin_fn_k':                    'hi5.cin_fn.k',
    'hi5_dur_cin_mean':                'hi5.dur_cin.mean',
    'hi5_dur_cin_std':                 'hi5.dur_cin.std',
    'ohr_cancer_fn_transform_prob':    'ohr.cancer_fn.transform_prob',
    'ohr_cin_fn_k':                    'ohr.cin_fn.k',
    'ohr_dur_cin_mean':                'ohr.dur_cin.mean',
    'ohr_dur_cin_std':                 'ohr.dur_cin.std',
}


def fmt(x):
    """Human-readable number formatting for markdown."""
    if x is None:
        return '—'
    if abs(x) < 5e-3 and x != 0:
        return f'{x:.2e}'
    return f'{x:.3f}'


import run_sims as rs
parsets = rs.load_top_parsets(N_PARSETS)
# Load shrunk calib for the mismatch value only
calib = sc.load('results/gabon_calib.obj')

rows = []
for name, key in KEY.items():
    prior = PRIORS[name]
    values = np.array([p[key] for p in parsets])
    rows.append({
        'parameter':    name,
        'prior_low':    prior[1],
        'prior_high':   prior[2],
        'best':         float(values[0]),
        'top10_min':    float(values.min()),
        'top10_median': float(np.median(values)),
        'top10_max':    float(values.max()),
    })

df = pd.DataFrame(rows)
df.to_csv('results/table1_parameter_uncertainty.csv', index=False)

# Markdown
lines = ['| Parameter | prior_low | prior_high | best | top10_min | top10_median | top10_max |',
         '|-----------|-----------|------------|------|-----------|--------------|-----------|']
for _, r in df.iterrows():
    lines.append(
        f"| `{r['parameter']}` | {fmt(r['prior_low'])} | {fmt(r['prior_high'])} | "
        f"{fmt(r['best'])} | {fmt(r['top10_min'])} | {fmt(r['top10_median'])} | {fmt(r['top10_max'])} |"
    )
md = '\n'.join(lines)

best_mismatch = float(calib.df.iloc[0]['mismatch']) if hasattr(calib, 'df') else None
mismatch_str = f' (mismatch = {best_mismatch:.3f})' if best_mismatch is not None else ''

with open('TABLE1.md', 'w') as f:
    f.write('# Table 1 — parameter uncertainty propagated in `run_scenarios.py`\n\n')
    f.write(f'The top {N_PARSETS} calibration parsets (ranked by mismatch, taken from '
            f'the top-{N_PARSETS} rows of `results/gabon_calib.obj`) drive the parameter '
            f'uncertainty band in the scenario ensemble. Prior bounds are the calibration '
            f'search intervals defined in `run_sims.make_calib_pars`. **best** is rank 1 '
            f'by mismatch{mismatch_str}.\n\n')
    f.write(md + '\n')

print(md)
print(f'\nSaved TABLE1.md and results/table1_parameter_uncertainty.csv ({N_PARSETS} parsets)')