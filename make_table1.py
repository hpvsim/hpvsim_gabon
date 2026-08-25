"""
Generate Table 1: parameter uncertainty propagated in run_scenarios.py.

Loads the top-N calibration parsets (from results/gabon_pars_top50.obj) and
compares them to the calibration search bounds defined in run_sims.run_calib.
Writes:
  - results/table1_parameter_uncertainty.csv
  - TABLE1.md
"""
import sciris as sc
import pandas as pd
import numpy as np

N_PARSETS = 10

# Priors from run_sims.run_calib, format [initial, low, high, step]
PRIORS = {
    'beta':                          [0.2,   0.1,    0.34,   0.02],
    'm_cross_layer':                 [0.3,   0.1,    0.7,    0.05],
    'f_cross_layer':                 [0.1,   0.05,   0.5,    0.05],
    'm_partners_c_par1':             [0.2,   0.1,    0.6,    0.02],
    'f_partners_c_par1':             [0.2,   0.1,    0.6,    0.02],
    'sev_dist_par1':                 [1.0,   0.5,    1.5,    0.01],
    'hi5_cancer_fn_transform_prob':  [1.5e-3, 0.5e-3, 2.5e-3, 2e-4],
    'hi5_cin_fn_k':                  [0.15,  0.1,    0.25,   0.01],
    'hi5_dur_cin_par1':              [4.5,   3.5,    5.5,    0.5],
    'hi5_dur_cin_par2':              [20.0,  16.0,   24.0,   0.5],
    'ohr_cancer_fn_transform_prob':  [1.5e-3, 0.5e-3, 2.5e-3, 2e-4],
    'ohr_cin_fn_k':                  [0.15,  0.1,    0.25,   0.01],
    'ohr_dur_cin_par1':              [4.5,   3.5,    5.5,    0.5],
    'ohr_dur_cin_par2':              [20.0,  16.0,   24.0,   0.5],
}

# Path into the parset dict returned by calib.trial_pars_to_sim_pars(which_pars=i)
PATH = {
    'beta':                          ['beta'],
    'm_cross_layer':                 ['m_cross_layer'],
    'f_cross_layer':                 ['f_cross_layer'],
    'm_partners_c_par1':             ['m_partners', 'c', 'par1'],
    'f_partners_c_par1':             ['f_partners', 'c', 'par1'],
    'sev_dist_par1':                 ['sev_dist', 'par1'],
    'hi5_cancer_fn_transform_prob':  ['genotype_pars', 'hi5', 'cancer_fn', 'transform_prob'],
    'hi5_cin_fn_k':                  ['genotype_pars', 'hi5', 'cin_fn', 'k'],
    'hi5_dur_cin_par1':              ['genotype_pars', 'hi5', 'dur_cin', 'par1'],
    'hi5_dur_cin_par2':              ['genotype_pars', 'hi5', 'dur_cin', 'par2'],
    'ohr_cancer_fn_transform_prob':  ['genotype_pars', 'ohr', 'cancer_fn', 'transform_prob'],
    'ohr_cin_fn_k':                  ['genotype_pars', 'ohr', 'cin_fn', 'k'],
    'ohr_dur_cin_par1':              ['genotype_pars', 'ohr', 'dur_cin', 'par1'],
    'ohr_dur_cin_par2':              ['genotype_pars', 'ohr', 'dur_cin', 'par2'],
}


def get(parset, path):
    for k in path:
        parset = parset[k]
    return parset


def fmt(x):
    """Human-readable number formatting for markdown."""
    if x is None:
        return '—'
    if abs(x) < 5e-3 and x != 0:
        return f'{x:.2e}'
    return f'{x:.3f}'


parsets = sc.loadobj('results/gabon_pars_top50.obj')[:N_PARSETS]

rows = []
for name, path in PATH.items():
    prior = PRIORS[name]
    values = np.array([get(p, path) for p in parsets])
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

with open('TABLE1.md', 'w') as f:
    f.write('# Table 1 — parameter uncertainty propagated in `run_scenarios.py`\n\n')
    f.write(f'The top {N_PARSETS} calibration parsets (ranked by mismatch, taken from '
            f'`results/gabon_pars_top50.obj`) drive the parameter uncertainty band in the '
            f'scenario ensemble. Prior bounds are the calibration search intervals defined '
            f'in `run_sims.run_calib`. **best** is rank 1 by mismatch (mismatch = 0.875).\n\n')
    f.write(md + '\n')

print(md)
print(f'\nSaved TABLE1.md and results/table1_parameter_uncertainty.csv ({N_PARSETS} parsets)')