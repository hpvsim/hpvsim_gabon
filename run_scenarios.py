'''
Run HPVsim scenarios (v3.1.0). Requires an HPC/VM to run with debug=False.
'''


# %% General settings

import os

os.environ.update(
    OMP_NUM_THREADS='1',
    OPENBLAS_NUM_THREADS='1',
    NUMEXPR_NUM_THREADS='1',
    MKL_NUM_THREADS='1',
)

# Standard imports
import numpy as np
import sciris as sc
import starsim as ss
import hpvsim as hpv

# Imports from this repository
import run_sims as rs


# What to run
debug = 0
n_parsets = [10, 1][debug]  # Top-N calibration parsets, for propagating parameter uncertainty
n_seeds = [1, 1][debug]     # Stochastic seeds per parset per scenario
end_year = 2100  # Simulation horizon; must match the vaccination scale-up schedule below

# Screening eligibility age window (30-50 = 20 years). See spec §3.2.
# Ported from hpvsim_pxv_younger's _annual_from_lifetime: interprets
# screen_coverage as lifetime coverage; per-year prob is
# p = 1 - (1-C)^(1/N) with N = full age range = 20.
SCREEN_AGE_LO = 30
SCREEN_AGE_HI = 50
SCREEN_AGE_YEARS = SCREEN_AGE_HI - SCREEN_AGE_LO  # 20


def _annual_from_lifetime(lifetime_cov, n_years=SCREEN_AGE_YEARS):
    """Convert lifetime screening coverage to per-year prob for routine_screening."""
    c = float(min(max(lifetime_cov, 0.0), 1.0))
    if c <= 0:
        return 0.0
    if c >= 1.0:
        return 1.0
    return 1.0 - (1.0 - c) ** (1.0 / n_years)


# %% Functions
def make_screen_treat(screen_coverage=0.15, triage_coverage=0.9, treat_coverage=0.75, start_year=2020):
    """Screening + treatment intervention list. `screen_coverage` is lifetime
    coverage over ages SCREEN_AGE_LO to SCREEN_AGE_HI (see spec §3.2)."""

    age_range = [SCREEN_AGE_LO, SCREEN_AGE_HI]
    model_annual_screen_prob = _annual_from_lifetime(screen_coverage)

    screening = hpv.routine_screening(
        prob=model_annual_screen_prob,
        start_year=start_year,
        product='hpv',
        age_range=age_range,
        label='screening',
        name='intv_screening',
    )

    screen_positive = lambda sim: sim.interventions['intv_screening'].outcomes['positive']
    assign_treatment = hpv.routine_triage(
        start_year=start_year,
        prob=triage_coverage,
        annual_prob=False,
        product='tx_assigner',
        eligibility=screen_positive,
        label='tx assigner',
        name='intv_tx_assigner',
    )

    ablation_eligible = lambda sim: sim.interventions['intv_tx_assigner'].outcomes['ablation']
    ablation = hpv.treat_num(
        prob=treat_coverage,
        product='ablation',
        eligibility=ablation_eligible,
        label='ablation',
        name='intv_ablation',
    )

    excision_eligible = lambda sim: np.union1d(
        sim.interventions['intv_tx_assigner'].outcomes['excision'],
        sim.interventions['intv_ablation'].outcomes['unsuccessful'],
    )
    excision = hpv.treat_num(
        prob=treat_coverage,
        product='excision',
        eligibility=excision_eligible,
        label='excision',
        name='intv_excision',
    )

    radiation_eligible = lambda sim: sim.interventions['intv_tx_assigner'].outcomes['radiation']
    radiation = hpv.treat_num(
        prob=treat_coverage / 4,  # assume 4x dropoff in coverage for cancer treatment vs pre-cancer treatment
        product=hpv.radiation(),
        eligibility=radiation_eligible,
        label='radiation',
        name='intv_radiation',
    )

    return [screening, assign_treatment, ablation, excision, radiation]


def make_screen_treat_scenarios():
    """Make screening & treatment scenarios, looping over screening coverage."""
    screen_treat_scenarios = dict()
    screen_coverages = [0.1, 0.4, 0.9]  # Low / medium / high lifetime screening coverage
    for scov in screen_coverages:
        label = f'Screen {int(scov*100)}%'
        screen_treat_scenarios[label] = make_screen_treat(screen_coverage=scov)
    return screen_treat_scenarios


def _ever_vaxed_uids(sim):
    """UIDs vaxed by any HPV vaccine intervention so far."""
    out = ss.uids()
    for intv in sim.interventions.values():
        v = getattr(intv, 'vaccinated', None)
        if v is None:
            continue
        out = out.union(v.uids)
    return out


def _never_vaxed(sim):
    """Eligibility: alive AND not yet vaxed by any HPV vaccine intervention."""
    return sim.people.alive.uids.remove(_ever_vaxed_uids(sim))


def make_vx_scenarios(product='bivalent', start_year=2025, end_year=end_year):
    """Vaccination scenarios: no vaccination vs a scale-up to 90% routine coverage."""

    routine_age = (9, 10)

    vx_scenarios = dict()
    vx_scenarios['No vaccination'] = []

    vx_years = np.arange(start_year, end_year + 1)
    scaleup = [0.3, 0.6, 0.9]
    final_cov = 0.9
    vx_cov = np.concatenate([scaleup + [final_cov] * (len(vx_years) - len(scaleup))])

    routine_vx = hpv.campaign_vx(
        prob=vx_cov,
        years=vx_years,
        product=product,
        age_range=routine_age,
        eligibility=_never_vaxed,
        label='Routine vx',
        name='routine_vx',
    )
    vx_scenarios['90% vax coverage'] = [routine_vx]

    return vx_scenarios


def make_sims(location='gabon', calib_pars=None, scenarios=None, stop=end_year):
    """Build all scenario sims for one ss.MultiSim run.

    calib_pars is a list of parsets (or None). Returns (msim, tags) where
    tags[i] = {'scenario': ..., 'parset_idx': ..., 'seed_idx': ...} for
    sim msim.sims[i].
    """
    parsets = calib_pars if calib_pars is not None else [None]

    sims = []
    tags = []
    for name, interventions in scenarios.items():
        for pi, parset in enumerate(parsets):
            for si in range(n_seeds):
                sim = rs.make_sim(location=location, calib_pars=parset, debug=debug,
                                  interventions=interventions, stop=stop,
                                  seed=pi * n_seeds + si, verbose=-1)
                sim.label = name
                sims.append(sim)
                tags.append({'scenario': name, 'parset_idx': pi, 'seed_idx': si})

    msim = ss.MultiSim(sims=sims)
    return msim, tags


def run_sims(location='gabon', calib_pars=None, scenarios=None, verbose=0.2):
    """Run scenarios in a single ss.MultiSim and aggregate per-scenario.

    Uses ss.Result.annualize() to collapse each sim's sub-annual output to
    yearly resolution, then aggregates across sims with pandas: for each
    scenario × metric we return a DataFrame indexed by year with columns
    ['median', 'low', 'high'] (min/max across the parset ensemble).
    """
    import pandas as pd

    msim, tags = make_sims(location=location, calib_pars=calib_pars, scenarios=scenarios)
    n_workers = None if not debug else 1
    msim.run(n_cpus=n_workers)

    for sim in msim.sims:
        sim.shrink()

    metrics = ['asr_cancer_incidence', 'new_cancers', 'new_cancer_deaths']
    msim_dict = sc.objdict()
    scen_labels = list(scenarios.keys())

    # Group sim indices by scenario tag
    by_scen = {name: [] for name in scen_labels}
    for i, tag in enumerate(tags):
        by_scen[tag['scenario']].append(i)

    for scen_label in scen_labels:
        indices = by_scen[scen_label]
        entry = sc.objdict()
        for m in metrics:
            frames = []
            for i in indices:
                ann = msim.sims[i].results['all_hpv'][m].annualize()
                df = ann.to_df()
                df['year'] = pd.to_datetime(df['timevec']).dt.year
                df['sim'] = i
                frames.append(df[['year', 'value', 'sim']])
            long = pd.concat(frames, ignore_index=True)
            agg = long.groupby('year')['value'].agg(median='median', low='min', high='max')
            entry[m] = agg
        msim_dict[scen_label] = entry

    return msim_dict


# %% Run as a script
if __name__ == '__main__':

    os.makedirs('results', exist_ok=True)
    os.makedirs('figures', exist_ok=True)

    T = sc.timer()
    do_run = True
    do_process = True
    location = 'gabon'

    scenarios = dict()
    scenarios['Baseline'] = []

    # Add combined scenarios
    screen_treat_scenarios = make_screen_treat_scenarios()
    vx_scenarios = make_vx_scenarios()
    for st_label, st_intvs in screen_treat_scenarios.items():
        for vx_label, vx_intvs in vx_scenarios.items():
            combined_label = f'{st_label} + {vx_label}'
            combined_intvs = st_intvs + vx_intvs
            scenarios[combined_label] = combined_intvs
    print(f'Set up {len(scenarios)} scenarios to run.')

    if do_run:
        print(f'Running scenarios for location: {location}')

        parsets = rs.load_top_parsets(n_parsets)
        msim_dict = run_sims(location=location, calib_pars=parsets, scenarios=scenarios, verbose=-1)

        if do_process:
            print('Post-processing results...')
            sc.saveobj(f'results/scens_{location}.obj', msim_dict)

    T.toc('Done')