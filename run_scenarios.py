'''
Run HPVsim scenarios
Note: requires an HPC to run with debug=False; with debug=True, should take 5-15 min
to run.
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
import hpvsim as hpv

# Imports from this repository
import run_sims as rs


# What to run
debug = 0
n_parsets = [10, 1][debug]  # Top-N calibration parsets, for propagating parameter uncertainty
n_seeds = [1, 1][debug]     # Stochastic seeds per parset per scenario
end_year = 2100  # Simulation horizon; must match the vaccination scale-up schedule below


# %% Functions
def make_screen_treat(screen_coverage=0.15, triage_coverage=0.9, treat_coverage=0.75, start_year=2020):
    """ Make screening & treatment intervention """

    age_range = [30, 50]  # WHO-recommended screening age range
    len_age_range = (age_range[1]-age_range[0])/2
    # Convert age-range coverage to an annual screening probability, assuming a uniform
    # annual hazard of screening over half the age range (i.e. each person is expected
    # to be screened at least once every len_age_range years)
    model_annual_screen_prob = 1 - (1 - screen_coverage)**(1/len_age_range)

    rescreen_interval_years = 5
    screen_eligible = lambda sim: np.isnan(sim.people.date_screened) | \
                                  (sim.t > (sim.people.date_screened + rescreen_interval_years / sim['dt']))
    screening = hpv.routine_screening(
        prob=model_annual_screen_prob,
        eligibility=screen_eligible,
        start_year=start_year,
        product='hpv',
        age_range=age_range,
        label='screening'
    )

    # Assign treatment
    screen_positive = lambda sim: sim.get_intervention('screening').outcomes['positive']
    assign_treatment = hpv.routine_triage(
        start_year=start_year,
        prob=triage_coverage,
        annual_prob=False,
        product='tx_assigner',
        eligibility=screen_positive,
        label='tx assigner'
    )

    ablation_eligible = lambda sim: sim.get_intervention('tx assigner').outcomes['ablation']
    ablation = hpv.treat_num(
        prob=treat_coverage,
        product='ablation',
        eligibility=ablation_eligible,
        label='ablation'
    )

    excision_eligible = lambda sim: np.union1d(sim.get_intervention('tx assigner').outcomes['excision'],
                                                sim.get_intervention('ablation').outcomes['unsuccessful'])
    excision = hpv.treat_num(
        prob=treat_coverage,
        product='excision',
        eligibility=excision_eligible,
        label='excision'
    )

    radiation_eligible = lambda sim: sim.get_intervention('tx assigner').outcomes['radiation']
    radiation = hpv.treat_num(
        prob=treat_coverage/4,  # assume an additional 4x dropoff in coverage for cancer treatment (radiation) vs pre-cancer treatment
        product=hpv.radiation(),
        eligibility=radiation_eligible,
        label='radiation'
    )

    screen_treat_intvs = [screening, assign_treatment, ablation, excision, radiation]

    return screen_treat_intvs


def make_screen_treat_scenarios():
    """ Make screening & treatment scenarios, looping over screening coverage """

    screen_treat_scenarios = dict()

    screen_coverages = [0.1, 0.4, 0.9]  # Low/medium/high screening coverage scenarios
    for scov in screen_coverages:
        label = f'Screen {int(scov*100)}%'
        screen_treat_scenarios[label] = make_screen_treat(screen_coverage=scov)

    return screen_treat_scenarios


def make_vx_scenarios(product='bivalent', start_year=2025, end_year=end_year):
    """ Make vaccination scenarios: no vaccination vs a scale-up to 90% routine coverage """

    routine_age = (9, 10)  # Routine HPV vaccination age
    eligibility = lambda sim: (sim.people.doses == 0)

    vx_scenarios = dict()

    vx_scenarios['No vaccination'] = []

    # Baseline vaccination scenarios
    vx_years = np.arange(start_year, end_year + 1)
    scaleup = [0.3, 0.6, 0.9]  # 3-year scale-up to the target coverage below

    # Maintain 90%
    final_cov = 0.9
    vx_cov = np.concatenate([scaleup+[final_cov]*(len(vx_years)-len(scaleup))])

    routine_vx = hpv.campaign_vx(
        prob=vx_cov,
        years=vx_years,
        product=product,
        age_range=routine_age,
        eligibility=eligibility,
        interpolate=False,
        annual_prob=False,
        label='Routine vx'
    )
    vx_scenarios['90% vax coverage'] = [routine_vx]

    return vx_scenarios


def make_sims(location='gabon', calib_pars=None, scenarios=None, end=end_year):
    """ Set up scenarios. calib_pars is a list of parsets (or None). """

    parsets = calib_pars if calib_pars is not None else [None]

    all_msims = sc.autolist()
    for name, interventions in scenarios.items():
        sims = sc.autolist()
        for pi, parset in enumerate(parsets):
            for si in range(n_seeds):
                sim = rs.make_sim(location=location, calib_pars=parset, debug=debug, interventions=interventions, end=end, seed=pi * n_seeds + si, verbose=-1)
                sim.label = name
                sims += sim
        all_msims += hpv.MultiSim(sims)

    msim = hpv.MultiSim.merge(all_msims, base=False)

    return msim


def run_sims(location='gabon', calib_pars=None, scenarios=None, verbose=0.2):
    """ Make and run the scenario simulations, in parallel unless debug is set """
    msim = make_sims(location=location, calib_pars=calib_pars, scenarios=scenarios)
    parallel = not debug
    msim.run(verbose=verbose, parallel=parallel)
    return msim


# %% Run as a script
if __name__ == '__main__':

    import os
    os.makedirs('results', exist_ok=True)
    os.makedirs('figures', exist_ok=True)

    T = sc.timer()
    do_run = True
    do_save = False
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

    # Run scenarios (usually on VMs, runs n_seeds in parallel over M scenarios)
    if do_run:
        print(f'Running scenarios for location: {location}')

        parsets = sc.loadobj(f'results/{location}_pars_top50.obj')[:n_parsets]
        msim = run_sims(location=location, calib_pars=parsets, scenarios=scenarios, verbose=-1)

        if do_save: msim.save(f'results/scens_{location}.msim')

        if do_process:
            print('Post-processing results...')

            metrics = ['year', 'asr_cancer_incidence', 'cancers', 'cancer_deaths']

            # Process results
            scen_labels = list(scenarios.keys())
            mlist = msim.split(chunks=len(scen_labels))

            msim_dict = sc.objdict()
            for si, scen_label in enumerate(scen_labels):
                reduced_sim = mlist[si].reduce(output=True)
                mres = sc.objdict({metric: reduced_sim.results[metric] for metric in metrics})
                msim_dict[scen_label] = mres

            sc.saveobj(f'results/scens_{location}.obj', msim_dict)

    T.toc('Done')
