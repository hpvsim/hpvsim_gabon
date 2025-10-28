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
n_seeds = [10, 1][debug]  # How many seeds to run per cluster


# %% Functions
def make_st(screen_coverage=0.15, treat_coverage=0.7, start_year=2020):
    """ Make screening & treatment intervention """

    age_range = [30, 50]
    len_age_range = (age_range[1]-age_range[0])/2
    model_annual_screen_prob = 1 - (1 - screen_coverage)**(1/len_age_range)

    screen_eligible = lambda sim: np.isnan(sim.people.date_screened) | \
                                  (sim.t > (sim.people.date_screened + 5 / sim['dt']))
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
        prob=1.0,
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

    excision_eligible = lambda sim: list(set(sim.get_intervention('tx assigner').outcomes['excision'].tolist() +
                                             sim.get_intervention('ablation').outcomes['unsuccessful'].tolist()))
    excision = hpv.treat_num(
        prob=treat_coverage,
        product='excision',
        eligibility=excision_eligible,
        label='excision'
    )

    radiation_eligible = lambda sim: sim.get_intervention('tx assigner').outcomes['radiation']
    radiation = hpv.treat_num(
        prob=treat_coverage/4,  # assume an additional dropoff in CaTx coverage
        product=hpv.radiation(),
        eligibility=radiation_eligible,
        label='radiation'
    )

    st_intvs = [screening, assign_treatment, ablation, excision, radiation]

    return st_intvs


def make_vx_scenarios(product='bivalent', start_year=2025):

    routine_age = (9, 10)
    eligibility = lambda sim: (sim.people.doses == 0)

    vx_scenarios = dict()

    # No vaccination
    vx_scenarios['No vaccination'] = []

    # Baseline vaccination scenarios
    vx_years = np.arange(start_year, 2100 + 1)
    scaleup = [0.3, 0.6, 0.9]

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
    vx_scenarios['90% routine coverage'] = [routine_vx]

    return vx_scenarios


def make_sims(location='gabon', calib_pars=None, vx_scenarios=None, end=2100):
    """ Set up scenarios """

    st_intv = make_st()

    all_msims = sc.autolist()
    for name, vx_intv in vx_scenarios.items():
        sims = sc.autolist()
        for seed in range(n_seeds):
            interventions = vx_intv + st_intv
            sim = rs.make_sim(location=location, calib_pars=calib_pars, debug=debug, interventions=interventions, end=end, seed=seed, verbose=-1)
            sim.label = name
            sims += sim
        all_msims += hpv.MultiSim(sims)

    msim = hpv.MultiSim.merge(all_msims, base=False)

    return msim


def run_sims(location='gabon', calib_pars=None, vx_scenarios=None, verbose=0.2):
    """ Run the simulations """
    msim = make_sims(location=location, calib_pars=calib_pars, vx_scenarios=vx_scenarios)
    parallel = ~(debug)
    msim.run(verbose=verbose, parallel=parallel)
    return msim


# %% Run as a script
if __name__ == '__main__':

    T = sc.timer()
    do_run = True
    do_save = False 
    do_process = True
    location = 'gabon'

    # Run scenarios (usually on VMs, runs n_seeds in parallel over M scenarios)
    if do_run:
        calib_pars = sc.loadobj(f'results/{location}_pars.obj')
        vx_scenarios = make_vx_scenarios()
        msim = run_sims(location=location, calib_pars=calib_pars, vx_scenarios=vx_scenarios)

        if do_save: msim.save(f'results/vs_{location}.msim')

        if do_process:

            metrics = ['year', 'asr_cancer_incidence', 'cancers', 'cancer_deaths']

            # Process results
            vx_scenarios = make_vx_scenarios()
            scen_labels = list(vx_scenarios.keys())
            mlist = msim.split(chunks=len(scen_labels))

            msim_dict = sc.objdict()
            for si, scen_label in enumerate(scen_labels):
                reduced_sim = mlist[si].reduce(output=True)
                mres = sc.objdict({metric: reduced_sim.results[metric] for metric in metrics})
                msim_dict[scen_label] = mres

            sc.saveobj(f'results/vx_scens_{location}.obj', msim_dict)

    print('Done.')
