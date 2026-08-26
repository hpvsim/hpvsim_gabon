"""
Define an HPVsim simulation for Gabon, including calibration
"""

# Standard imports
import os
import numpy as np
import sciris as sc
import starsim as ss
import hpvsim as hpv

# Imports from this repository
import utils as ut

# %% Settings and filepaths

# Debug switch
debug = 0  # Run with smaller population sizes and in serial
do_shrink = True  # Do not keep people when running sims (saves memory)

# Run settings (v3 JournalStorage handles ~100 workers safely; see spec §4.2)
n_trials    = [1000, 2][debug]
n_workers   = [100, 1][debug]

# Save settings
do_save = True
save_plots = True

# Top-N calibration trials to keep in the shrunken (committable) calib object
N_KEEP = 50


def make_calib_pars():
    """v3-native calibration priors for Gabon. Format: [best, low, high] per param.

    Beta is fixed in make_sim (see spec §4.1). Priors follow kazakhstan's
    nesting convention (network sub-dict, cross_immunity.rel_sev.loc, dur_cin
    mean/std). See docs/superpowers/specs/2026-08-25-gabon-v3-port-design.md §4.1.
    """
    pars = dict(
        m_cross_layer=[0.15, 0.1, 0.7],
        f_cross_layer=[0.1, 0.05, 0.5],
        network=dict(
            m_partners_casual=[0.2, 0.1, 0.6],
            f_partners_casual=[0.2, 0.1, 0.6],
        ),
        cross_immunity=dict(rel_sev=dict(loc=[1.0, 0.5, 1.5])),
    )
    for g in ['hi5', 'ohr']:
        pars[g] = dict(
            cancer_fn=dict(transform_prob=[1.5e-3, 0.5e-3, 2.5e-3]),
            cin_fn=dict(k=[0.15, 0.1, 0.25]),
            dur_cin=dict(mean=[4.5, 3.5, 5.5], std=[20, 16, 24]),
        )
    return pars


# %% Simulation creation functions
def make_sim(location='gabon', calib_pars=None, debug=0, interventions=None, analyzers=None, seed=1, stop=2020, verbose=0.1):
    """Define parameters, analyzers, and interventions for the Gabon simulation."""
    pars = sc.objdict(
        beta=0.120,  # fixed at v2 best (dropped from priors on v3, see spec §4.1)
        verbose=verbose,
        rand_seed=seed,
    )

    # Sexual debut (fitted to 2019-21 Gabon DHS)
    pars.debut_f = ss.lognorm_ex(mean=17.29, std=2.54)
    pars.debut_m = ss.lognorm_ex(mean=17.65, std=3.15)

    # Marital layer age-participation (fitted to 2019-21 DHS)
    pars.layer_probs_marital = np.array([
        [0, 5, 10,   15,   20,   25,   30,   35,   40,   45,   50,   55,   60,   65,    70,   75],
        [0, 0,  0,  0.1,  0.1, 0.15, 0.15, 0.15,  0.2,  0.3,  0.4,  0.4,  0.2, 0.07, 0.035, 0.007],
        [0, 0,  0,  0.1,  0.1, 0.15, 0.15,  0.2,  0.2,  0.4,  0.4,  0.4,  0.2,  0.1,  0.05, 0.01 ],
    ])

    # Casual layer age-participation (fitted to 2019-21 DHS)
    pars.layer_probs_casual = np.array([
        [0, 5, 10,  15,  20,  25,  30,  35,  40,  45,  50,  55,   60,   65,   70,   75],
        [0, 0, 0.2, 0.4, 0.5, 0.6, 0.5, 0.4, 0.4, 0.4, 0.3, 0.2, 0.10, 0.02, 0.02, 0.02],
        [0, 0, 0.2, 0.4, 0.4, 0.5, 0.6, 0.5, 0.5, 0.4, 0.4, 0.2, 0.02, 0.02, 0.02, 0.02],
    ])

    # Number of concurrent partners per layer
    pars.m_partners_marital = 0.01
    pars.m_partners_casual = 0.2
    pars.f_partners_marital = 0.01
    pars.f_partners_casual = 0.2

    if analyzers is None:
        analyzers = []

    sim = hpv.Sim(
        location=location,
        n_agents=[10e3, 1e3][debug],
        dt=[0.25, 1.0][debug],
        start=[1960, 1980][debug],
        stop=stop,
        genotypes=[16, 18, 'hi5', 'ohr'],
        ms_agent_ratio=100,
        pars=pars,
        interventions=interventions,
        analyzers=analyzers,
    )

    if calib_pars is not None:
        hpv.route_pars(sim, dict(calib_pars))

    return sim


# %% Simulation running functions
def run_sim(calib_pars=None, analyzers=None, debug=debug, seed=1, verbose=.1, do_shrink=do_shrink, do_save=do_save, stop=2020):
    """Run a single Gabon simulation."""
    sim = make_sim(debug=debug, seed=seed, analyzers=analyzers, calib_pars=calib_pars, stop=stop)
    sim.label = f'Sim-{seed}'
    sim['verbose'] = verbose
    sim.run()
    if do_shrink:
        sim.shrink()
    if do_save:
        sim.save(f'results/gabon.sim')
    return sim


def load_top_parsets(n, calib_path='results/gabon_calib.obj'):
    """Load top-n calibration parsets by mismatch as a list of flat dicts with
    dotted keys (e.g. ``network.m_partners_casual``). Apply to a sim via
    ``hpv.route_pars(sim, parset)``. Matches the pxv_younger/kazakhstan idiom.
    """
    calib = sc.load(calib_path)
    top = calib.df.nsmallest(n, 'mismatch').reset_index(drop=True)
    par_cols = [c for c in top.columns if c not in ('index', 'mismatch', 'rand_seed')]
    return [{c: row[c] for c in par_cols} for _, row in top.iterrows()]


def run_calib(n_trials=None, n_workers=None, do_save=True, do_plot=False, filestem=''):
    """Calibrate the model. Saves full calib to raw_results/, shrunk (top-N) to results/, best_pars to results/."""
    sim = make_sim()
    data = [
        'data/gabon_cancer_cases.csv',
        'data/gabon_asr_cancer_incidence.csv',
    ]

    calib = hpv.Calibration(
        sim,
        calib_pars=make_calib_pars(),
        data=data,
        total_trials=n_trials,
        n_workers=n_workers,
        reseed=False,
    )
    try:
        calib.calibrate()
    except Exception as e:
        print(f'calibrate() raised: {e}; saving partial results anyway')

    if do_save:
        os.makedirs('raw_results', exist_ok=True)
        os.makedirs('results', exist_ok=True)
        sc.saveobj(f'raw_results/gabon_calib{filestem}.obj', calib)
        shrunk = calib.shrink(n_results=N_KEEP)
        sc.saveobj(f'results/gabon_calib{filestem}.obj', shrunk)
        sc.saveobj(f'results/gabon_pars{filestem}.obj', calib.best_pars)

    if do_plot:
        os.makedirs('figures', exist_ok=True)
        fig = hpv.plot_calibration(calib)
        fig.savefig(f'figures/gabon_calib{filestem}.png')

    if getattr(calib, 'best_pars', None) is not None:
        print(f'Best pars: {calib.best_pars}')
    return sim, calib


def plot_calib(filestem=''):
    """Load the shrunk calib from results/ and render the calibration figure."""
    ut.set_font()
    calib = sc.load(f'results/gabon_calib{filestem}.obj')
    fig = hpv.plot_calibration(calib)
    os.makedirs('figures', exist_ok=True)
    fig.savefig(f'figures/gabon_calib{filestem}.png', dpi=120)
    return fig


def run_parsets(debug=False, verbose=.1, analyzers=None, save_results=True, stop=2040, **kwargs):
    """Run the top-N calibrated parsets in parallel via ss.MultiSim."""
    parsets = load_top_parsets(N_KEEP)

    sims = []
    for i, parset in enumerate(parsets):
        sim = make_sim(debug=debug, seed=i, analyzers=analyzers, calib_pars=parset, stop=stop, verbose=verbose)
        sims.append(sim)

    msim = ss.MultiSim(sims=sims)
    msim.run(n_cpus=None)

    for sim in msim.sims:
        sim.shrink()

    if save_results:
        sc.saveobj('results/gabon_msim.obj', [s.results for s in msim.sims])

    return msim


# %% Run as a script
if __name__ == '__main__':

    os.makedirs('results', exist_ok=True)
    os.makedirs('figures', exist_ok=True)

    # List of what to run
    to_run = [
        # 'run_sim',
        # 'age_pyramids',
        # 'run_calib',
        'plot_calib'
        # 'run_parsets'
    ]

    T = sc.timer()

    if 'run_sim' in to_run:
        calib_pars = None
        sim = run_sim(calib_pars=calib_pars, do_save=False, do_shrink=False)
        sim.plot()

    if 'age_pyramids' in to_run:
        calib_pars = sc.loadobj('results/gabon_pars.obj')
        ap = hpv.age_pyramid(
            timepoints=['2025', '2050', '2075', '2100'],
            datafile='data/gabon_age_pyramid.csv',
            edges=np.linspace(0, 100, 21),
        )
        sim = run_sim(stop=2100, calib_pars=calib_pars, analyzers=[ap], do_save=True, do_shrink=True)

    if 'run_calib' in to_run:
        sim, calib = run_calib(n_trials=n_trials, n_workers=n_workers, filestem='', do_save=True)

    if 'plot_calib' in to_run:
        plot_calib(filestem='')

    if 'run_parsets' in to_run:
        msim = run_parsets()

    T.toc('Done')
