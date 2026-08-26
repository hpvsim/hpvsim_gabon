"""Smoke test for the Gabon baseline sim under v3."""
import run_sims as rs


def test_make_sim_debug_runs():
    sim = rs.make_sim(debug=1, stop=2000)
    sim.run()
    r = sim.results['all_hpv']
    assert float(r['cum_infections'].values[-1]) > 0, \
        'sim ran but produced no infections — network or seeding broken'