"""Smoke test for the Gabon baseline sim.

No calibrated parameter set is committed yet (gabon awaits recalibration under
hpvsim==2.2.6), so this only checks that an uncalibrated debug-mode sim builds
and runs -- it does not validate against a target ASR.
"""
import run_sims as rs


def test_make_sim_debug_runs():
    sim = rs.make_sim(debug=1)
    sim.run()
    assert sim.results['year'][-1] >= sim['end']
