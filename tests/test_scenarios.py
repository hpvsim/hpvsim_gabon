"""Smoke test: screening + treatment + vaccination scenarios build and run."""
import run_sims as rs
import run_scenarios as rsc


def test_screen_treat_and_vx_scenarios_build_and_run():
    st_intvs = rsc.make_screen_treat(screen_coverage=0.4)
    vx_intvs = rsc.make_vx_scenarios()['90% vax coverage']
    sim = rs.make_sim(interventions=st_intvs + vx_intvs, stop=2030, debug=1)
    sim.run()
    assert 'asr_cancer_incidence' in sim.results['all_hpv']