# Table 1 — parameter uncertainty propagated in `run_scenarios.py`

The top 10 calibration parsets (ranked by mismatch, taken from the top-10 rows of `results/gabon_calib.obj`) drive the parameter uncertainty band in the scenario ensemble. Prior bounds are the calibration search intervals defined in `run_sims.make_calib_pars`. **best** is rank 1 by mismatch (mismatch = 2.393).

| Parameter | prior_low | prior_high | best | top10_min | top10_median | top10_max |
|-----------|-----------|------------|------|-----------|--------------|-----------|
| `m_cross_layer` | 0.100 | 0.700 | 0.699 | 0.635 | 0.662 | 0.700 |
| `f_cross_layer` | 0.050 | 0.500 | 0.487 | 0.435 | 0.489 | 0.500 |
| `network_m_partners_casual` | 0.100 | 0.600 | 0.422 | 0.422 | 0.527 | 0.571 |
| `network_f_partners_casual` | 0.100 | 0.600 | 0.589 | 0.426 | 0.576 | 0.593 |
| `cross_immunity_rel_sev_loc` | 0.500 | 1.500 | 1.398 | 1.263 | 1.445 | 1.484 |
| `hi5_cancer_fn_transform_prob` | 5.00e-04 | 2.50e-03 | 2.04e-03 | 1.79e-03 | 2.20e-03 | 2.50e-03 |
| `hi5_cin_fn_k` | 0.100 | 0.250 | 0.116 | 0.116 | 0.199 | 0.225 |
| `hi5_dur_cin_mean` | 3.500 | 5.500 | 4.234 | 3.515 | 4.044 | 4.810 |
| `hi5_dur_cin_std` | 16.000 | 24.000 | 20.830 | 16.691 | 20.545 | 22.416 |
| `ohr_cancer_fn_transform_prob` | 5.00e-04 | 2.50e-03 | 1.98e-03 | 1.20e-03 | 1.38e-03 | 1.98e-03 |
| `ohr_cin_fn_k` | 0.100 | 0.250 | 0.129 | 0.108 | 0.188 | 0.232 |
| `ohr_dur_cin_mean` | 3.500 | 5.500 | 5.276 | 3.500 | 4.994 | 5.276 |
| `ohr_dur_cin_std` | 16.000 | 24.000 | 21.613 | 16.330 | 18.996 | 21.942 |
