# Table 1 — parameter uncertainty propagated in `run_scenarios.py`

The top 10 calibration parsets (ranked by mismatch, taken from `results/gabon_pars_top50.obj`) drive the parameter uncertainty band in the scenario ensemble. Prior bounds are the calibration search intervals defined in `run_sims.run_calib`. **best** is rank 1 by mismatch (mismatch = 0.875).

| Parameter | prior_low | prior_high | best | top10_min | top10_median | top10_max |
|-----------|-----------|------------|------|-----------|--------------|-----------|
| `beta` | 0.100 | 0.340 | 0.120 | 0.120 | 0.120 | 0.140 |
| `m_cross_layer` | 0.100 | 0.700 | 0.150 | 0.150 | 0.200 | 0.250 |
| `f_cross_layer` | 0.050 | 0.500 | 0.050 | 0.050 | 0.050 | 0.050 |
| `m_partners_c_par1` | 0.100 | 0.600 | 0.100 | 0.100 | 0.140 | 0.560 |
| `f_partners_c_par1` | 0.100 | 0.600 | 0.140 | 0.100 | 0.100 | 0.140 |
| `sev_dist_par1` | 0.500 | 1.500 | 1.360 | 1.360 | 1.400 | 1.500 |
| `hi5_cancer_fn_transform_prob` | 5.00e-04 | 2.50e-03 | 1.10e-03 | 9.00e-04 | 1.30e-03 | 1.50e-03 |
| `hi5_cin_fn_k` | 0.100 | 0.250 | 0.220 | 0.190 | 0.200 | 0.220 |
| `hi5_dur_cin_par1` | 3.500 | 5.500 | 5.000 | 5.000 | 5.000 | 5.500 |
| `hi5_dur_cin_par2` | 16.000 | 24.000 | 24.000 | 24.000 | 24.000 | 24.000 |
| `ohr_cancer_fn_transform_prob` | 5.00e-04 | 2.50e-03 | 1.30e-03 | 5.00e-04 | 5.00e-04 | 1.30e-03 |
| `ohr_cin_fn_k` | 0.100 | 0.250 | 0.130 | 0.100 | 0.110 | 0.130 |
| `ohr_dur_cin_par1` | 3.500 | 5.500 | 4.000 | 4.000 | 5.000 | 5.000 |
| `ohr_dur_cin_par2` | 16.000 | 24.000 | 16.500 | 16.000 | 17.000 | 17.500 |
