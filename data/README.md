# data/

Calibration targets and demographic inputs for the Gabon HPV model.

| File | Content |
|---|---|
| `gabon_cancer_cases.csv` | Cervical cancer cases by 5-year age bin (year 2020). Columns: `year, name, age, sex, genotype, value`. Consumed by `hpv.Calibration` as an `age_results` target. |
| `gabon_asr_cancer_incidence.csv` | Age-standardized cervical cancer incidence rate for 2020 (per 100k women). Columns: `year, name, genotype, value`. Consumed by `hpv.Calibration` as a scalar target. |
| `gabon_age_pyramid.csv` | Population age structure by sex (2025). Used by the `hpv.age_pyramid` analyzer when the `age_pyramids` stage is enabled in `run_sims.py`. |

## Provenance

Cancer targets are from **GLOBOCAN 2020** (Gabon national estimates); the ASR
value (30.8/100k) matches the GLOBOCAN 2020 estimate for Gabon.

Sexual-behavior parameters (debut age distributions and partnership layer
probabilities used in `run_sims.make_sim`) are derived by fitting to the 2019–21
Gabon Demographic and Health Survey; see
https://www.researchsquare.com/article/rs-3074559/v1.
