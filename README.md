# hpvsim_gabon

An [HPVsim](https://hpvsim.org) model of cervical cancer for Gabon, calibrated to
national cancer case and age-standardized incidence data. Built on **hpvsim v2.x**
(not yet migrated to v3.x).

## Install

```bash
pip install -r requirements.txt
```

Requires `hpvsim==2.2.6`.

## What's here

| File | Purpose |
|------|---------|
| `run_sims.py` | Defines the Gabon simulation (`make_sim`), plus calibration (`run_calib`) and single/multi-run helpers. |
| `run_scenarios.py` | Combined screening & treatment / vaccination scenarios, run in parallel over seeds. |
| `plot_fig1_residual.py` | Plots residual cervical cancer burden across scenarios. |
| `utils.py` | Fonts, DHS-derived sexual-debut distributions, calibration/plotting helpers. |
| `data/` | Calibration targets (cancer cases, ASR incidence, age pyramid). |
| `tests/` | Smoke tests: baseline sim and scenario builders run in debug mode. |

## Calibration status

No calibrated parameter set is committed yet — a recalibration under hpvsim==2.2.6
is pending. `run_calib`/`plot_calib` in `run_sims.py` produce `results/gabon_pars.obj`
once run.

## How to run

Each script has a `to_run` (or `do_run`/`do_save`) list near its `__main__` block —
edit that list to select which stage to run.

```bash
python run_sims.py          # calibration + single-sim helpers; see to_run in the file
python run_scenarios.py     # scenario runs (debug=True: ~5-15 min; full run needs an HPC/VM)
python plot_fig1_residual.py
```

Calibration (`run_calib` in `run_sims.py`) and full scenario runs (`run_scenarios.py`
with `debug=0`) are long-running and should be run on a VM, not locally.

## Testing

```bash
pytest tests/
```

## Data provenance

Sexual debut and partnership parameters are derived by fitting to the 2019-21 Gabon
DHS; see https://www.researchsquare.com/article/rs-3074559/v1.
