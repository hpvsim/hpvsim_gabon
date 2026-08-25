# hpvsim_gabon

An [HPVsim](https://hpvsim.org) model of cervical cancer for Gabon, calibrated
to national cancer case and age-standardized incidence data.

Built on **hpvsim v2.2.6** (not yet migrated to v3.x).

## Install

```bash
conda create -n hpvsim_gabon python=3.11 -y
conda activate hpvsim_gabon
pip install -r requirements.txt
```

Requires `hpvsim==2.2.6`.

## What's here

| File / dir | Purpose |
|---|---|
| `run_sims.py` | `make_sim` for Gabon, plus `run_calib` and single/multi-run helpers. |
| `run_scenarios.py` | Screening × vaccination scenarios; iterates over the top-N calibration parsets for parameter uncertainty. |
| `plot_fig1_residual.py` | Plots ASR trajectory + cumulative cancers; `plot_fig1(filestem='')` selects the source `scens_gabon{filestem}.obj`. |
| `make_table1.py` | Regenerates `TABLE1.md` + `results/table1_parameter_uncertainty.csv` from the committed top-N parsets. |
| `utils.py` | Font loader, DHS-derived debut helpers, plotting helpers. |
| `data/` | Calibration targets. See [`data/README.md`](data/README.md). |
| `results/` | Committed calibration parsets + Table 1 CSV. See [`results/README.md`](results/README.md). |
| `figures/` | `gabon_calib.png` (calibration fit), `gabon_vax_screening.png` / `_top10.png` (scenarios). |
| `assets/` | `LibertinusSans-Regular.otf` — plot font. |
| `tests/` + `conftest.py` | Smoke tests for the baseline sim and scenario builders. |

## How to reproduce

Each script has a `to_run` / `do_run` block near its `__main__` that toggles stages.

```bash
python run_sims.py          # calibration + single-sim helpers
python run_scenarios.py     # scenario runs (uses top-N parsets from results/)
python plot_fig1_residual.py
```

Calibration (`run_calib`, 2500 trials × 50 workers by default) and full scenario
runs need an HPC / VM. Debug mode (`debug=1` at the top of each file) drops
population size and disables parallelism, taking ~5–15 min on a laptop.

## Parameter uncertainty

`run_scenarios.py` propagates parameter uncertainty by iterating over the top-N
calibration parsets (`n_parsets=10` by default) with `n_seeds=1` per parset. The
scenario uncertainty band reflects parameter variation across the best calibration
draws, not stochastic replicates.

[`TABLE1.md`](TABLE1.md) lists each calibrated parameter with its prior bounds
and the min/median/max/best across the top-10 parsets — this is what the
ensemble propagates through the scenarios. Regenerate after re-calibrating with:

```bash
python make_table1.py
```

## Calibration

The committed calibration was run under `hpvsim==2.2.6` (2500 trials × 50 workers,
best mismatch = 0.875). The full `gabon_calib.obj` (unshrunk) is not committed;
`results/gabon_pars_all.obj` and `results/gabon_pars_top50.obj` are.

**Storage caveat**: the default optuna backend is a local SQLite file. With ≥50
concurrent workers this can hit `sqlite3.OperationalError: database is locked`
and kill the run. If that happens, reduce `n_workers` (~20–30 is safe) or
configure a non-SQLite optuna storage backend.

## Testing

```bash
pytest tests/
```

The tests build a debug-mode sim and scenario intervention list and confirm they
run to completion; they do not validate the calibration fit.
