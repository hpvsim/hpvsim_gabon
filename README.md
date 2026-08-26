# hpvsim_gabon

An [HPVsim](https://hpvsim.org) model of cervical cancer for Gabon, calibrated
to national cancer case and age-standardized incidence data.

Built on **hpvsim v3.1.0** (starsim-based). See
[`docs/superpowers/specs/2026-08-25-gabon-v3-port-design.md`](docs/superpowers/specs/2026-08-25-gabon-v3-port-design.md)
for the migration record.

## Install

```bash
conda create -n hpvsim_gabon python=3.11 -y
conda activate hpvsim_gabon
pip install -r requirements.txt
```

`hpvsim==3.1.0` is not yet on PyPI; `requirements.txt` installs it from the
`starsimhub/hpvsim@rc3.1.0` git branch. Switch to `hpvsim==3.1.0` once the
release lands.

## What's here

| File / dir | Purpose |
|---|---|
| `run_sims.py` | `make_sim` for Gabon (v3 `hpv.Sim`), `make_calib_pars`, `run_calib`, `plot_calib`, `run_parsets` (via `ss.MultiSim`). |
| `run_scenarios.py` | Screening × vaccination scenarios; iterates over the top-N calibration parsets for parameter uncertainty. Screening prob is derived via `_annual_from_lifetime` (spec §3.2). |
| `plot_fig1_residual.py` | Plots ASR trajectory + cumulative cancers; `plot_fig1(filestem='')` selects `results/scens_gabon{filestem}.obj`. |
| `make_table1.py` | Regenerates `TABLE1.md` + `results/table1_parameter_uncertainty.csv` from the shrunk `results/gabon_calib.obj`. |
| `utils.py` | Font loader, DHS-derived debut helpers, plotting helpers. |
| `data/` | Calibration targets. See [`data/README.md`](data/README.md). |
| `results/` | Committed shrunk calib + best pars + Table 1 CSV. See [`results/README.md`](results/README.md). |
| `raw_results/` | Full (unshrunk) calibration object (gitignored). Populated by `run_calib` on the VM. |
| `figures/` | `gabon_calib.png` (calibration fit), `gabon_vax_screening.png` / `_top10.png` (scenarios). |
| `assets/` | `LibertinusSans-Regular.otf` — plot font. |
| `tests/` + `conftest.py` | Smoke tests for the baseline sim and scenario builders. |
| `docs/superpowers/specs/` + `plans/` | Design spec and implementation plan for the v3 port. |

## How to reproduce

Each script has a `to_run` / `do_run` block near its `__main__` that toggles stages.

```bash
python run_sims.py          # calibration + single-sim helpers
python run_scenarios.py     # scenario runs (uses top-N parsets from results/gabon_calib.obj)
python plot_fig1_residual.py
```

Calibration (`run_calib`, 1000 trials × 100 workers by default) and full
scenario runs need an HPC / VM. Debug mode (`debug=1` at the top of each file)
drops population size and disables parallelism, taking a few minutes on a
laptop.

## Parameter uncertainty

`run_scenarios.py` propagates parameter uncertainty by iterating over the top-N
calibration parsets (`n_parsets=10` by default) with `n_seeds=1` per parset. The
scenario uncertainty band reflects parameter variation across the best
calibration draws, not stochastic replicates.

[`TABLE1.md`](TABLE1.md) lists each calibrated parameter with its prior bounds
and the min/median/max/best across the top-10 parsets — this is what the
ensemble propagates through the scenarios. Regenerate after re-calibrating with:

```bash
python make_table1.py
```

## Calibration

- Priors are defined in `run_sims.make_calib_pars` — 13 parameters in a nested
  dict following the kazakhstan v3 convention (`network.*`,
  `cross_immunity.rel_sev.loc`, per-genotype `cancer_fn`/`cin_fn`/`dur_cin`).
- `beta` is fixed at `0.120` (v2 best) inside `make_sim` and is NOT calibrated
  on v3 (spec §4.1).
- Under v3, `hpv.Calibration` uses JournalStorage by default; 100 workers is
  safe (unlike v2's SQLite backend, which contended above ~30).
- Post-calibration, the full object is saved to `raw_results/gabon_calib.obj`
  and the shrunk (top-50) copy to `results/gabon_calib.obj` via
  `calib.shrink(n_results=50)`.

## Testing

```bash
pytest tests/
```

The tests build a debug-mode sim and scenario intervention list and confirm
they run to completion; they do not validate the calibration fit.