# hpvsim_gabon v2.2.6 → v3.1.0 port design

**Date:** 2026-08-25
**Author:** Robyn Stuart (with Claude assistance)
**Status:** Design — awaiting review

## 1. Context

`hpvsim_gabon` targets `hpvsim==2.2.6`. The collaborator will eventually run
these analyses on v3.x, and the two existing v3 references in adjacent repos
give us a well-trodden path:

- **`hpvsim_pxv_younger`** — v2 → v3 migration (spec at `docs/superpowers/specs/2026-08-14-v3-migration-design.md`, approved 2026-08-14). Similar migration surface to gabon; source of the "drop `beta`" finding and the raw_results/ + results/ two-tier layout.
- **`hpvsim_kazakhstan`** — fresh-start v3 build; the reference for 3.1.0 patterns (nested calib_pars with `network`/`cross_immunity`, HPVTotal ASR, `hpv.plot_calibration`, `calib.shrink()`, `raw_results/` + `results/`).

Gabon's v2 code has been fully engineering-uplifted and its Part-1 / Part-2
uncertainty work is committed on `main`. Port scope is now clear: bring gabon
to v3.1.0 on a `v3-port` branch, including a fresh v3 recalibration and
regenerated figures.

## 2. Scope

**In scope:**
- Port every code path exercised today to v3.1.0.
- Recalibrate under v3 with a v3-native prior structure.
- Regenerate the calibration and scenario figures, plus `TABLE1.md`.
- Update all documentation (root `README.md`, `data/README.md`, `results/README.md`) to describe the v3 pipeline.

**Out of scope:**
- Any scientific change (targets stay the same, scenarios stay the same, DHS-derived debut and layer_probs stay the same).
- Adopting sti_notification-style LHS ensemble methodology (not part of gabon's current workflow).
- Any manuscript-facing artifacts beyond `TABLE1.md`.

**Explicit deletions during the port:**
- `utils.shrink_calib` — v3 has `calib.shrink(n_results=N)` built in.
- `conftest.py` — keep. Nothing about v3 changes the sys.path discovery pattern.
- The `sev_dist.par1` prior — renamed to `cross_immunity.rel_sev.loc` (see §4.1).
- The `beta` prior — fix at v2 best (0.120), do not calibrate on v3 (see §4.1).

## 3. Migration surface

Total footprint: ~600 LOC across 4 python files. Concrete API changes required
(from the recon of `hpvsim_kazakhstan` and `hpvsim_pxv_younger` and the v3 hpvsim
source):

### 3.1 `run_sims.py`

**Sim construction:**
- `hpv.Sim(pars=objdict, interventions=..., analyzers=...)` → `hpv.Sim(location=..., n_agents=..., dt=..., start=..., stop=..., genotypes=..., ms_agent_ratio=..., pars=merged_dict, interventions=..., analyzers=..., rand_seed=seed)`.
- `end=2020` → `stop=2020`.
- Network pars (`debut`, `layer_probs`, `m_partners`, `f_partners`) stay in the `pars=` dict; merge via `sc.mergedicts(network_pars_fn(), user_pars)`.

**Debut distributions:** current v2 form is `dict(f=dict(dist='lognormal', par1=..., par2=...))`. v3 wants Starsim distributions — but check `hpvsim_kazakhstan/model.py:13-45` for the exact idiom; the safer path is to keep the flat dict shape v3 still accepts and only touch it if the port fails.

**Calibration:**
- `hpv.Calibration(sim, calib_pars=..., genotype_pars=..., datafiles=[...], total_trials=..., n_workers=..., storage=...)` → `hpv.Calibration(sim, calib_pars=merged_dict, data=[...], total_trials=..., n_workers=..., reseed=False)`.
- `datafiles=` → `data=`.
- Single `calib_pars` dict now (no separate `genotype_pars`); genotype params nest inside the same dict keyed by genotype label. See §4.1 for the exact structure.
- `reseed=False` is the 3.1.0 hpv override default per `feedback_calibration_reseed.md`; do not need to pass explicitly, but pass it anyway for clarity.
- Post-run: `sc.saveobj('raw_results/gabon_calib.obj', calib)`; `shrunk = calib.shrink(n_results=50)`; `sc.saveobj('results/gabon_calib.obj', shrunk)`; `sc.saveobj('results/gabon_pars.obj', calib.best_pars)`.
- `hpv.plot_calibration(calib)` returns the figure — replaces the custom `calib.plot(...)` + font munging block.

**MultiSim / run_parsets:**
- `hpv.MultiSim(sims)` → `ss.MultiSim(sims=sims)`.
- `msim.reduce()` gone. `run_parsets` in v2 used it to build an aggregated Sim; in v3 iterate `msim.sims` post-run and shrink each.

**age_pyramid analyzer:** `hpv.age_pyramid(timepoints=..., datafile=..., edges=...)` unchanged in v3 as a call signature; verify at port time.

### 3.2 `run_scenarios.py`

**Screening probability conversion — semantic change from v2.**
Gabon's v2 `make_screen_treat` converts the input `screen_coverage` to a per-year
prob using `1 - (1 - C)^(1 / (age_range/2))` (i.e. `N = 10` for the 30–50
window), interpreting the input as coverage-per-rescreen-cycle. Adopt
pxv_younger's [`_annual_from_lifetime`](../../../../hpvsim_pxv_younger/run_scenarios.py#L224-L232)
instead: use the **full** age range (`N = 20`), interpreting the input as
**lifetime coverage**. Formula:

```python
SCREEN_AGE_LO, SCREEN_AGE_HI = 30, 50
SCREEN_AGE_YEARS = SCREEN_AGE_HI - SCREEN_AGE_LO   # 20, not 10

def _annual_from_lifetime(lifetime_cov, n_years=SCREEN_AGE_YEARS):
    c = float(min(max(lifetime_cov, 0.0), 1.0))
    if c <= 0: return 0.0
    if c >= 1.0: return 1.0
    return 1.0 - (1.0 - c) ** (1.0 / n_years)
```

`make_screen_treat(screen_coverage=X)` now interprets `X` as lifetime coverage
over ages 30–50, matching pxv_younger. Note this shifts the screening-coverage
axis of the scenario figure — the 10%/40%/90% labels now mean lifetime
probabilities. Document in the run_scenarios docstring.

**Scenario builder:**
- `hpv.routine_screening(prob, eligibility, start_year, product='hpv', age_range=[30, 50], label)` — call signature unchanged; only the per-year `prob` value changes per the conversion above.
- `hpv.routine_triage(prob, product='tx_assigner', eligibility=<lambda>, start_year, annual_prob=False)` — **`annual_prob=False` is required** (gabon v2 already sets it; make sure it survives the port). pxv_younger uses this pattern verbatim ([run_scenarios.py:363-369](../../../../hpvsim_pxv_younger/run_scenarios.py#L363-L369)).
- `hpv.treat_num(prob, product='ablation'|'excision'|hpv.radiation(), eligibility=<lambda>, label)` — same signature. Verify `hpv.radiation()` still exists as a bare class in v3.
- `hpv.campaign_vx(prob, years, product='bivalent', age_range, eligibility, interpolate=False, annual_prob=False, label)` — **still exists in v3**, name unchanged. Product string `'bivalent'` should also survive; fallback is `product=hpv.vx(name='bivalent', sterilizing_p=0.95)`.

**Eligibility callbacks:** rewrite `sim.get_intervention('name')` → `sim.interventions['name']`. Registration order still matters — screening must precede triage which must precede treatment.

**MultiSim + reduction:**
- `hpv.MultiSim(sims)` → `ss.MultiSim(sims=sims)`.
- `MultiSim.merge(list, base=False)` and `msim.split(chunks=N)` both **gone** in v3.
- Replace with: single `ss.MultiSim(sims=all_sims)`, then post-run iterate `msim.sims` grouped by scenario tag, extract each sim's `results['all_hpv']['<metric>']`, and compute `median/min/max` across the group with numpy. Wrap in a small helper that returns `sc.objdict({metric: Result_like_shim(values, low, high)})` so `plot_fig1_residual.py` needs no changes to its consumption pattern.
- `.reduce(output=True)` gone — replaced by the aggregation helper above.

### 3.3 `plot_fig1_residual.py`

Current code reads `msim_dict[scen]['metric'].values/.low/.high`. Keep this
external interface. The Result-like shim from §3.2 delivers exactly this
shape, so `plot_fig1(filestem='')` needs zero changes beyond the aggregation
happening upstream.

Font path (`utils.set_font`) unchanged — the `assets/LibertinusSans-Regular.otf`
committed in Part 2 works on v3.

### 3.4 `utils.py`

- `set_font` — no change.
- `logn_percentiles_to_pars` + `get_debut` — pure math + scipy, no change.
- `plot_single` — no change (matplotlib only).
- **Delete** `shrink_calib` — replaced by `calib.shrink()`.

### 3.5 `make_table1.py`

Parameter names change under v3's calib_pars structure. Update `PRIORS` and
`PATH` mappings:
- Drop `beta` (fixed at 0.120, not calibrated).
- Rename `sev_dist_par1` → `cross_immunity_rel_sev_loc`.
- Rename `m_partners_c_par1` → `network_m_partners_casual`.
- Rename `f_partners_c_par1` → `network_f_partners_casual`.
- Rename `hi5_dur_cin_par1` → `hi5_dur_cin_mean`; `par2` → `std`.
- Same for `ohr`.
- `m_cross_layer` / `f_cross_layer` unchanged.

### 3.6 `tests/`

Rewrite smoke tests against v3 API:
- `test_baseline.py`: `sim = rs.make_sim(debug=1); sim.run(); assert sim.results['all_hpv'].cum_infections.values.sum() > 0`.
- `test_scenarios.py`: `sim = rs.make_sim(interventions=st_intvs + vx_intvs, stop=2030, debug=1); sim.run(); assert 'asr_cancer_incidence' in sim.results['all_hpv']`.
- Remove `conftest.py` if v3 test discovery no longer needs it; add back if bare `pytest tests/` still fails.

### 3.7 `requirements.txt`

- `hpvsim==2.2.6` → `hpvsim==3.1.0`.
- Add `starsim>=1.5` if not pulled transitively.

## 4. Recalibration strategy

### 4.1 Prior structure

v3-native `calib_pars` for gabon (single dict, `[best, low, high]` per param, no step):

```python
def make_calib_pars():
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
```

Search dim = **13** (was 14 on v2; `beta` dropped).

**Beta**: fixed at 0.120 (v2 best) in `make_sim` pars, not calibrated. If v3
prevalence is unrealistically off with beta=0.120, revisit as an open question
(see §9).

**Priors that hit v2 bounds** (per current TABLE1.md): `f_cross_layer` at 0.050,
`hi5_dur_cin_par2` (→ `std`) at 24.0, `ohr_cancer_fn_transform_prob` at 5e-4.
These bounds are widened in the spec above to give the v3 fit room:
- `f_cross_layer` lower bound stays at 0.05 but upper is now 0.5 (was 0.5 on v2 — already generous).
- `hi5.dur_cin.std` upper stays at 24 (widening further risks unrealistic durations).
- `ohr.cancer_fn.transform_prob` lower stays at 5e-4 (widening below invites near-zero cancer output).

If the v3 fit still pushes bounds, revisit these in a follow-up commit.

### 4.2 Trial budget

- **Pilot:** `n_trials=1000, n_workers=100` on **zebra** (IDM Azure). Matches pxv_younger.
- **Storage:** v3 default is JournalStorage — no SQLite lock issues at 100 workers.
- **N_KEEP:** 50 (matches gabon's Part-2 top-50 convention).
- **Assess quality via `hpv.plot_calibration(calib)`.** If fit obviously bad (targets outside model band), scale up trials or adjust priors before the freeze commit.

### 4.3 Reproducibility

- Full `calib` object goes to `raw_results/gabon_calib.obj` (gitignored, ~MB scale).
- Shrunk (`n_results=50`) goes to `results/gabon_calib.obj` (~KB, committed).
- Best pars goes to `results/gabon_pars.obj` (dict, ~KB, committed).
- Top-N parsets file (`results/gabon_pars_top50.obj`) is built from the shrunk calib via `[calib.trial_pars_to_sim_pars(which_pars=i) for i in range(50)]` in the same commit.

## 5. Scenario runs

- Same 7 scenarios as v2 (Baseline + 3 screen × 2 vax).
- Same `n_parsets=10, n_seeds=1` default from Part 2.
- One `ss.MultiSim` of `10 × 1 × 7 = 70` sims. Tag each sim with `(scenario, parset_idx, seed_idx)` before assembling.
- After `msim.run(n_cpus=all)`, iterate sims, group by scenario tag, and compute median / min / max time series across the group.
- Write `results/scens_gabon_top10.obj` in the same shape as Part 2 (`sc.objdict({scen: sc.objdict({metric: Result_shim})})`).

## 6. Verification

**Scoped acceptance criterion: visual equivalence.**

| Figure | v2 baseline | v3 verdict targets |
|---|---|---|
| `gabon_calib.png` | v2 fit at mismatch = 0.875 | Data points inside the model boxplot IQR at every age; ASR inside the boxplot IQR. |
| `gabon_vax_screening_top10.png` | Part 2 top-10 figure (committed) | Rank order of scenarios preserved. **Expect a systematic shift**: the pxv_younger-style lifetime-coverage conversion (§3.2) roughly halves the per-year screening prob compared to v2, so the 10/40/90% scenarios will show weaker screening effects than the v2 figure. That is intentional; still expect Screen 90% + 90% vax to eliminate before 2100 in most parsets. |
| `TABLE1.md` | current v2 table | Similar magnitudes; parameter names updated to v3 idioms; no priors stuck at bounds beyond those noted in §4.1. |

If any figure fails the criterion, treat as a blocker and iterate on priors /
recalibrate before the freeze commit.

**Numeric exactness is not expected** — v3 changes multiscale + network
annualization + genotype prior structure; small shifts are the point.

## 7. Directory layout post-migration

```
data/                          # unchanged
results/                       # committed
├── gabon_calib.obj            # shrunk calib (top-50), ~KB
├── gabon_pars.obj             # best_pars dict
├── gabon_pars_top50.obj       # top-50 parsets list (Part 2 consumer)
├── gabon_pars_all.obj         # top-100 parsets list (kept for compat)
├── table1_parameter_uncertainty.csv
├── scens_gabon_top10.obj      # (new, gitignored — regenerate)
└── README.md
raw_results/                   # gitignored, VM-side only
└── gabon_calib.obj            # full calib
figures/                       # tracked
├── gabon_calib.png            # regenerated via hpv.plot_calibration
└── gabon_vax_screening_top10.png
assets/                        # unchanged
├── LibertinusSans-Regular.otf
docs/superpowers/specs/
└── 2026-08-25-gabon-v3-port-design.md   # this file
run_sims.py                    # v3-ported
run_scenarios.py               # v3-ported
plot_fig1_residual.py          # unchanged interface, new upstream shim
utils.py                       # shrink_calib deleted
make_table1.py                 # parameter names updated
requirements.txt               # hpvsim==3.1.0
tests/                         # rewritten smoke tests
TABLE1.md                      # regenerated
README.md, data/README.md, results/README.md   # updated for v3
```

## 8. Commit sequence

Atomic port: **3 commits** on the `v3-port` branch, each self-contained and
testable.

1. **`Port to hpvsim v3.1.0`**
   - All code changes (§3): `run_sims.py`, `run_scenarios.py`, `plot_fig1_residual.py` (upstream shim), `utils.py`, `make_table1.py`, `requirements.txt`, `tests/`, `.gitignore` (add `raw_results/`).
   - README updates.
   - Keep `conftest.py` (v3 doesn't change the sys.path discovery issue).
   - No calibration artifacts, no regenerated figures — those are the next two commits.
   - Smoke test: `pytest tests/` passes.

2. **`Add v3 recalibration artifacts`**
   - `raw_results/gabon_calib.obj` (gitignored, present locally only).
   - `results/gabon_calib.obj` (shrunk, committed).
   - `results/gabon_pars.obj`, `results/gabon_pars_top50.obj`, `results/gabon_pars_all.obj` (committed).
   - Update `results/README.md` with v3-run details (trial count, best mismatch, date).

3. **`Regenerate figures and TABLE1 under v3`**
   - `figures/gabon_calib.png` — from `hpv.plot_calibration(calib)`.
   - `figures/gabon_vax_screening_top10.png` — from the v3 `plot_fig1(filestem='_top10')`.
   - `TABLE1.md` and `results/table1_parameter_uncertainty.csv` — regenerated by the updated `make_table1.py`.

## 9. Branch strategy

Single feature branch `v3-port` off `main`. When merged, v2 code is gone —
no dual paths, no `_v2.py` files, no `if HPVSIM_VERSION < 3:` shims. Merge
happens when the three commits are complete and the verification in §6 has
all three verdicts of MATCH or SHIFT-BUT-STORY-INTACT.

The current v2 code lives on `main` until merge, so anyone needing the v2
reference can `git checkout <sha-before-merge>`.

## 10. Open questions

Non-blocking; noted so future-me isn't surprised:

- **Debut distribution format**: whether v3 needs the v2 dict form rewritten as `ss.lognorm_ex(...)` or accepts the flat dict as-is. Empirical check at port time.
- **`m_cross_layer` / `f_cross_layer` identification**: under v3's dt-correct network annualization, these may become uninformative. Leave in the prior for the pilot; prune from calib_pars if flat across the top-10.
- **`hpv.radiation()` in v3**: verify the class still exists at the top-level namespace. If renamed, update `run_scenarios.make_screen_treat`.
- **Beta at 0.120 vs hpvsim v3 default**: if v3 prevalence is unrealistically low or high with beta=0.120, revisit — could either widen back into priors or fix at v3 default.
- **Split `run_calib` out of `run_sims.py` into `run_calibration.py`?** Both kazakhstan and pxv_younger keep calibration in a separate file. Gabon's simplicity doesn't force this, but adopting the convention would improve consistency across the three repos. Decide at port time.
- **Whether Part 2's `_top10` filestem survives the port**: worth checking if kazakhstan has a naming convention we should adopt instead. Reserve judgment until the port.