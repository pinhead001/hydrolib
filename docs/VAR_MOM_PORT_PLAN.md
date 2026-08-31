# Plan — closing the open numerical defects (TODO.md P3)

Status: **Phase 0 landed**; Phases 1-7 not started.
Written 2026-08-31 against `claude/todo-open-items-var-mom-709p3m`.

TODO.md P3 lists three open numerical items — weighted skew, the ADJE at-site skew MSE,
and confidence-interval shape — and states that all three bottom out in `var_mom`. That is
still broadly right, but two things measured while writing this plan change the order of
the work and shrink it.

---

## 1. What was measured, and what it changes

### 1.1 The censored-path defect is not in the truncated-moment code

TODO.md says:

> `_ema_iteration` reproduces `moms_p3` **exactly** on uncensored rows [...] and diverges
> only where intervals are censored [...] the residual is in the truncated-P3 moment code
> for censored intervals — that is what the port has to replace.

The first half is right; the conclusion is not. `moms_p3` reads its bias-correction factor
from `common /tac002/`, and the blockdata at `emafit.f:3898` sets `data bcf/1997/`. Under
`bcf = 1997` the factors are built from `n_bcf = n_e`, the count of **exactly observed**
peaks (`emafit.f:1407-1410`):

```fortran
if(bcf .eq. 1997) then
  n_bcf = n_e                                    ! exact peaks only
  c2    = n_bcf/(n_bcf-1.d0)
  c3    = (n_bcf**2)/((n_bcf-1.d0)*(n_bcf-2.d0))
```

`_ema_iteration` uses `n_bcf = n`, the total interval count — the `bcf = 2004` branch, which
is commented out one line below the default (`emafit.f:3899`). When nothing is censored
`n_e == n` and the two agree identically, which is exactly why Powder River is exact to
machine precision and why `TestMomentIterationOracle::test_exact_on_uncensored_rows` passes.
No existing test can distinguish the two.

Backing the per-interval non-central moments out of `moms_p3` and comparing them against
`_compute_ema_moments`, with `n_bcf = n_e` in the back-out:

| censored interval (Cains Coulee moments) | E[X] | E[X²] | E[X³] |
|---|---:|---:|---:|
| `[1e-20, 332]` (left-censored PILF) | 1.0e-10 | 2.3e-10 | 3.9e-10 |
| `[800, 1200]` | 1.6e-12 | 3.3e-12 | 5.2e-12 |
| `[1500, 3000]` | 8.5e-13 | 1.7e-12 | 2.6e-12 |

Relative difference, hydrolib against Fortran. **The truncated-P3 moment code is already
exact.** Nothing in it needs replacing. The same back-out with `n_bcf = n` reproduces the
0.15%/0.5% errors that were being attributed to it.

### 1.2 What the one-line correction is worth

Prototyped by monkey-patching `_ema_iteration` to use `sums.n_exact`; no other change.

| case | quantity | current | `n_bcf = n_e` | peakfq 8.1.0 |
|---|---|---:|---:|---:|
| Cains Coulee | at-site skew | −0.82996 | **−0.70791** | −0.70789 |
| Cains Coulee | std_log | 3.15% off | 3.05% off | 0.380882 |
| Big Sandy | at-site skew | +0.00196 | **+0.00675** | +0.00660 |
| Big Sandy | std_log | 0.634% off | **0.218% off** | 0.291043 |
| Big Sandy | mean_log | 0.006% off | 0.015% off | 3.717508 |
| Powder River | every moment | exact | exact | — |

The Cains Coulee at-site skew gap goes from 0.122 to **2e-5**. That is the `skew_at_site`
half of the `xfail(strict=True)` in `test_wymt_vs_golden.py` — it will start passing, and
strict means the build goes red until the parametrisation is split.

Two things get slightly worse and should not be argued with: Big Sandy's weighted skew
(−0.1187 → −0.1157 against −0.1563) and its Q100 (0.88% → 1.6%). Both are downstream of the
at-site skew through `_b17b_skew_mse` → `nG`, and the weighting path is still missing ADJE
and `Wd`. They are Phases 4 and 5 below. The at-site fit is what `bcf` governs, and the
at-site fit is now right.

**This is why Phase 0 comes before the port.** `mseg_all`, `detrat`, `var_mom` and `regmoms`
are all evaluated *at the at-site moments*. Calibrating a port against an at-site skew that
is 0.12 wrong means calibrating against a poisoned input.

### 1.3 The oracle surface can cover the entire dependency tree

`build_fortran/_emafort.pyf` exposes four routines today, and its header explains the choice
by noting that all four are already callable from peakfqr's R code via `.Fortran()`. That
constraint is R's, not f2py's: f2py will wrap any symbol present in the compiled objects,
including plain Fortran `function`s that `.Fortran()` cannot reach.

Verified this session — `mse_ema`, `mn2mvarb`, `regmoms` and `qP3sub` added to the pyf, the
extension rebuilt from the unmodified vendored sources, all four called successfully, and
`mseg_all_sub` still reproducing the golden `as_G_mse_o` bit for bit. **No edit to
`vendor/`.** The scratch change was reverted; Phase 1 does it properly.

That means every phase below can be diffed against the Fortran *at its own level* rather
than through `emafitpr`'s ~1e13-conditioned fixed point.

---

## 2. Phases

Estimates are for someone who has read the Fortran. Each phase ends green, so the sequence
can stop at any phase boundary without leaving the tree broken.

### Phase 0 — bias-correction factor (½ day) — **DONE**

**Change.** In `_ema_iteration`, build `c2`/`c3` from `sums.n_exact` rather than `n`. Keep
the divisor `n` — `moms(2) = (c2*s_e(2) + s_c(2))/n` uses the total. Guard `n_exact < 3`,
which `moms_p3` does not (it would divide by zero); log and fall back to no correction, the
Fortran's `bcf = 0` branch.

**Oracle.** `moms_p3`, already wired. Extend `TestMomentIterationOracle` with a censored
case asserting equality to ~1e-9 instead of today's "0.7% variance, 4.9% skew" bounds.

**Acceptance.**
- `test_censored_rows_are_where_it_diverges` inverts: it becomes an equality test. Rename it.
- Split the Cains Coulee `test_skew_matches` parametrisation — `skew_at_site` passes (2e-5
  against a 0.02 bound), `skew_weighted` stays `xfail(strict=True)` pending Phases 4-5.
- Re-measure and re-record every tolerance in `test_native_vs_golden.py`,
  `test_wymt_vs_golden.py` and `tests/validation/`. Several are pinned near current values;
  a *better* fit can trip an upper bound that was calibrated as a ceiling. Do not widen a
  bound to make this pass without saying why in the same commit.
- `hydrolib/validation/data/big_sandy_03606500.json` tolerances are measured values; regenerate.

**Trap.** `_ExpectedSums.central()` already splits exact from censored correctly and needs no
change. The bug is only in the two coefficients.

**What actually happened.** The change landed as `_bias_correction_factors`, a module-level
helper, so the `n_exact < 3` guard and its provenance have somewhere to live. Exactly two
tests failed, both predicted: the oracle's censored-divergence assertion (now an equality
test at 1e-8, measured 1.7e-10) and the Cains Coulee `skew_at_site` XPASS. Nothing else in
the suite moved.

The tolerance re-measurement was the larger half of the work, and it found the thing worth
recording. Big Sandy's `std_log` error fell 4x and its at-site skew is now 1.4e-4 off, but
**three quantile bounds tightened their margin instead of loosening it** — Big Sandy's Q100
went 0.88% → 1.59% and Cains Coulee's 2.7% → 4.7%. Both are downstream of the weighted skew,
which still lacks ADJE and `Wd`: a *correct* at-site skew feeds a still-wrong `nG`, so the
fit gets further from peakfq's quantiles before Phases 4-5 bring it back. No bound was
widened. Three now sit at 1.06-1.19x headroom and were left there deliberately.

Measured, for comparison against later phases:

| | Big Sandy | Powder River | Cains Coulee |
|---|---:|---:|---:|
| mean_log | 0.015% | 0.000% | 0.0064 abs |
| std_log | 0.218% | 0.000% | 3.05% |
| at-site skew | 1.4e-4 | 0.0 | 2.0e-5 |
| weighted skew | 0.0407 | 0.0 | 0.1243 |
| worst quantile | 2.82% @ 0.002 | 0.10% @ 0.002 | 21.1% @ 0.995 |
| Q100 | 1.59% | 0.04% | 4.69% |

### Phase 1 — extend the oracle surface (½ day)

**Change.** Add to `_emafort.pyf`: `mse_ema`, `mn2mvarb`, `mc2mnvb`, `regmoms`, `qP3sub`,
`detrat`, and `var_emab`. Rewrite the header comment — the current text says the exposed
routines are the ones R can already call, and that is no longer the selection rule.

**Traps.**
- `mn2mvarb(mn, s_mn, mc, s_mc)` takes `mc` as `intent(in,out)`: `mn2m_var` writes it.
- `var_emab`'s `cv_yp_syp(2,2,*)` is 3-D with the quantile index last; declare it
  `dimension(2,2,nq)`.
- `mse_ema` pseudo-orthogonalises internally (`mc(1) = 0`, thresholds shifted by `mc(1)`).
  Pass raw at-site moments and raw thresholds; do not pre-shift, or the shift is applied twice.
- Threshold grouping must be on **exact** values. `_threshold_groups` in
  `test_fortran_oracles.py` already documents why rounding to 12 decimals moves Big Sandy's
  skew MSE by 2.2e-4.

**Acceptance.** One oracle test per new symbol pinning its value on all three parity cases,
`make parity` green. Nothing in `hydrolib/` changes in this phase.

### Phase 2 — P3 primitives (1–2 days)

**New module** `hydrolib/ema/p3.py`. Pure functions, no analyzer state.

| function | Fortran | notes |
|---|---|---|
| `m2p`, `p2m`, `p2mn`, `m2mn`, `mn2m`, `mn2p` | `probfun.f:752-838` | trivial algebra |
| `mp3(tl, tu, m, k)` | `emafit.f:2983` | **generalise the existing code to k ≤ 6** |
| `pp3`, `dp3`, `qp3` | `emafit.f:3210/3154/3266` | Wilson–Hilferty blend |
| `whlp2z`, `whz2lp`, `choose` | `emafit.f:3127/3141` | |

`mp3` is not new work. `_compute_ema_moments`'s truncated-gamma path is `mp3` for k ≤ 3 and
§1.1 shows it is exact; lift it out, generalise the binomial expansion to k = 1..6, and have
`_compute_ema_moments` call it. `var_mom` needs k = 6.

**Traps.**
- The gamma/Wilson–Hilferty blend weight is `w = max(0, (|g| − 0.0007)/(0.0010 − 0.0007))`,
  `wwh = (1 + cos(πw))/2` for `w < 1`, else 0. hydrolib's current hard switch at
  |g| = 0.001 agrees *above* the band but not inside it. Port the cosine taper.
- Off-support intervals: `mP3` returns `t**k` at whichever endpoint is nearer the mean when
  `cdfu == cdfl`. hydrolib uses the midpoint. Different, and only reachable on pathological
  input — port the Fortran's rule anyway.
- `qP3` has explicit `q <= 0` / `q >= 1` branches returning `parms(1)` (the bound τ) or
  ±1e31 depending on the sign of the skew. `log_pearson3_ppf` does not do this.

**Oracle.** `qP3sub` from Phase 1 for `qp3`; `moms_p3` for `mp3` at k ≤ 3. For k = 4..6
there is no direct oracle — check against `scipy.integrate.quad` of `x**k · dp3(x)` over the
interval, which is independent of the transcription.

**Acceptance.** `mp3` k ≤ 3 within 1e-12 of the `moms_p3` back-out on every parity interval;
k = 4..6 within 1e-8 of quadrature; `qp3` within 1e-12 of `qP3sub` across the standard AEPs
and both skew signs.

### Phase 3 — `var_mom` (2–3 days)

**New module** `hydrolib/ema/varmom.py`, transcribing `emafit.f:2299-2513`.

```
var_mom(nthresh, n, tl, tu, mc) -> (3,3)
  ├── clamp mc: mc[0]=0, |g| into [0.06324555, 1.41] keeping sign
  ├── per threshold group: pP3 at tl,tu -> p1,p2,p3
  │                        mP3 three times (k=6) -> mu_x, e_x
  │                        varb -> mu_x @ vn @ mu_x.T          (23 lines)
  │                        varc -> mu_n*(e_x[i+j] - e_x[i]e_x[j])  (23 lines)
  │                        d_est -> expmomderiv at qP3(pa), qP3(pb)  (39 lines)
  └── a = inv(I - d_t/n_t); varm = (a @ (vb_t+vc_t) @ a.T)/n_t**2
```

`d_est` needs `expmomderiv` = `dexpect ∘ dpdm` (`probfun.f:841/1136/648`), which in turn
needs `ddgam` — ∂P(α,x)/∂α (`probfun.f:1038`).

**Traps.**
- `ddgam` **deliberately returns 0** when `|x − α|/√(α+1) > 7`. Reproduce the cutoff; a more
  accurate implementation will not match.
- `mc[0]` is forced to 0 and thresholds shifted by `mc_in[0]` — but the `tl > tu` check
  (`emafit.f:2350`) uses the *shifted* values and hard-errors. Raise, don't clamp: it fires
  when MGBT lifts `tl` above an input `tu`, which is real bad input.
- `dmxytf(A,B)` is `A @ B.T`; `dmxtyf(A,B)` is `A.T @ B`. Two characters apart, used in
  adjacent lines of `gridmake`/`var_mom`.
- `d_est` skips a tail whose probability is below `pmin = 1e-6`, leaving that Jacobian zero.
- Fortran arrays are column-major; `cmoms(3,3)`, `varm(3,3)` and `mu_x(3,3)` all come back
  transposed relative to the naive read. `mu_x[j][col]` is `col ∈ {below, inside, above}`.

**Oracle.** `var_mom`, already wired.

**Acceptance.** Within 1e-10 relative of the Fortran on all three parity cases, plus a
random sweep over skews in [−1.5, 1.5], 1–4 threshold groups and censoring fractions
0–90%. Retain the existing symmetry / positive-definiteness assertions.

### Phase 4 — `mse_ema`, `mn2mvarb`, and ADJE (3–4 days)

The hard phase. `mn2mvarb` inverts `mc2mnvb` by Newton iteration with a finite-difference
Jacobian and step-halving.

```
mn2mvarb(mn, s_mn, mc) -> s_mc            emafit.f:2514, 77 lines
  ├── mn2m_var          first-order start                    emafit.f:2864
  ├── mc2mnvb           s_mc -> s_mn by 2x2x2 quadrature      emafit.f:2594
  │     └── gridmake    normquad + gammaquad + chol33         emafit.f:2039
  ├── jmc2mnvb          6x6 Jacobian by central difference    emafit.f:2654
  │     └── diff2/jacf  step shrinks x10 until chol33 succeeds
  └── Newton on the 6 unique elements, halving until chol33 succeeds
```

Then `mse_ema` = `s_mc[kmom,kmom]`, and `mseg_all` under ADJE is
`bias_adj * mseg(min(n,150), g)` with
`bias_adj = mse_ema(actual groups) / mse_ema(one group of min(n,150), thresholds ±99)`.

**Quadrature.** `normquad(n)` is Gauss–Hermite rescaled by √2 and 1/√π —
`numpy.polynomial.hermite.hermgauss` gives it directly. `gammaquad(n, α, β)` is generalised
Gauss–Laguerre with parameter `α−1`, nodes scaled by β, weights by `1/Γ(α)`;
`scipy.special.roots_genlaguerre` gives it, and for `α > 160` the Fortran caps α and
rescales β — port that branch.

**Traps.**
- `chol33` is *not* Cholesky. It is upper-triangular only after swapping rows 1↔2 and
  columns 1↔2, because the middle variate is gamma and the outer two normal. Transcribe it
  literally; `numpy.linalg.cholesky` is a different matrix and will silently give a wrong
  quadrature grid.
- The convergence test is `err = D(6)**2 <= 1e-8` **and** `I == 1` — converged *and* the
  last step was a full step. Both conditions.
- `jmc2mnvb`'s step is `delta * s**sterm(i)` with `sterm = (2,3,1,4,2,0)` and
  `delta = 1e-4`, scaling each covariance element by its own power of σ. Not a uniform step.
- `diff2` calls `jacf` three times per `(i,j)`; `jacf` for fixed `i` and varying `j` returns
  elements of the *same* `s_mn`. Compute each perturbation once and fill the whole column —
  18 `mc2mnvb` evaluations per Jacobian, not 108.
- `mseg_all`'s second `mse_ema` call passes `nthresh = 1` with `INF = (−99, 99)` and
  `n_adj = min(n, 150)`, while the first passes the real groups with the real `n`. The
  asymmetry is deliberate and commented in the Fortran (SAS 12/20/2021).
- The `bias_adj < 1` clamp is commented out upstream (SAS 9/13/2023). Leave it out.
- `at_site_option` is `ADJE` **except** when MGBT ran and found PILFs, where `emafitpr:707`
  sets `B17B`. So Cains Coulee — 11 PILFs from MGBT — takes the **B17B** path and Big Sandy
  takes ADJE. Getting this backwards makes Cains Coulee look broken.

**Oracle.** `mn2mvarb`, `mc2mnvb`, `mse_ema` from Phase 1; `mseg_all_sub` already wired.

**Acceptance.** `mse_ema` within 1e-9 relative for kmom = 1,2,3 on all three cases;
`mseg_all` reproducing `as_G_mse_o` (0.09437 Big Sandy, 0.22119 Cains Coulee, 0.07154 Powder
River) within 1e-9. `_b17b_skew_mse`'s n > 150 divergence — 31% high at n = 200 — closes,
and `test_b17b_formula_diverges_above_150_observations` inverts into an equality test.
`as_G_PRL_o` (54.373 Big Sandy) follows from two more `mseg_all` calls.

### Phase 5 — `detrat` and `Wd` (1–2 days)

`detrat` (`emafit.f:3615`, 140 lines) is self-contained given Phase 2 and `expmomcderiv`
(`probfun.f:882`), which is `expmomderiv` plus Greg Schwarz's EQ 23 central-moment
conversion. Wire into `_regional_skew_equivalent_years`, replacing the `wd = 1.0` and its
`logger.debug`. Add the `ERL` option (`erlg`, `emafit.f:1781`) while there; it is 20 lines
and completes `wght_opt_n`.

**Trap.** `detrat` takes **post-MGBT** thresholds. Cains Coulee's input thresholds carry no
censoring and give `Wd = 1.0`; only after MGBT lifts the lower threshold to log10(332) does
it give 0.184. `TestDeterminantRatioOracle::test_needs_the_post_mgbt_thresholds` already
pins this. The native path must group thresholds *after* MGBT censoring, from
`self._intervals`, not from the constructor arguments.

**Oracle.** `detratsub`, already wired.

**Acceptance.** `Wd` within 1e-9 of `detratsub` on all three cases. With Phases 0/4/5 in
place the weighted skew should close: **this is the phase where the two remaining
`xfail(strict=True)` skew assertions — Big Sandy's `test_weighted_skew` and Cains Coulee's
`skew_weighted` — are expected to start passing.** TODO.md records that feeding peakfq's own
`as_G_mse` (0.09437) through the current code already gives −0.1592 against −0.1563, a 1.9%
gap inside the 0.02 bound, so Big Sandy is bounded by ADJE alone; Cains Coulee needs `Wd`.

### Phase 6 — `VAR_EMAB` and CI shape (2–3 days)

```
var_emab(nth, nobs, tl, tu, mc, pq, eps, r_s2, r_m_mse, r_s2_mse, r_g_mse)
  ├── regmoms -> s_mc      var_mom + m2mn + mn2mvarb, then regional blending
  ├── gridmake -> 8 outer nodes; regmoms+gridmake again -> 8x8 inner nodes
  ├── per quantile: qP3 on every node
  │     cv_yp_syp = [[cov(yp,yp), cov(yp,syp)], [.., cov(syp,syp)]]
  └── ci_ema_m3b -> beta1, nu, t, asymmetric bounds
```

`ci_ema_m3b` is 50 lines and TODO.md already transcribes it. `regmoms` is 110 lines of
regional blending on top of Phase 3 and 4 output. The quadrature is 8 outer × 8 inner = 64
`regmoms` calls per fit — each a `mn2mvarb` Newton solve, so this is the phase to profile.

Add `compute_confidence_limits(..., method=...)`, defaulting to the new path with the
existing symmetric formula retained and reachable. `run_analysis` should carry
`var_est` (= `cv_yp_syp[0,0]`) and `nu` per quantile onto `FrequencyResults`.

**Traps.**
- `fp_tnc_icdf` inverts `fp_tnc_cdf` by bisection on [−32, 32], and `fp_tnc_cdf` uses a
  **normal approximation for ν > 20**, not the exact t. `scipy.stats.t.ppf` will not match;
  `nu` from `ci_ema_m3b` frequently exceeds 20. Port the branch.
- `ci_ema_m3b` returns `nu = cv_yp_syp(1,1)` on its last line — it overwrites the degrees of
  freedom with the variance. That is how `emafitpr` gets `var_est`. Deliberate, keep it.
- `nu` floors at 5 and the denominators floor at `c_min = 0.5`.
- `emafitpr` interpolates two `var_emab` calls at ±0.06324555 when |skew| is below that
  (`emafitpr:826-852`), and monotonicity is enforced on the bounds afterwards.
- `regmoms`'s `wG` calls `mseg_all` rather than using `s_mc(3,3)` — the vendored source
  carries a comment saying this makes no sense but that changing it moves the intervals.
  Port what is there.

**Acceptance.** Asymmetry ratio within 0.05 of the reference at AEP 0.1/0.02/0.01 —
1.03 / 1.31 / 1.41 on Big Sandy — retiring the last two `xfail(strict=True)`.
`test_native_intervals_are_exactly_symmetric` must be **deleted**, not inverted: it pins
behaviour this phase removes. Add a runtime budget assertion; a single `run_analysis()`
costs about two seconds today.

### Phase 7 — the small follow-ups

Independent of everything above, and each already specified in TODO.md:

- `FrequencyComparator` comparing skew by percent difference. Give skew an absolute
  tolerance in skew units instead of leaning on `parameter_scale_floor`.
- `plot_peak_flows_with_thresholds` unused — move the app's return-period lines and
  max-peak annotation into `hydrolib/freq_plot.py`, delete the app copy.
- `MethodOfMoments` dropping `user_low_outlier_threshold`. Real B17B work (censoring plus the
  conditional-probability adjustment); do it after Phase 6 or explicitly decline it.
- `origin/dev` deletion is the owner's call. Tip `86cb147`, recorded for recovery.
- Tag pushes still 403 from the agent environment.

---

## 3. Sequencing and risk

```
Phase 0  bcf ────────────────┐
Phase 1  oracles ────────────┤
Phase 2  P3 primitives ──────┴──> Phase 3 var_mom ──> Phase 4 mse_ema/ADJE ──┬──> Phase 6 CI
                                                   └─> Phase 5 detrat/Wd ────┘
```

Phases 0 and 1 are independent of each other and of everything else; do both first. Phase 5
depends on Phase 2 but not on Phases 3-4, so it can run in parallel with Phase 4 — and
`detratsub` is already wired, so it needs nothing from Phase 1.

**Where this is most likely to go wrong:** Phase 4's Newton solver. It is the only iterative
piece with no closed form, its convergence test has a subtle second condition, and its
Jacobian depends on a non-standard Cholesky. Budget for it going long. If it does, Phases 0,
2, 3 and 5 still land the skew work; only ADJE and the CI shape are gated on it.

**What could invalidate the plan:** nothing in Phases 0-3, which are measured or transcribed
from readable source. Phase 6's cost is the real unknown — 64 `mn2mvarb` solves per fit is a
guess at the shape of the work, not a measurement, and if it is too slow the answer is to
cache `regmoms` across quantiles (`s_mc` does not depend on `pq`), which the Fortran does not
bother to do.

## 4. Conventions

- Every phase ends with `black`, `isort`, `pytest tests/` and `make parity` green.
- Reproduce CI: `PYTHONSAFEPATH=1 python -m pytest`, and remember the matrix is 3.9–3.12.
- Known-failing stays `xfail(strict=True)`. When a phase retires one, delete the marker in
  the same commit — strict means the build goes red the moment it starts passing.
- Do not edit anything under `vendor/`. Every oracle in this plan is reachable by declaring
  the symbol in `_emafort.pyf`; that was verified, not assumed.
- Golden files are regenerated only by `python tools/gen_fortran_golden.py`, never by hand.
