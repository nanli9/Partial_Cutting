# Nodal ring: fast per-tet method vs solving the whole ring directly

Date: 2026-11-03. Machine: Intel Core Ultra 9 386H, 16 threads, Windows 11, MSVC Release build of 2026-11-03 18:39:45.
All numbers come from runs in `ring_perf/` (files starting with `final_`), checked by an independent re-run (files starting with `audit5_`).
Anything not measured is marked **derived** or **estimate**.

---

## 1. The answer

1. **Faster than before.** On ring 2, crack K1, our method used to take 69.8 s per ring at 1.400% error. It now takes 0.295 s at 1.394% error. That is **237x faster at the same accuracy**. The cheapest setting takes 0.114 s at 1.94% error (612x). The goal (10x faster, error at most about 2%) is met.
2. **Against solving the whole ring directly** (same solver, same points, same threads): to reach a value error of 3%, 2% and 1.5%, ours needs 0.11, 0.11 and 0.19 s. The direct solve needs 2.5, 5.7 and 10.3 s. Ours is **22-54x faster per ring**.
3. **Ours has a floor.** Its error does not go below about 1.4% (K1), 1.8% (K2) and 1.1% (K3), however many walks it runs. This floor comes from the glue, not from noise. The direct solve has no floor. It passes ours after about **12 s per ring** and keeps improving (0.83% at 36 s).
4. **Worst-point error and vertex values are a different story.** On the largest single-point error and on the 7 vertex values, the direct solve matches or beats ours much sooner: after about 1-4 s per ring for K1, 0.5-2 s for the ramp data, 8-43 s for K2.
5. **Gradient.** Ours has a much smaller point-by-point gradient error at every time we measured (31% vs 52%, even when the direct gets 36 s). But ours' gradient error has a floor of 21-30%. The direct one does not.
6. **K3 (crack through the whole ring):** the direct solve is exact at any cost. Ours keeps a 1.1% error (vertex v3 = 0.773 instead of 1).
7. **Verdict.** If 1.4-1.8% average error is good enough, ours is far cheaper: 22-54x per ring, and about 87-217x per mesh node (**derived**). If you need better than that, or good worst-point and vertex values, the direct solve wins after a few seconds to tens of seconds per ring.

---

## 2. Words used in this report

- **Walk**: one random walk from a sample point. The solver averages many walks per point. More walks = less noise, more time.
- **Noise**: the random part of the error. It shrinks when you run more walks.
- **Bias** (or **floor**): the part of the error that stays when you run more walks.
- **GT** (reference): one very long direct run on the whole ring, 25,600 walks per point. We use it as the "true" answer.
- **Direct method**: one walk solve of the whole ring O in physical space, like the GT, but with fewer walks and different random seeds.
- **Our method**: solve each tet on its own (in its own reference-tet coordinates xi), then glue the 8 tets with one 7x7 linear system for the 7 vertex values.
- **Given-value face** (Dirichlet face): a face where the solution's value is prescribed. All faces in our per-tet solves are like this.
- **No-flow face** (Neumann face): a face where only "nothing flows through" is prescribed. The outer faces of O in the direct solve and the GT are like this. Walks have to bounce off them, so their walks are longer.
- **Corner function P_k**: the per-tet solution with value 0 on the crack and the linear "hat" lambda_k of corner k on the faces.
- **H0**: the per-tet solution with value +-g on the crack and the formula J on the faces (Level C).
- **Gram table**: the 4x4 table of energy products between the corner functions, plus 4 products with H0. The glue needs only these numbers.
- **Control variate**: if a known exact function is close to the answer, solve only the difference by walks and add the known function back exactly. The difference is small, so it needs far fewer walks. The result is still correct on average.
- **Q_k and R**: the two differences we solve by walks. Q_k = P_k - lambda_k. R = H0 - J_xi, where J_xi is J built in the tet's own xi coordinates (exactly harmonic in xi).
- **Stop distance** (shell): a walk stops when it gets this close to the boundary. Smaller = more accurate, more steps.
- **H RMS**: the average error of the glued H over O (weighted by volume), divided by 2g, in %.
- **H nc**: the same, after removing the noise of both the result and the GT. What is left is the bias.
- **H max**: the largest error at any single point, divided by 2g.
- **grad raw**: average gradient error, relative to the GT's average gradient size. It includes the GT's own 15-26% gradient noise.
- **E err**: error of the energy (1/2) int |grad H|^2 over the ring (minus a thin tube around the crack), in %. **E lines**: the same, only within 0.1 of the crack-face lines.
- **vert max**: the largest error over the 7 vertex values.
- **Ramp**: crack data g = min(1, sqrt(r/L)) instead of g = 1. It is smoother, so walks are less noisy.
- **Setting names**: "Q32/R16/B8" = 32 walks per point for Q, 16 for R, in 8 batches. "N256/B16" = 256 walks per point, 16 batches.

---

## 3. What was changed to make it faster

All changes are opt-in. Check 4's existing modes are untouched. The per-tet solve stays in each tet's own xi coordinates.

| # | Change | Measured gain | Changes accuracy? |
|---|---|---|---|
| 1 | **Lean walk step.** Drop the walk's bounding box (never read). When every face has given values, also drop the face ray tests (they can never hit). | Ours 2.09x (0.618 s -> 0.295 s). Direct 1.52x (5.404 s -> 3.548 s), which keeps its face tests because O's faces are no-flow. | No. Bit-identical walks (checked with memcmp). |
| 2 | **Control variate for the corner functions:** P_k = lambda_k + Q_k. Walks find only Q_k; lambda_k and its gradient are added exactly. | Part of the combined gain below. It is what allows 32 walks per point instead of 2560. | No (correct on average; z tests pass). |
| 3 | **Control variate for H0:** H0 = J_xi + R. Walks find only R. | Part of the combined gain. In ring 2 K1, the R walks take about 0.09-0.10 s of the 0.295 s. Ring 1 needs no R walks at all. | No (correct on average; z tests pass). |
| 4 | **Skip fields that are known exactly.** No Q walks in uncut tets (P_k = lambda_k there). No R walks where J_xi equals J (all of ring 1). | Q walks run in 7/4/7 of 8 tets in ring 2 (K1/K2/K3) and 6/3/7 in ring 1. R walks run in 8 tets in ring 2, 0 in ring 1. | No. |
| 5 | **One walk run per tet instead of two.** The two Level-A fields that nothing used are no longer computed. | Together with #2-#4: 1.43 M walks per ring instead of 166 M (116x fewer). | No. |
| 6 | **Fewer walks per point:** 2560 in 64 batches -> Q 32 / R 16 in 8 batches. Allowed by #2 and #3. | Included in the 116x above. | Only the noise. The floor stays (1.394% vs 1.400% before). |
| 7 | **All tets in one parallel loop** (fused). | None measurable (0.271 s vs 0.272 s). | No. Kept because it is free. |
| 8 | **Looser stop distance** 1e-4 sqrt(2) instead of 1e-6 sqrt(2). Opt-in (`--shell`), not used in any headline number. | 29-37% faster, both methods (walks 34 -> 19 steps for ours, 46 -> 31 for the direct). | Slightly biased. No visible H change (1.388% vs 1.394%). |

**Total, ring 2 K1:** 69.8 s -> 0.295 s (237x). Rough split: 2.1x from the lean step, the rest from running 116x fewer walks.
**Cutting walks alone does not do it:** check 4's old path at 256 walks takes 7.1 s (10x faster), but its gradient error grows to 84%.

---

## 4. Comparison tables

How everything was measured:

1. Wall clock on 16 threads. Median of 3 repeats; runs over 30 s ran once.
2. Same evaluation points for every method: check 4's points, about 32,000 per ring. Both H and grad H are computed there.
3. Same solver code and same stop distance (1e-6 sqrt(2)) for both methods. Same lean step where it applies.
4. Each method is charged for all its walks. The direct also pays for 7 extra vertex points. Ours gets its vertex values from the glue, without walks.
5. The GT has 25,600 walks per point, at least 10x any compared run, with independent random seeds.
6. Background load from other programs was 14-32% (Chrome, Claude). No other partial_cut.exe was running.

### 4.1 Equal time: best average value error (H RMS, %) each method reaches within a time budget

"-" = no measured setting fits the budget. This is **average value error only**. See 4.3 for worst-point and vertex errors.

| config | method | <= 0.25 s | <= 1 s | <= 5 s | <= 20 s | <= 70 s |
|---|---|---|---|---|---|---|
| ring 2 K1 | ours | 1.442 (Q16/R16, 0.205 s) | 1.382 (Q32/R64, 0.571 s) | 1.382 | 1.382 | 1.382 |
| ring 2 K1 | direct | 9.976 (N16, 0.228 s) | 4.967 (N64) | 2.500 (N256) | 1.280 (N1024, 14.3 s) | 0.828 (N2560, 36.2 s) |
| ring 2 K1 | old path (check 4) | | | | | 1.400 (69.8 s) |
| ring 2 K2 | ours | 1.819 (Q8/R8, 0.085 s) | 1.806 (Q128/R64, 0.73 s) | 1.805 | 1.805 | 1.805 |
| ring 2 K2 | direct | - (N16 takes 0.268 s: 12.09) | 8.571 (N32) | 3.035 (N256) | 1.532 (N1024, 17.0 s) | 0.994 (N2560, 43.5 s) |
| ring 2 K2 | old path (check 4) | | | | | 1.820 (64.3 s) |
| ring 2 K3 | ours | 1.105-1.110 | same | same | same | same |
| ring 2 K3 | direct | 0.000 (exact from 0.196 s) | 0.000 | 0.000 | 0.000 | 0.000 |
| ring 1 K1 | ours | 1.367 | 1.367 | 1.367 | 1.367 | 1.367 |
| ring 1 K1 | direct | - | 4.981 | 2.496 | 1.280 | 0.832 |
| ring 1 K2 | ours | 1.759 | 1.743 | 1.743 | 1.743 | 1.742 |
| ring 1 K2 | direct | - | 8.512 | 3.038 | 1.544 | 1.005 |
| ring 1 K1 ramp | ours | 1.162 | 1.158 | 1.158 | 1.158 | 1.158 |
| ring 1 K1 ramp | direct | - | 1.726 | 0.859 | 0.441 | 0.288 |
| ring 1 K2 ramp | ours | 1.529 | 1.529 | 1.529 | 1.529 | 1.529 |
| ring 1 K2 ramp | direct | - | 2.209 | 0.791 | 0.399 | 0.260 |

Notes:
- Ours' error after removing noise is only 0.03-0.1 points lower than shown: ours' error is almost all floor.
- The direct's error after removing noise is 0-0.3%: its error is almost all noise.
- Each value is the best setting of the grid, from one set of random seeds. Across 3 seed sets the values move by about +-0.05 points (section 4.6). This favours ours slightly at small walk counts. It changes no ranking.

**When the direct catches up on average value error** (time where the direct's log-log interpolated curve reaches ours' floor):

| config | ours' floor | direct reaches it at |
|---|---|---|
| ring 2 K1 | 1.39% | about 12.0 s |
| ring 2 K2 | 1.81% | about 12.1 s |
| ring 1 K1 | 1.37% | about 12.6 s |
| ring 1 K2 | 1.74% | about 13.2 s |
| ring 1 K1 ramp | 1.16% | about 2.0 s |
| ring 1 K2 ramp | 1.53% | about 1.1 s |
| K3 (both rings) | 1.1% | always (direct is exact) |

### 4.2 Equal accuracy: time to reach an average value error (H RMS)

First measured setting at or below the threshold, with the log-log interpolated time in [brackets].

| config | method | 3% | 2% | 1.5% | 1% |
|---|---|---|---|---|---|
| ring 2 K1 | ours | 0.114 s (Q8/R8/B4) | 0.114 s | 0.205 s [0.190] | not reached (floor 1.36%) |
| ring 2 K1 | direct | 3.548 s [2.476] | 7.247 s [5.667] | 14.336 s [10.321] | 36.160 s [24.324] |
| ring 2 K1 | old path at 256 walks | 7.146 s | 7.146 s | 7.146 s | not reached |
| ring 2 K1 | old path (check 4) | 69.8 s | 69.8 s | 69.8 s | not reached |
| ring 2 K2 | ours | 0.085 s | 0.085 s | not reached (floor 1.73-1.78%) | not reached |
| ring 2 K2 | direct | 8.464 s [4.379] | 16.990 s [9.839] | 43.481 s [18.377] | 43.481 s [42.959] |
| ring 2 K2 | old path at 256 walks | 7.103 s | 7.103 s | not reached | not reached |
| ring 2 K3 | ours | 0.101 s | 0.101 s | 0.101 s | not reached (floor 1.10%) |
| ring 2 K3 | direct | 0.196 s (exact) | 0.196 s | 0.196 s | 0.196 s |
| ring 1 K1 | ours | 0.045 s | 0.045 s | 0.045 s | not reached |
| ring 1 K1 | direct | 3.620 [2.490] | 7.237 [5.763] | 14.391 [10.426] | 36.131 [24.604] |
| ring 1 K2 | ours | 0.027 s | 0.051 s [0.041] | not reached (floor 1.71%) | not reached |
| ring 1 K2 | direct | 8.408 [4.377] | 16.825 [9.929] | 42.674 [18.096] | not reached (1.004% at N2560) |
| ring 1 K3 | ours | 0.052 s | 0.052 s | 0.052 s | not reached |
| ring 1 K3 | direct | 0.194 s (exact) | 0.194 s | 0.194 s | 0.194 s |
| ring 1 K1 ramp | ours | 0.050 s | 0.050 s | 0.050 s | not reached (floor 1.16%) |
| ring 1 K1 ramp | direct | 0.468 [0.321] | 0.912 [0.682] | 1.799 [1.195] | 3.620 [2.660] |
| ring 1 K2 ramp | ours | 0.030 s | 0.030 s | not reached (floor 1.53%) | not reached |
| ring 1 K2 ramp | direct | 0.531 [0.293] | 1.052 [0.648] | 2.119 [1.159] | 4.269 [2.651] |

Speed-up of ours where ours reaches the threshold (interpolated times):
- ring 2 K1: 3% 2.476 / 0.114 = 22x; 2% 5.667 / 0.114 = 50x; 1.5% 10.321 / 0.190 = 54x.
- ring 2 K1 at ours' default Q32/R16 accuracy (1.394%): about 12.0 s vs 0.295 s = 41x.
- ring 2 K2: 2% 9.839 / 0.085 = 116x.
- ring 1 K1 ramp: 1.5% 1.195 / 0.050 = 24x.

Ours' times here are upper bounds: its cheapest grid point (Q8/R8/B4) already meets most thresholds.

### 4.3 Worst-point error and vertex values (the direct catches up sooner here)

Ours' worst-point error (11-15%) and vertex errors (0.06-0.23) are floors: more walks do not lower them.

| config | ours' best H max | direct matches or beats it at | ours' best vert max | direct matches or beats it at |
|---|---|---|---|---|
| ring 2 K1 | 13.82% | N256: 11.45%, 3.55 s | 0.1172 (worst: v5) | N256: 0.1089, 3.5 s; N512: 0.074, 7.2 s |
| ring 2 K2 | 11.37% | N512: 9.82%, 8.46 s | 0.0559 | N1024: 0.0577, 17 s (about equal); N2560: 0.0271, 43.5 s |
| ring 1 K1 | 13.34% | N256: 11.73%, 3.62 s | 0.0992 | N128: 0.055, 1.8 s |
| ring 1 K2 | 12.03% | N512: 8.84%, 8.41 s | 0.0634 | N1024: 0.0454, 17.0 s |
| ring 1 K1 ramp | 6.47% | N128: 5.2%, 1.80 s | 0.079 | N64: 0.053, 0.91 s |
| ring 1 K2 ramp | 7.18% | N64: 6.7%, 1.05 s | 0.1469 | N32: 0.110, 0.53 s |
| K3 (both) | 9.9-10.7% | always (exact) | 0.226-0.227 (v3 = 0.773 vs 1) | always (exact) |

Vertex values in ring 2 K1, g = 1 (v0..v6):

| source | v0 | v1 | v2 | v3 | v4 | v5 | v6 |
|---|---|---|---|---|---|---|---|
| GT | +0.984 | +0.702 | -0.448 | +0.055 | -0.830 | +0.893 | -0.531 |
| ours Q32/R16 | +0.970 | +0.680 | -0.445 | +0.004 | -0.815 | +0.776 | -0.444 |
| ours Q2560/R640 | +0.974 | +0.683 | -0.442 | +0.011 | -0.812 | +0.771 | -0.440 |
| direct N2560 | +0.984 | +0.706 | -0.433 | +0.036 | -0.840 | +0.875 | -0.521 |

### 4.4 Gradient

Raw average gradient error, relative to the GT's average gradient. It includes the GT's own 15-22% gradient noise.

| config | method | <= 1 s | <= 5 s | <= 20 s | <= 70 s |
|---|---|---|---|---|---|
| ring 2 K1 | ours | 30.9% (Q32/R64) | 27.3% | 26.7% | 26.7% |
| ring 2 K1 | direct | 352.6% | 151.6% | 82.9% | 51.6% |
| ring 2 K2 | ours | 41.0% | 37.9% | 37.5% | 37.5% |
| ring 2 K2 | direct | 706% | 227% | 111% | 74.9% |

| config | threshold | ours | direct |
|---|---|---|---|
| ring 2 K1 | 50% | 0.205 s (Q16/R16) | not reached at N2560 (51.6%, 36.2 s) |
| ring 2 K1 | 40% | 0.571 s (Q32/R64, 30.9%) | not reached |
| ring 2 K1 | 30% | 1.89 s (Q64/R256, 27.3%) | not reached |
| ring 2 K2 | 50% | 0.397 s (Q64/R32) | not reached (74.9% at N2560) |
| ring 2 K2 | 40% | 1.43 s (Q256/R128) | not reached |

What this means, in plain words:
1. Point by point, ours' gradient is much less noisy at every budget we measured.
2. But after removing noise, ours' gradient error is a floor of about 21% (K1) and 30% (K2). The direct's is about 0: it is pure noise.
3. Integrated gradient quantities show the floor. Energy near the crack-face lines ("E lines"): ours is +5.0 +- 0.2% (ring 2 K1) and +4.9-6.4% (ring 2 K2) off. The direct is consistent with 0 already at N256 in K1 (0.32 +- 3.54%), and clearly so at N1024: -0.23 +- 1.06% (14.3 s, K1) and -0.78 +- 1.57% (17 s, K2).
4. Ours' total energy is also off: +2% (K1), +9% (K2), +13% (ramp K1) and +36% (ramp K2).
5. Extrapolation (**estimate**): the direct's raw gradient error falls like 1/sqrt(N). Matching ours' 30.9% (0.571 s) would need about N = 8,600 walks, about 120 s per ring.

### 4.5 Speed-up against check 4's old path at its own setting

Old path = 4096 points per tet, 2 runs per tet, 2560 walks per point in 64 batches, check 4's seeds. Measured here on all three ring-2 cracks.

| config | old path | fast path, same H error | speed-up | fast path, same or better gradient | speed-up |
|---|---|---|---|---|---|
| ring 2 K1 | 69.81 s, H 1.400%, grad 37.55%, E +2.47 +- 0.83%, vert max 0.1218 | Q32/R16/B8: 0.295 s, H 1.394%, E +2.00 +- 0.72%, vert max 0.1172 | 237x | Q32/R64/B8: 0.571 s, grad 30.9% | 122x |
| ring 2 K2 | 64.32 s, H 1.820%, grad 47.39%, E +9.36 +- 0.99%, vert max 0.0637 | Q64/R32: 0.397 s, H 1.809% | 162x | Q32/R64: 0.438 s, grad 41.0% | 147x |
| ring 2 K3 | 66.57 s, H 1.111%, abs grad 7.84%, vert max 0.2271 | Q32/R16: 0.274 s, H 1.105% | 243x | Q32/R64: 0.531 s, abs grad 7.57% | 125x |

Cheapest settings: ring 2 K1 Q8/R8/B4 0.114 s at 1.938% (612x); ring 2 K2 Q8/R8 0.085 s at 1.819% (757x); ring 2 K3 Q8/R8 0.101 s at 1.110% (659x).

### 4.6 Random-seed spread (H RMS %, 3 independent seed sets)

| config | ours Q32/R16 | ours Q64/R32 | ours Q256/R128 | direct N256 | direct N1024 |
|---|---|---|---|---|---|
| ring 2 K1 | 1.394 / 1.395 / 1.414 | 1.408 / 1.386 / 1.389 | 1.397 / 1.385 / 1.384 | 2.500 / 2.482 / 2.536 | 1.280 / 1.259 / 1.286 |
| ring 2 K2 | 1.831 / 1.817 / 1.863 | 1.809 / 1.817 / 1.816 | 1.808 / 1.816 / 1.812 | 3.035 / 3.037 / 3.023 | 1.532 / 1.525 / 1.538 |
| ring 2 K3 | 1.105 / 1.110 / 1.111 | - | 1.109 / 1.112 / 1.110 | 0 / 0 / 0 | 0 / 0 / 0 |
| ring 1 K1 | 1.367 / 1.386 / 1.463 | 1.404 / 1.371 / 1.404 | 1.406 / 1.376 / 1.389 | 2.496 / 2.509 / 2.513 | 1.280 / 1.285 / 1.277 |
| ring 1 K2 | 1.759 / 1.751 / 1.754 | 1.786 / 1.756 / 1.747 | 1.743 / 1.746 / 1.741 | 3.038 / 3.025 / 3.022 | 1.544 / 1.540 / 1.542 |

Timings across seed sets: ours Q32/R16 ring 2 K1 0.29-0.33 s; direct N1024 14.4-14.8 s.

### 4.7 The GT's own noise (how exact the reference is)

| config | H noise, average (max), / 2g | gradient noise / gradient size | energy (ring minus crack tube) | energy near crack-face lines |
|---|---|---|---|---|
| ring 2 K1 | 0.250% (0.313%) | 15.30% | 1.51176 +- 0.00746 (0.49%) | 0.47591 +- 0.00091 |
| ring 2 K2 | 0.301% (0.313%) | 21.98% | 0.89436 +- 0.00525 | 0.27550 +- 0.00092 |
| ring 2 K3 | 0 (exact +-1) | not meaningful (gradient is 0) | 0 | 0 |
| ring 1 K1, g = 1 | 0.250% (0.313%) | 16.82% | 1.51485 +- 0.01243 | 0.73121 +- 0.00114 |
| ring 1 K1, ramp | 0.086% (0.101%) | 15.75% | 0.21675 +- 0.00250 | 0.06687 +- 0.00007 |
| ring 1 K2, g = 1 | 0.301% (0.313%) | 25.99% | 0.91960 +- 0.04610 (5%) | 0.30924 +- 0.00089 |
| ring 1 K2, ramp | 0.078% (0.083%) | 21.73% | 0.08774 +- 0.00251 | 0.02971 +- 0.00004 |
| ring 1 K3 | 0 | not meaningful | 0 | 0 |

Every "nc" column removes both the GT's noise and the result's own noise.

### 4.8 Full sweeps per configuration

Columns: time (for ours: Q walks / R walks / Gram / field), walks in millions, H RMS, H nc, H max, gradient raw / nc, E err (+- its noise), E lines (+- its noise), vert max. The glue takes 13-200 us in every row and setup 0.7-2.2 ms, so they are not shown.

#### Ring 2 K1, g = 1

| setting | time s (Q / R / Gram / field) | walks M | H RMS % | H nc % | H max % | grad raw / nc % | E err % | E lines % | vert max |
|---|---|---|---|---|---|---|---|---|---|
| ours Q8/R8/B4 | 0.114 (0.050/0.050/0.0053/0.0076) | 0.49 | 1.938 | 1.838 | 15.39 | 58.9 / 19.1 | 1.04 (1.32) | 1.87 (2.73) | 0.1523 |
| ours Q16/R16/B8 | 0.205 (0.094/0.096/0.0086/0.0053) | 0.97 | 1.442 | 1.365 | 13.88 | 48.0 / 21.0 | 1.92 (0.75) | 3.77 (1.13) | 0.1240 |
| ours Q32/R16/B8 | 0.295 (0.181/0.098/0.0092/0.0048) | 1.43 | 1.394 | 1.359 | 13.89 | 47.8 / 21.1 | 2.00 (0.72) | 4.06 (1.12) | 0.1172 |
| ours Q32/R64/B8 | 0.571 (0.181/0.373/0.0097/0.0048) | 2.98 | 1.382 | 1.351 | 13.82 | 30.9 / 21.1 | 2.06 (0.53) | 5.28 (0.43) | 0.1205 |
| ours Q64/R32/B16 | 0.580 (0.358/0.193/0.0164/0.0097) | 2.85 | 1.408 | 1.380 | 14.01 | 37.5 / 20.9 | 1.91 (0.63) | 4.90 (1.01) | 0.1244 |
| ours Q64/R256/B16 | 1.891 (0.361/1.502/0.0178/0.0089) | 10.12 | 1.394 | 1.368 | 13.96 | 27.3 / 21.3 | 1.91 (0.50) | 5.00 (0.23) | 0.1244 |
| ours Q128/R64/B16 | 1.122 (0.713/0.378/0.0195/0.0096) | 5.71 | 1.391 | 1.365 | 14.01 | 30.5 / 21.2 | 2.18 (0.53) | 5.58 (0.45) | 0.1229 |
| ours Q256/R128/B16 | 2.167 (1.392/0.745/0.0183/0.0093) | 11.41 | 1.397 | 1.373 | 14.06 | 28.4 / 21.1 | 1.99 (0.51) | 5.06 (0.33) | 0.1229 |
| ours Q1024/R256/B16 | 7.140 (5.611/1.499/0.0180/0.0100) | 37.35 | 1.388 | 1.364 | 13.98 | 27.3 / 21.3 | 1.93 (0.50) | 5.02 (0.23) | 0.1225 |
| ours Q2560/R640/B64 | 17.900 (14.039/3.764/0.0646/0.0299) | 93.36 | 1.386 | 1.363 | 13.95 | 26.7 / 21.3 | 1.94 (0.50) | 5.12 (0.21) | 0.1225 |
| direct N16/B8 | 0.228 | 0.52 | 9.976 | 0.000 | 52.14 | 854.8 / 0.0 | noise | -31.3 (99.1) | 0.1697 |
| direct N32/B8 | 0.446 | 1.04 | 7.100 | 0.571 | 35.18 | 523.7 / 155.2 | noise | -5.17 (17.19) | 0.2022 |
| direct N64/B16 | 0.902 | 2.08 | 4.967 | 0.000 | 23.57 | 352.6 / 0.0 | noise | -8.30 (8.88) | 0.1122 |
| direct N128/B16 | 1.798 | 4.15 | 3.529 | 0.000 | 15.51 | 235.8 / 46.0 | noise | -6.12 (5.65) | 0.1245 |
| direct N256/B16 | 3.548 | 8.31 | 2.500 | 0.000 | 11.45 | 151.6 / 0.0 | -31.2 (31.4) | 0.32 (3.54) | 0.1089 |
| direct N512/B16 | 7.247 | 16.62 | 1.779 | 0.000 | 8.80 | 112.0 / 0.0 | -8.80 (20.23) | -0.48 (2.42) | 0.0737 |
| direct N1024/B16 | 14.336 | 33.24 | 1.280 | 0.122 | 6.06 | 82.9 / 22.3 | 3.08 (13.00) | -0.23 (1.06) | 0.0341 |
| direct N2560/B64 | 36.160 | 83.10 | 0.828 | 0.000 | 3.84 | 51.6 / 0.0 | -1.00 (3.33) | -0.06 (0.63) | 0.0191 |
| old path N256/B16 | 7.146 | 16.62 | 1.476 | 1.376 | 13.99 | 84.0 / 26.7 | 4.64 (4.31) | 15.09 (12.07) | 0.1208 |
| old path 2560/64 (check 4) | 69.812 | 166.16 | 1.400 | 1.369 | 14.00 | 37.5 / 22.4 | 2.47 (0.83) | 5.33 (0.87) | 0.1218 |

"noise" = the direct's total-energy error is dominated by noise at small N (noise from tens to thousands of %). Mean walk length: ours 34.1-34.3 steps, direct 46.2 (its no-flow faces make walks longer). The direct's field extraction takes 0.001-0.005 s; the rest is walks.

#### Ring 2 K2, g = 1

| setting | time s | walks M | H RMS % | H nc % | H max % | grad raw / nc % | E err % | E lines % | vert max |
|---|---|---|---|---|---|---|---|---|---|
| ours Q8/R8/B4 | 0.085 | 0.39 | 1.819 | 1.696 | 11.37 | 73.9 / 31.2 | 10.80 (1.19) | 7.01 (1.08) | 0.0678 |
| ours Q16/R16/B8 | 0.145 | 0.77 | 1.840 | 1.757 | 12.01 | 58.6 / 29.7 | 9.43 (0.83) | 6.42 (1.60) | 0.0564 |
| ours Q32/R16/B8 | 0.196 | 1.02 | 1.831 | 1.730 | 11.90 | 58.6 / 29.7 | 9.45 (0.80) | 6.54 (1.57) | 0.0575 |
| ours Q32/R64/B8 | 0.438 | 2.57 | 1.831 | 1.745 | 11.91 | 41.0 / 29.6 | 8.96 (0.65) | 5.48 (0.74) | 0.0559 |
| ours Q64/R32/B16 | 0.397 | 2.05 | 1.809 | 1.760 | 11.57 | 47.8 / 29.7 | 9.27 (0.71) | 6.23 (1.02) | 0.0633 |
| ours Q128/R64/B16 | 0.730 | 4.10 | 1.806 | 1.775 | 11.42 | 41.5 / 29.7 | 8.97 (0.62) | 5.37 (0.83) | 0.0668 |
| ours Q256/R128/B16 | 1.434 | 8.19 | 1.808 | 1.779 | 11.46 | 39.0 / 29.6 | 8.92 (0.63) | 4.94 (0.69) | 0.0658 |
| ours Q2560/R640/B64 | 11.029 | 61.30 | 1.807 | 1.781 | 11.46 | 37.5 / 29.7 | 9.03 (0.59) | 5.26 (0.38) | 0.0647 |
| direct N64/B16 | 1.093 | 2.07 | 6.055 | 0.586 | 29.82 | 524.0 / 145.0 | noise | 15.2 (12.7) | 0.1730 |
| direct N256/B16 | 4.279 | 8.26 | 3.035 | 0.238 | 11.48 | 227.0 / 0.0 | -25.6 (36.5) | -1.94 (4.36) | 0.0948 |
| direct N512/B16 | 8.464 | 16.52 | 2.152 | 0.097 | 9.82 | 153.6 / 0.0 | -30.3 (25.4) | -2.56 (2.54) | 0.0636 |
| direct N1024/B16 | 16.990 | 33.04 | 1.532 | 0.000 | 6.45 | 111.1 / 0.0 | -4.86 (12.31) | -0.78 (1.57) | 0.0577 |
| direct N2560/B64 | 43.481 | 82.61 | 0.994 | 0.000 | 4.55 | 74.9 / 14.2 | 2.66 (4.86) | 1.32 (1.05) | 0.0271 |
| old path N256/B16 | 7.103 | 16.52 | 1.896 | 1.792 | 11.80 | 99.7 / 26.7 | 6.80 (3.10) | -3.87 (6.44) | 0.0703 |
| old path 2560/64 (check 4) | 64.315 | 165.19 | 1.820 | 1.787 | 11.57 | 47.4 / 30.2 | 9.36 (0.99) | 4.64 (1.84) | 0.0637 |

Direct N16 / N32 / N128: 0.268 / 0.535 / 2.164 s, H RMS 12.089 / 8.571 / 4.290%. Mean walk length: ours 34.3, direct 60.0. Ours runs Q walks in 4 tets and R walks in 8.

#### Ring 2 K3, g = 1 (the crack cuts through the whole ring; the GT is exactly +-1 with zero gradient, so absolute gradient and energy are shown)

| setting | time s | walks M | H RMS % | H max % | abs grad / 2g % | E (absolute) | vert max |
|---|---|---|---|---|---|---|---|
| ours Q8/R8/B4 | 0.101 | 0.48 | 1.110 | 10.65 | 14.05 | 0.0112 | 0.2269 |
| ours Q32/R16/B8 | 0.274 | 1.40 | 1.105 | 10.62 | 11.23 | 0.0120 | 0.2262 |
| ours Q32/R64/B8 | 0.531 | 2.94 | 1.109 | 10.64 | 7.57 | 0.0117 | 0.2267 |
| ours Q256/R128/B16 | 2.064 | 11.23 | 1.109 | 10.64 | 7.17 | 0.0117 | 0.2267 |
| ours Q2560/R640/B64 | 16.684 | 91.83 | 1.110 | 10.65 | 6.83 | 0.0118 | 0.2269 |
| direct, every N (0.196 s at N16 to 30.70 s at N2560) | - | 0.51-81.9 | 0.000 | 0.00 | 0.00 | 0.0000 | 0.0000 |
| old path 2560/64 (check 4) | 66.571 | 163.70 | 1.111 | 10.66 | 7.84 | 0.0117 | 0.2271 |

Ours' glued v3 = +0.773 against an exact +1, at every setting and in the old path too.

#### Ring 1 K1, g = 1 (ring 1 needs no R walks)

| setting | time s | H RMS % | H nc % | H max % | grad raw / nc % | E err % | E lines % | vert max |
|---|---|---|---|---|---|---|---|---|
| ours Q8/R8/B4 | 0.045 | 1.389 | 1.333 | 13.34 | 26.8 / 20.2 | 2.40 (0.82) | 3.98 (0.21) | 0.0992 |
| ours Q32/R16/B8 | 0.149 | 1.367 | 1.338 | 13.55 | 26.5 / 20.4 | 2.39 (0.82) | 3.99 (0.16) | 0.1075 |
| ours Q256/R128/B16 | 1.145 | 1.406 | 1.382 | 13.70 | 26.5 / 20.4 | 2.42 (0.82) | 4.03 (0.16) | 0.1115 |
| ours Q2560/R640/B64 | 12.242 | 1.386 | 1.363 | 13.73 | 26.4 / 20.4 | 2.41 (0.82) | 4.02 (0.16) | 0.1119 |
| direct N128/B16 | 1.799 | 3.520 | ~0 | 18.0 | 251 raw | noise | -0.2 | 0.055 |
| direct N256/B16 | 3.620 | 2.496 | 0.000 | 11.73 | 179.5 / 45.6 | 17.0 (46.0) | -1.05 (2.47) | 0.0306 |
| direct N512/B16 | 7.237 | 1.794 | 0.177 | 8.73 | 121.2 / 0 | -17.4 (31.8) | 0.10 (1.11) | 0.0345 |
| direct N1024/B16 | 14.391 | 1.280 | 0.121 | 6.74 | 81.7 / 0 | -5.44 (14.31) | 0.61 (0.67) | 0.0177 |
| direct N2560/B64 | 36.131 | 0.832 | 0.072 | 3.76 | 51.9 / 0 | -3.24 (5.49) | 0.60 (0.55) | 0.0139 |
| old path N256/B16 | 6.805 | 1.458 | 1.371 | 13.51 | 81.5 / 21.6 | 3.05 (2.79) | 7.25 (4.52) | 0.1090 |

#### Ring 1 K2, g = 1

| setting | time s | H RMS % | H nc % | H max % | grad raw / nc % | E err % | E lines % | vert max |
|---|---|---|---|---|---|---|---|---|
| ours Q8/R8/B4 | 0.027 | 2.365 | 2.161 | 13.61 | 41.0 / 30.3 | 5.36 (5.02) | 3.17 (0.42) | 0.1059 |
| ours Q32/R16/B8 | 0.085 | 1.759 | 1.707 | 12.44 | 40.4 / 30.9 | 5.49 (5.01) | 3.22 (0.29) | 0.0728 |
| ours Q256/R128/B16 | 0.579 | 1.743 | 1.713 | 12.16 | 40.4 / 30.9 | 5.49 (5.01) | 3.21 (0.29) | 0.0664 |
| ours Q2560/R640/B64 | 5.705 | 1.742 | 1.715 | 12.03 | 40.4 / 30.9 | 5.49 (5.01) | 3.19 (0.29) | 0.0634 |
| direct N256/B16 | 4.269 | 3.038 | 0.264 | 14.19 | 281.5 / 85.7 | noise | 2.62 (3.44) | 0.0920 |
| direct N512/B16 | 8.408 | 2.170 | 0.292 | 8.84 | 180.7 / 0 | noise | 1.98 (1.95) | 0.0659 |
| direct N1024/B16 | 17.010 | 1.544 | 0.163 | 7.70 | 122.0 / 0 | -10.1 (72.1) | 1.11 (0.93) | 0.0454 |
| direct N2560/B64 | 42.674 | 1.004 | 0.108 | 4.82 | 78.0 / 0 | -17.6 (26.8) | 0.58 (0.85) | 0.0453 |
| old path N256/B16 | 6.485 | 1.834 | 1.712 | 12.42 | 131.1 / 0 | -48.4 (70.0) | 0.28 (3.93) | 0.0685 |

#### Ring 1 K1, ramp

| setting | time s | H RMS % | H nc % | H max % | grad raw / nc % | E err % | E lines % | vert max |
|---|---|---|---|---|---|---|---|---|
| ours Q8/R8/B4 | 0.050 | 1.200 | 1.184 | 6.73 | 42.2 / 36.6 | 13.35 (1.16) | -3.75 (0.16) | 0.0848 |
| ours Q32/R16/B8 | 0.154 | 1.162 | 1.158 | 6.47 | 40.5 / 36.9 | 13.27 (1.15) | -3.54 (0.11) | 0.0791 |
| ours Q2560/R640/B64 | 11.973 | 1.161 | 1.158 | 6.47 | 40.1 / 36.9 | 13.28 (1.15) | -3.59 (0.10) | 0.0790 |
| direct N64 / N128 | 0.912 / 1.799 | 1.726 / 1.213 | ~0 | 8.9 / 5.2 | 358 / 231 raw | noise | -11.4 / -4.7 | 0.053 / 0.032 |
| direct N256/B16 | 3.620 | 0.859 | 0.000 | 4.02 | 165.7 / 44.6 | 0.95 (43.09) | -0.37 (2.18) | 0.0359 |
| direct N1024/B16 | 14.391 | 0.441 | 0.049 | 1.92 | 84.7 / 27.6 | 5.17 (19.67) | 0.90 (0.67) | 0.0133 |
| direct N2560/B64 | 36.131 | 0.288 | 0.039 | 1.44 | 53.1 / 0 | 0.08 (8.45) | 0.83 (0.34) | 0.0091 |

#### Ring 1 K2, ramp

| setting | time s | H RMS % | H nc % | H max % | grad raw / nc % | E err % | E lines % | vert max |
|---|---|---|---|---|---|---|---|---|
| ours Q8/R8/B4 | 0.030 | 1.577 | 1.565 | 7.27 | 62.9 / 58.5 | 36.28 (2.87) | -8.97 (0.37) | 0.1482 |
| ours Q32/R16/B8 | 0.089 | 1.529 | 1.526 | 7.18 | 62.7 / 58.7 | 36.38 (2.86) | -8.69 (0.16) | 0.1469 |
| ours Q2560/R640/B64 | 5.664 | 1.532 | 1.530 | 7.23 | 62.6 / 58.7 | 36.42 (2.86) | -8.60 (0.15) | 0.1479 |
| direct N32 / N64 | 0.531 / 1.052 | 2.209 / 1.572 | 0.1-0.2 | 9.1 / 6.7 | 702 / 479 raw | noise | -0.9 / 1.4 | 0.110 / 0.053 |
| direct N256/B16 | 4.269 | 0.791 | 0.102 | 3.34 | 241.6 / 92.6 | noise | 1.39 (2.03) | 0.0193 |
| direct N1024/B16 | 17.010 | 0.399 | 0.026 | 1.83 | 102.6 / 0 | -5.22 (40.70) | -0.20 (0.81) | 0.0121 |
| direct N2560/B64 | 42.674 | 0.260 | 0.022 | 1.17 | 68.3 / 0 | -8.55 (19.11) | -0.10 (0.54) | 0.0127 |

#### Ring 1 K3, g = 1 (exact +-1 GT)

- Ours: Q8 0.052 s, Q32/R16 0.186 s, Q256 1.288 s, Q2560 12.747 s. H RMS 1.085-1.089%, max 9.9%, vert max 0.226-0.227 (v3 = 0.773 vs 1).
- Direct: exact at every N (0.194 s at N16 to 30.06 s at N2560).

#### Other runs (ring 2 K1 unless noted)

- **Without the lean step** (both methods): ours Q32/R16 0.618 s vs 0.295 s (2.09x); direct N256 5.404 s vs 3.548 s (1.52x). Errors identical.
- **Looser stop distance** 1e-4 sqrt(2) (biased, opt-in):

| setting | time at 1e-6 | time at 1e-4 | change | H RMS at 1e-4 | mean walk length |
|---|---|---|---|---|---|
| ours Q32/R16 | 0.295 s | 0.209 s | -29% | 1.388% | 34.2 -> 19.2 |
| ours Q64/R32 | 0.580 s | 0.368 s | -37% | 1.407% | - |
| direct N256 | 3.548 s | 2.456 s | -31% | 2.515% | 46.2 -> 31.0 |
| direct N1024 | 14.34 s | 9.94 s | -31% | 1.268% | - |

  Ring 2 K2: ours Q32/R16 0.125 s (vs 0.196 s); direct N1024 13.05 s (vs 16.99 s).
- **Independent re-run** (audit, same build): all error numbers matched bit for bit; times within 5-11% (ours Q32/R16 0.280 s vs 0.295 s; direct N256 3.281 s vs 3.548 s).

---

## 5. What it means for a real mesh (DERIVED, not measured)

Nothing in this section was run on a real mesh. It is arithmetic on the measured 8-tet ring times.

Assumptions:
1. A real nodal ring has n_r = 20-24 tets. Each tet lies in 4 rings.
2. Ours: a tet's local solves (Q, R, Gram) are shared by every ring that contains it. A ring of n_r tets then pays for n_r/4 tets' worth of solves, i.e. (n_r/4)/8 = n_r/32 of the measured 8-tet ring time.
3. Direct: a ring solve is not shared. It costs about n_r/8 times the measured 8-tet direct time, at the same walk length per point.
4. Glue: about 2e-5 s per node (measured 13-200 us per ring).
5. Same evaluation-point density as here (about 4,000 points per tet). Both methods scale linearly with points.
6. Result: direct/ours per node = 4 x (direct ring time / ours ring time).
7. The direct is treated as if the ring were convex, like the test octahedron. Real rings are star-shaped but usually not convex, so this direct implementation (it uses a convex-domain solver) could not run on them as is. A non-convex version with no-flow faces would cost more per step.
8. The share of cut tets is taken from the 8-tet test ring (7 of 8). A real 20-24-tet ring probably has fewer cut tets, and uncut tets skip Q walks.

Assumptions 3, 7 and 8 all favour the direct, so the real ratio is probably larger than shown.

Ring 2 K1, g = 1:

| accuracy (H RMS) | ours per ring | ours per node (n_r = 20 / 24) | direct per ring | direct per node (n_r = 20 / 24) | direct / ours per node |
|---|---|---|---|---|---|
| 3% | 0.114 s | 0.071 / 0.085 s | 2.476 s, interpolated | 6.2 / 7.4 s | about 87x |
| 2% | 0.114 s (Q8/R8) | 0.071 / 0.085 s | 5.667 s, interpolated | 14.2 / 17.0 s | about 199x |
| 1.5% | 0.190 s, interpolated | 0.119 / 0.143 s | 10.32 s, interpolated | 25.8 / 31.0 s | about 217x |
| 1.394% (ours Q32/R16) | 0.295 s | 0.184 / 0.221 s | about 12.0 s, interpolated | 30.0 / 36.0 s | about 163x |

Other cases:
- Ring 2 K2 at 2%: ours 0.085 s, direct 9.84 s interpolated: about 463x per node. At ours' floor of 1.81%: about 120-560x.
- Gradient, ring 2 K1 (**estimate**): matching ours' 30.9% raw gradient error would take the direct about N = 8,600 walks, about 120 s per ring: about 210x per ring, about 840x per node.
- Ramp data: the direct is much less noisy, so the gap shrinks. From the log-log interpolated times: about 26-96x per node at 3% to 1.5%. (The bench's own DERIVED lines, which use the first measured setting at or below the threshold, give 38-145x.)
- K3 (crack through the ring): the direct is exact and ours keeps 1.1%. The direct wins. The bench prints a per-node ratio of 7.8x for K3, but it is meaningless.
- Below ours' floor (better than about 1.4% K1, 1.8% K2, 1.2-1.5% ramp), and for worst-point and vertex errors (section 4.3), the direct is the only option of the two.

---

## 6. How to run it

Run from the project folder. Build first with the project's build command:

```
"C:\Program Files\CMake\bin\cmake.exe" --build build --config Release --parallel
```

1. **The full benchmark** (all 6 configurations, default grids; estimate about 25 min, since the main run below with two extra grid points took 1615 s):
   ```
   build\Release\partial_cut.exe --ring-bench
   ```
2. **The exact main run of this report:**
   ```
   build\Release\partial_cut.exe --ring-bench --ours 8/8/4,16/16/8,32/16/8,32/64/8,64/32/16,64/256/16,128/64/16,256/128/16,1024/256/16,2560/640/64 --direct 16,32,64,128,256,512,1024,1024/64,2560 --current 256 --repeat 3 --tag final_main
   ```
3. **A quick check** (ring 2 K1, ours 32/16/8 and 64/32/16, direct 64 and 256, one repeat; estimate a few seconds of solving from the times above):
   ```
   build\Release\partial_cut.exe --ring-bench --quick
   ```
4. **One configuration:** add `--ring 1|2 --config 0|1|2` (config 0 = K1, 1 = K2, 2 = K3).
5. **Options:**
   - `--ours "Q/R[/B],..."`, `--direct "N[/B],..."`, `--current "N[/B],..."`, or `none` to skip a method.
   - `--current 2560/64` runs check 4's old path at its own setting (about 70 s).
   - `--shell 1e-4` uses the looser stop distance (biased, opt-in).
   - `--repeat R` (timing repeats, default 3), `--seeds S` (independent seed sets), `--no-lean`, `--no-fuse`, `--tag NAME`.
6. **Output:** `ring_perf/<tag>_rows.csv` (one row per method, setting, data type and seed set), `ring_perf/<tag>_summary.csv` (equal time and equal accuracy), `ring_perf/<tag>_log.txt`. It never overwrites: a taken name gets `_2`, `_3`. It never writes to `gate_cache/` or `gate_report.md`.
7. **Self-test** (lean walks bit for bit, J_xi, which fields need walks, z tests, determinism; about 10-20 s):
   ```
   build\Release\partial_cut.exe --ring-bench --selftest
   ```
8. **The viewer with the fast method:**
   ```
   build\Release\partial_cut.exe --gate coupling --view --check 4 --ring 2 --config 0 --fast
   ```
   - Without `--fast`, the viewer is unchanged (check 4's old path).
   - The panel has a checkbox "fast method (per-tet control variates + lean walks)", off by default. `--fast` ticks it at start, at Q 32 / R 16.
   - When ticked: the per-tet slider becomes Q's walks per point, a second slider sets R's, and a "fast defaults (Q 32 / R 16)" button appears.
   - The GT loads from `gate_cache/` whenever there are 4096 points per tet and 25,600 reference walks, whatever the per-tet walks are.
9. **The viewer's computation without a window** (prints the panel lines and a hash of every field):
   ```
   build\Release\partial_cut.exe --ring-bench --view --fast --ring 2 --config 0
   ```

---

## 7. Caveats

1. **Ours has a fixed floor that more walks do not lower.** It comes from the glue with Level-C face data:
   - H RMS 1.36% (K1), 1.73-1.78% (K2), 1.10% (K3) after removing noise; worst point 11-15%; ring 2 K1 worst vertex v5 (error 0.12); K3 v3 = 0.773 instead of 1.
   - Energy +1.9-2.2% (K1), +8.9-9.5% (ring 2 K2; check 4's old run gives +9.36% too), +13.3% (ramp K1), +36.4% (ramp K2).
   - Gradient floor after removing noise: 21% (K1), 30% (K2), 37-59% (ramp).
   The direct has none of these floors. With enough walks it beats ours on every metric.
2. **"Ours wins" is about average value error and pointwise gradient noise.** On worst-point error and on vertex values the direct catches up after about 1-4 s (K1), 0.5-2 s (ramp), 8-43 s (K2). See section 4.3.
3. **K3 relative metrics are undefined.** The GT is +-1 with zero gradient, so relative gradient and energy errors divide by zero. The bench prints nan or huge numbers there. That is a reporting artifact, not a solver bug. Use the absolute columns.
4. **The direct's noise-corrected energy and gradient are unreliable at N <= 1024.** Its noise is tens to thousands of %, so the corrected value swings between 0 and 150%. Only the raw gradient trend and the E-lines region are informative there.
5. **Gradient errors include the GT's own 15-26% noise,** and the denominator is not corrected for it. This affects both methods the same way.
6. **One-sided differences:**
   - Only ours gets the control variates and the dropped face tests (its faces all have given values).
   - A J-based control variate for the direct was not built. It would need walks with nonzero no-flow data.
   - The direct always computes both data types (g = 1 and ramp) from the same walks, even in ring 2 where only g = 1 is used. This costs it a little extra (not measured; **estimate**: a few %, no conclusion changes). In ring 1 the opposite holds: one direct run serves both data types, while ours runs twice.
   - The stop distance has the same number but not the same length in space: ours applies it in each tet's xi frame, the direct in physical space.
7. **Best of grid, one seed.** Each equal-time value is the best grid setting from one seed set. Across 3 seed sets values move by about +-0.05 points. The direct's runs at different N share seeds, so they are not independent either. No ranking changes.
8. **Equal-accuracy times for ours are upper bounds.** Its cheapest setting already meets most thresholds. Interpolated times are marked [ ].
9. **Not charged to either method:** building the point sets and weights (0.010-0.014 s). The evaluation uses about 32,000 points per ring, which is the error-measurement density, not necessarily what XFEM needs.
10. **Loose stop distance rows are biased and opt-in.** Not used in any headline number.
11. **Timing noise:** 14-32% background load on a hybrid CPU (4 fast cores, 12 slower ones). Absolute times jitter by a few %, up to about 12% across seed sets. Ratios are stable. The old-path baselines (over 30 s) ran once each.
12. **Section 5 is derived,** not measured on a real mesh.

## 8. What was not tried

1. GPU walks.
2. Reusing walks across re-cuts for the per-tet solves (the cache exists for the single-tet path only).
3. Fewer or better-placed sample points per tet for the glue integrals (the glue needs only 4x4 + 4 numbers per tet). All runs here use check 4's points.
4. Further per-step speed-ups beyond the lean step: a faster random-direction generator, a faster closest-point search on the crack, removing the per-walk value callbacks.
5. A control variate for the direct solve.
6. Lowering ours' floor (better face data than Level C). This is the only way to make ours more accurate; more walks do not help.
7. A real mesh with 20-24-tet, non-convex rings.

Raw data: `ring_perf/final_*` (main sweep, old-path baselines, no-lean, loose, seeds) and `ring_perf/audit5_*` (independent re-run).
