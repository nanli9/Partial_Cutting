# Ring stress test: our method vs the ring GT on many crack types

Date: 2026-11-04. Exe: build\Release\partial_cut.exe (built 2026-11-03 22:39). Total compute: about 3.5 h.

## The short answer

- Our method holds up on flat cracks, complete cuts (flat and curved), a gently curved front, curved crack surfaces and a badly distorted ring. In these cases the error is 1-2%, the same as the flat baseline (1.4%).
- It degrades (2-5%) on a strongly curved front, a kinked crack, two parallel cracks, a crack running close to an interior face, and a front close to the outer boundary.
- It fails on two crossing cracks (borderline), two small cracks in one tet, a tiny crack in one tet, and a crack 0.001 from an interior face.
- In the cases that degrade or fail, more walks do not help. The leftover error comes from the allowed face values (J plus a linear function per face), not from noise.

## Terms

- **GT**: the reference answer. One walk-on-stars solve on the whole ring, 25,600 walks per point.
- **Direct**: the same solver as the GT, run with fewer walks so it takes about as long as ours.
- **avg / worst point**: error over the whole ring, as a % of the jump 2g (the RMS, then the largest single point).
- **worst vertex**: the largest error of the 7 glued vertex values (H runs from -1 to +1).
- **floor**: our error with about 30x more walks (Q1024/R256/B16). It is the error that is left once the noise is gone.
- **default**: our fast setting Q32/R16/B8 (median of 3 seed sets; 6 for D2, E3, F3).

Ring 2 (node off-centre), g = 1, 2048 points per tet unless noted.

## Main table

| Case | Ours avg / worst point / worst vertex | Our floor | Our time | GT time | Direct at our time | Direct time to match our floor | Verdict |
|---|---|---|---|---|---|---|---|
| A1 flat crack K1 (baseline) | 1.39% / 14% / 0.12 | 1.38% | 0.40 s | 570 s | 8.1% | > 7.9 s | works |
| A1s same, at stress-test density | 1.41% / 13% / 0.11 | 1.39% | 0.23 s | 204 s | 7.5% | > 4.0 s | works |
| A2 flat crack K2 | 1.82% / 11% / 0.07 | 1.80% | 0.22 s | 686 s | ~14% | > 9.0 s | works |
| A3 flat complete cut | 1.11% / 11% / 0.23 | 1.11% | 0.40 s | exact (+-1) | exact | - | works |
| B1 curved complete cut (cylinder) | 0.95% / 4.2% / 0.04 | 0.91% | 0.77 s | exact (+-1) | exact | - | works |
| B2 flat complete cut 0.001 from the node | 1.10% / 9.1% / 0.23 | 1.10% | 0.20 s | exact (+-1) | exact | - | works |
| C1 flat crack, gently curved front | 1.28% / 7.7% / 0.09 | 1.24% | 0.88 s | 735 s | 7.6% | > 15 s | works |
| C2 flat crack, strongly curved front | 3.49% / 17% / 0.27 | 3.47% | 0.83 s | 1025 s | 9.7% | 6.4 s | degrades |
| D1 curved crack (sphere, radius 2) | 1.52% / 15% / 0.24 | 1.26% | 1.3 s | 849 s | 6.4% | > 17 s | works |
| D2 strongly curved bowl (radius 0.6) | 1.6% / 8.7% / 0.10 | 1.59% | 2.1 s | 1289 s | 6.3% | > 26 s | works (noisy) |
| D4 sine-wave crack (wavelength 1, tightest bend radius 0.25) | 1.83% / 19% / 0.20 | 1.83% | 0.94 s | 836 s | 7.5% | ~16 s | works |
| D5 sine-wave crack (wavelength 0.6, tightest bend radius 0.09) | 2.06% / 18% / 0.11 | 1.99% | 1.3 s | 1202 s | 8.1% | ~19 s | degrades (just over 2%) |
| D3 kinked crack (40 degrees) | 2.20% / 18% / 0.06 | 2.20% | 0.36 s | 359 s | 6.1% | 2.8 s | degrades |
| E1 two parallel cracks (g = 1 and 0.5) | 3.81% / 22% / 0.29 | 3.51% | 0.29 s | 269 s | 5.4% | 0.67 s | degrades |
| E2 two small cracks in one tet | 7.0% / 15% / 0.24 | 5.81% | 0.20 s | 780 s * | ~19% | 2.2 s | fails |
| E3 two crossing cracks (X junction) | 5.4% / 30% / 0.31 | 4.71% | 0.25 s | 248 s | 6.2% | 0.43 s | fails (borderline) |
| F1 front through the node | 1.61% / 14% / 0.07 | 1.64% | 0.09 s | 191 s | ~13% | > 3.7 s | works |
| F2 crack 0.001 from an interior face | 13.2% / 21% / 0.39 | 12.9% | 0.09 s | 216 s | ~16% | < 0.14 s | fails |
| F2b same crack 0.03 from the face | 4.27% / 12% / 0.18 | 4.16% | 0.10 s | 223 s | ~15% | 1.2 s | degrades |
| F3 tiny disk inside one tet | 4.7% / 6.7% / 0.11 | 4.14% | 0.13 s | 1154 s * | ~30% | 6.7 s | fails ** |
| F4 front 0.03 from the outer boundary | 3.61% / 21% / 0.33 | 3.61% | 0.21 s | 182 s | 1.4% | < 0.12 s | degrades |
| G1a K1 on a distorted ring | 1.83% / 14% / 0.12 | 1.79% | 0.22 s | 199 s | 7.6% | > 3.9 s | works |
| G1b sphere crack on a distorted ring | 1.72% / 15% / 0.25 | 1.84% | 1.3 s | 837 s | 6.4% | 16 s | works |

Notes:
- E1's errors are a % of 2 x the larger g (= 2).
- "~" = extrapolated: ours was faster than the cheapest direct run (16 walks per point).
- "> x s" = direct had still not reached our floor at its largest run (512 walks per point). "< x s" = direct's cheapest run already beat our floor.
- A1, A2, A3 use check 4's 4096 points per tet and its cached GTs (their GT times are check 4's original runs). * E2 and F3 use 1024 points per tet.
- ** F3: the automatic rule says "degrades" (avg <= 5%). I call it a failure because our error is about 3 times the size of the field itself. The six seeds range from 1.5% to 8.2%.
- D2: the six seeds give 1.5-2.6%. One seed in six had a heavy tail (2.6%).
- E3: the default (5.4%) is just over the 5% line, but the floor (4.7%) is just under it.
- "Our time" covers all face-function variants of a case, solved together. Geometry and point setup (at most 0.2 s) is not charged to either method.

## What works, what breaks, and why

- **Flat, complete and gently curved cuts work.** The new general crack code reproduces the existing fast path (A1: 1.39% vs 1.394%). With the plain side indicator as J, the complete cuts (A3, B1, B2) come out exactly +-1, as they must. Slivers at the node, a front through the node and a distorted ring (+30% error) cause no trouble.
- **In every degrading case, the floor equals the default.** So more walks do not help. For C2 and D3 a plain solve without control variates gives the same error. The cause is the allowed face values (J + linear per face), not noise.
  - C2: the error sits near the outer boundary, where the crack meets the outer faces along a curve. The theta-type J does better here (2.60%) than the solid-angle J (3.49%).
  - F4: the narrow uncut strip (the ligament) is the problem (9.6% near the front).
  - E1, E3: the blended J for several cracks is the weak point. E3 is worst at the junction (7.0%).
- **A crack close to and parallel to an interior face** (F2, F2b) shifts all vertex values, so the error is even over the whole ring. Shrinking the left-out tube around the crack lowers F2 from 13% to 7.5%. That is still a failure.
- **Only one tet cut** (E2, F3): the 7 vertex values are pinned only weakly, so they come out offset, with a large seed spread. F3 is very sensitive to how the energy near the small crack is sampled: a smaller tube moved its floor from 4.1% to 14.1%, so a smaller tube is not a fix. Likely cause (not verified): the vertex values are pinned only through that near-crack energy.
- **Error compared with the field's own size.** % of 2g flatters cases with small fields. Our error is about 300% of the field in F3, 111% in F2, 60% in E2 and 32% in F2b. It is 15% or less in all other cases.
- **Ours vs direct.** At equal time, ours wins everywhere except F4 (direct 1.4%, because almost every walk ends quickly on the nearby crack) and F2 (direct's cheapest real run, 12.2% at 0.14 s, already beats our 13.2%). But where our error stops at a floor, a little more direct time passes it: E3 0.4 s, E1 0.7 s, F2b 1.2 s, E2 2.2 s, D3 2.8 s, C2 6.4 s, F3 6.7 s. In the cases that work, direct needs more than 4-26 s.
- **Curved cracks are noisier and slower at the default.** The curved control variate leaves more noise. More R walks (Q32/R64) bring D1 to 1.21% and D2 to 1.62% at about 2.5x the time. Curved cases take 0.8-2.1 s, against 0.1-0.4 s for flat ones.
- **Gradients and energy.** In the cases that work, the gradient error is 17-36% (noise removed) and the energy error is 1-13%. With the ramp g, the value error is lower (A1 1.22%, A2 1.56%, E1 3.28%) but the energy error is higher (+15%, +40%, +84%).

## How sure the GT is

- GT noise is 0.05-0.31% of 2g in every case, with no walks capped or escaped.
- Checked independently:
  - flat: A1 and A2 match check 4's cached GTs;
  - exact: A3, B1 and B2 give exactly +-1;
  - Zombie: an independent solver agrees for A1s, D2, D5, E1 and E3.
- No independent check: C1, C2, D1, D3, D4, E2, F1-F4, G1a and G1b. They rely on the same new GT code that passed the checks above. D1's re-run with new seeds checks only its noise.
- Near the curved and multiple cracks (D2, E1, E3), Zombie and the GT scatter apart more than their noise allows. The excess is real, but it is below 0.2% of 2g (0.5% at the E3 junction). That is 20-30x smaller than our errors there, so no verdict changes. The cause was not found.

## What to do next

- Build a better face function for multiple cracks and junctions. The blend is the weak point (E1, E3).
- Allow richer face values than linear near the outer boundary and the ligament (C2, F4).
- Find a way to pin the vertex values when only one tet is cut (E2, F3).
- Use Q32/R64 as the default for curved cracks, and the theta-type J for strongly curved fronts.
- Consider a fallback to the direct solve where our floor is known to be high: a crack near an interior face, a tiny crack in one tet, a front near the boundary.

## How to rerun

From the project folder (the cached GTs are reused, so only our runs, the direct runs and the Zombie checks are redone; new outputs get a _2 suffix):

```
build\Release\partial_cut.exe --ring-stress --selftest
build\Release\partial_cut.exe --ring-stress --case A1,A1s,A2,A3 --tag stress_A
build\Release\partial_cut.exe --ring-stress --case B1,B2,C1,C2 --tag stress_BC
build\Release\partial_cut.exe --ring-stress --case D1,D2,D3 --zombie D2,E3 --tag stress_D
build\Release\partial_cut.exe --ring-stress --case E1,E2,E3 --zombie D2,E3,E1 --tag stress_E
build\Release\partial_cut.exe --ring-stress --case F1,F2,F3,F4,F2b --tag stress_F
build\Release\partial_cut.exe --ring-stress --case G1a,G1b --tag stress_G
build\Release\partial_cut.exe --ring-stress --case D4,D5 --zombie D5 --tag stress_D45
build\Release\partial_cut.exe --ring-stress --case E1,D2,F3,E3 --ours 32/16/8 --direct 16 --seeds 6 --repeat 1 --tag audit_seeds
build\Release\partial_cut.exe --ring-stress --case F2,F3 --eps 0.0005 --tag diag_small_eps
```

To look at one case in the viewer: `build\Release\partial_cut.exe --gate coupling --view --check 4 --case D2`.

Raw data, in ring_stress\:
- stress_all_table.csv: all cases in one table.
- stress_{A,BC,D,E,F,G,D45}_{rows,summary}.csv, plus _log.txt and _geometry.txt: per group (D45: the sine-wave cases, added later).
- audit_seeds_*: the extra seeds.
- audit_rerun_*: the reproduction run.
- diag_small_eps_* and diag_zombie_A1s_*: the diagnostics.
- Ignore the dev_*, smoke_*, nonreg_* and review_* files there; they come from earlier runs at reduced settings.

New GT caches are in gate_cache\ref_stress_*.bin (20 files). The existing caches and reports were not touched.
