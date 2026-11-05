# XFEM: the nodal ring, cracked and pulled open

A small CPU XFEM solid (neo-Hookean, lumped mass, backward Euler, Newton) that runs
partial_cut's partial-cut enrichment on check 4's nodal ring and shows it in Polyscope.
It is part of the partial_cut solution: the example calls partial_cut's code directly
(`partial_cut_core`, `buildXfemRing` in `src/coupling_gate.h`), so the ring and the glued
harmonic field h are computed in the same process, about a second per crack. The XFEM is
ported from the motus library (not a build dependency) and simplified: one CPU, no buffers
or devices, no VBD, no gravity. Every ported piece names the motus file it comes from.

## Build

With partial_cut itself, from the partial_cut folder (see its README):

```
"C:/Program Files/CMake/bin/cmake.exe" --build build --config Release --parallel
```

This builds the `xfem` library, `xfem_tests` and `xfem_nodal_ring` (folder "XFEM" in
`build/partial_cut.sln`) into `build/Release/`, next to `partial_cut.exe`. The `xfem`
library depends on Eigen only; the example also links `partial_cut_core` and Polyscope.

## Export files (optional)

The example needs no file. The same data can be written to a file, for `--data` and for
the self-tests' real cases. From the partial_cut folder:

```
build/Release/partial_cut --gate coupling --export-xfem --check 4 --ring 2 --config 0 --out XFEM/data/ring2_K1.txt
build/Release/partial_cut --gate coupling --export-xfem --check 4 --ring 2 --config 1 --out XFEM/data/ring2_K2.txt
build/Release/partial_cut --gate coupling --export-xfem --check 4 --ring 2 --config 2 --out XFEM/data/ring2_K3.txt
```

`--config 0 | 1 | 2` is K1 | K2 | K3. `--quad-level L` (default 1) sets how often every
sub-tet is refined before the front refinement. The three runs take 1.5, 0.9 and 0.4 s on
16 threads (256 walks per export point; `--walks-tet 2560 --walks-r 2560` is the old
setting, 10x the walks for a 0.3% change in K1's opening). The format is in `xfem/ring_export.hpp`.

## Run

From the partial_cut folder:

```
build/Release/xfem_nodal_ring                     # the window, K1 on ring 2, solved here by partial_cut
build/Release/xfem_nodal_ring --case K2           # K1 | K2 | K3
build/Release/xfem_nodal_ring --ring 1 --quad-level 0              # check 4's ring 1; the coarsest quadrature
build/Release/xfem_nodal_ring --data XFEM/data    # the export files instead of solving
build/Release/xfem_nodal_ring --no-window --case K1 --frames 480   # stats per frame, no window
build/Release/xfem_nodal_ring --case K1 --frames 140 --screenshot k1.png [--with-ui]
build/Release/xfem_nodal_ring --mode 2 | --mode 3 | --free  --profile ramp|step|sine  --load N
                              --ramp s --period s --dt s --E Pa --nu v --rho kg/m^3 --tol t
                              --front-jump zero|free  --data XFEM/tests/data (the synthetic rings)
build/Release/xfem_tests                          # the self-tests (exit code = failures)
```

Each crack is solved by partial_cut when it is first shown (K1 about 1.5 s, K2 0.9 s, K3 0.4 s
at level 1, printed on the console with partial_cut's self-checks) and kept for the rest of
the run. The live solve and `--data XFEM/data` give the same numbers: same seeds, same
solve.

**The panel.** Pick the crack, the load mode (I, II, III, or free equal-and-opposite
forces), the profile (ramp and hold, step, sine), the load in N and the ramp time or the
period. Play / pause, step, reset, and **crack view**. Polyscope's own Reset View goes to
its home view; crack view goes back to the view along the front. Under "material and
solver", applied on reset: **jump = 0 on the front**, dt, E, nu, rho and the Newton
tolerance. The stats show:

- Newton's iterations and residuals, and the energies;
- the largest crack opening, and the smallest and mean gap along n (negative means the
  faces overlap: there is no contact);
- the jump on the front and the smallest det F;
- a plot of the opening over time.

What is drawn: the deformed ring surface coloured by h (blue 0, white, red 1), the crack's +
sheet in red and its - sheet in blue, the front in yellow, the rest ring's edges, the nodes
(pinned purple, loaded yellow) and the load arrows. Every vertex is carried by the field of
its own tet on its own side, so the crack shows as the sheets parting.

## The method

`xfem/xfem_enrichment.hpp` has it in full.

- **Tet classes** from the export: 0 the crack misses the tet, 1 it cuts the tet in two,
  2 the front is inside it.
- **Node types per node**: harmonic if a class-2 tet is in the node's support
  (psi_i = h - h(X_i), h = (H/g + 1)/2 from partial_cut); else heaviside if a class-1 tet
  is (psi_i = H(phi) - H(phi_i), motus's +1/-1); else not enriched.
- **Why per node, not per tet.** Next to the front, h is not a +1/-1 step even in a tet the
  crack cuts through. Switching the function tet by tet would make a node's enrichment jump
  across the face between a fully cut tet and a front tet, which is a fake crack along that
  face. One function per node can only jump on the crack. A fully cut tet whose corners are
  all heaviside still runs motus's two-piece path, and in K3 that is every cut tet. In K1 and
  K2 every tet holds the node v0, which is harmonic, so all 8 tets take the quadrature path.
- **Three element paths.**
  - plain: one F;
  - two-piece: motus's +3/-1 and 2+/2- pieces, a constant F on each, energy =
    sum of V_piece W(F_piece);
  - quadrature: F = sum x_a gradN_a^T + sum a_a (psi_a gradN_a + N_a grad psi_a)^T at every
    export point, energy = sum w_q W(F_q). The export puts about 1400 points per tet in K1,
    and up to 3200 in a front tet.
- **Mass**: lumped. rho V / 4 per node, unchanged by cutting; rho times the integral of
  N_a psi_a^2 per enrichment.
- **Time stepping**: backward Euler's incremental potential, minimised by Newton with a
  sparse LDLT. The line search starts at the largest step that keeps det F > 0 at every point
  and piece.
- **The jump is held at zero on the front.** On the front h jumps by its full 1, so
  |grad h| ~ 1 / (2 pi r) and the row N_a grad h makes any jump that is nonzero there cost
  an energy that grows like log(1 / r). In the continuum that energy is infinite. Left free,
  the jump is held only by where the nearest quadrature point falls: the answer then depends
  on the quadrature level, and the faces overlap near the front. So each front tet gets
  sum_a N_a a_a = 0 at both ends of its piece of the front. That sum is linear in the tet, so
  the jump is zero along the whole front. Newton runs in the null space of these rows. This
  is the limit that refining the quadrature tends to; `--front-jump free` turns it off.
- **Loads**: external forces only, on the outer vertices' u_i. Mode I pins u_i (motus's
  fix_node) of the outer vertices on the - side and pulls the + side along n. Each step takes
  the load averaged over the step.

## Results

Defaults: E = 2e4 Pa, nu = 0.3, rho = 1000 kg/m^3, dt = 1/60 s, mode I. K1 and K2 use ramp
and hold over 1 s, at 400 N for K1 and 3000 N for K2: K2's crack is small, so the uncracked
part carries most of the pull. K3 uses the sine at 400 N. Over 8 s (480 steps), with the
ring's size 2:

| crack | profile | opening after 2 s (% of 2) | peak | smallest gap | min det F | Newton | step |
|---|---|---|---|---|---|---|---|
| K1 | ramp and hold | 1.9 to 9.1 | 8.8 at 2.3 s | 0 | 0.70 | 1-2 | 110 ms |
| K1 | sine, 3 s | 1.1 to 8.3 | | -0.023 | 0.72 | 1-2 | 112 ms |
| K2 | ramp and hold | 3.3 to 12.4 | 12.8 at 1.3 s | 0 | 0.88 | 1-2 | 55 ms |
| K2 | sine, 4 s | 0.5 to 8.8 | | -0.006 | 0.88 | 1-2 | 57 ms |
| K3 | sine, 3 s | 0.9 to 13.5 | | -0.075 | 0.99 | 2 | < 1 ms |
| K3 | ramp and hold | 24 at 2 s, then it flies off | | 0 | 0.995 | 1-2 | < 1 ms |

Under ramp and hold the crack opens from the front as a wedge and keeps swinging, since
nothing damps it but backward Euler. The sine pulses from 0 to the full load where there is
a front, so the crack opens and shuts every period. It swings a little past shut (mean gap
down to -0.0014 for K1, -0.0020 for K2) because there is no contact. With no front, the sine
is a cosine about 0, which keeps K3's free + half near.

Other loads on K1:

| run | opening | min det F |
|---|---|---|
| mode II | 24.9% (mostly slip) | 0.37 |
| mode III | 4.7% | 0.96 |
| free | 13.6% | 0.62 |
| step | 9.5% | 0.68 |

Mode III twists the + half, and part of the faces pass through each other (mean gap down
to -0.027).

**Held or free jump on the front.** The demo's K1 run with the jump free:

- Newton takes up to 4 iterations;
- min det F drops to 0.048;
- the faces overlap by up to 0.011 at the front;
- the jump on the front is 0.019 (held: 2e-16).

The quasi-static mode I solve (10 load steps) at three quadrature levels, mean gap along n
beyond 0.4 from the front:

| | level 0 | level 1 (default) | level 2 |
|---|---|---|---|
| K1 points | 2212 | 11312 | 63112 |
| K1 held | 0.0816 | 0.0811 | 0.0811 |
| K1 free | 0.0947 | 0.0919 | 0.0907 |
| K2 held | 0.1462 | 0.1441 | 0.1436 |
| K2 free | 0.1376 | 0.1348 | 0.1331 |

Held, the levels agree, and no sheet pair overlaps. Free, the answer keeps moving, about 24
of K1's 2615 sheet pairs overlap (by up to 0.0056), K2's faces overlap by 0.047 to 0.063,
and K1's min det F falls from 0.21 to 0.06.

**Tests**: `xfem_tests` passes 52 of 52. It covers:

- the export format;
- plain FEM with no enrichment;
- finite differences of every path and of the incremental potential;
- the exact Heaviside on the quadrature path, which equals the two-piece path;
- affine and rigid states, the masses, and a free run that only loses energy;
- (8) no seams along tet faces: at most 7.4e-14 of the jump;
- (9) mode I: K1's mean opening is 0.006 within 0.05 of the front and 0.081 beyond 0.4;
- (10) the front's constraint holds along the whole front to 2e-15 of the opening;
- (11) quadrature-level convergence on the synthetic ring: 0.14% held, 0.94% free.

## Limits

- **Flat cracks with a straight front only.** phi = n.X - d, the front is one segment, and
  the front constraint clips that segment to the front tets.
- **The 8-tet ring only.** Assembly is sparse and the constraint basis is local, so a
  bigger mesh should work, but only this ring has been run.
- **Monte Carlo noise in h and grad h.** Both are walk estimates at the quadrature points.
  For K1, the s.e. of h is 1.3e-4 (median) and 3.6e-4 (max). The s.e. of |grad h| is 0.016
  (median) and 0.27 (99th percentile), against a median |grad h| of 1.2. So h is continuous
  across tet faces only to that level. h on the drawn surfaces is exact face data.
- **The crack mouth.** Where the crack meets the outer surface, the surface carries the face
  data and the sheets exactly 1 / 0. The two edges part by 0.75% (K1) and 0.58% (K2) of
  the jump; in K3 they meet.
- **Linear opening at the front.** The opening grows linearly from the front inside a tet,
  not like sqrt(r), since that is what linear shape functions allow.
- **No contact and no damping.** Faces pushed together pass through each other. The panel
  says so when the mean gap is negative.
- **K3's free half turns slowly under the sine**, because the loads keep their direction
  while it rocks. One edge dips into the - half by 0.075 within 8 s.
- **Cost.** Level 1 costs about 110 ms per step for K1. Level 0 gives the held answers
  within 0.6% (K1) and 1.8% (K2) of level 2, with a fifth of level 1's points.

## Layout

| path | what |
|---|---|
| `xfem/types.hpp`, `tet_mesh.hpp`, `tet4.hpp`, `neo_hookean.hpp` | scalars, the mesh, the linear tet, the material |
| `xfem/xfem_cutter.*` | the +3/-1 and 2+/2- pieces, piece volumes, shape integrals |
| `xfem/xfem_enrichment.*` | tet classes, the per-node rule, psi, the element paths |
| `xfem/xfem_system.*` | the element paths, the mass, energy / gradient / hessian, the inversion filter, the front constraint |
| `xfem/backward_euler.*` | the incremental potential, Newton (in the constraints' null space), the line search |
| `xfem/ring_export.*`, `synthetic_ring.*` | the export reader, writer and checks; a stand-in export for the tests |
| `examples/nodal_ring/` | the Polyscope example; `live_ring.*` takes the ring from partial_cut_core (buildXfemRing) |
| `tests/` | the self-tests; `tests/data/` the synthetic rings (`xfem_tests --write-synthetic`) |
| `data/` | partial_cut's exports, `ring2_K1.txt`, `ring2_K2.txt`, `ring2_K3.txt` (for `--data` and the self-tests) |
