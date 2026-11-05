# XFEM: the nodal ring, cracked and pulled open

A small CPU XFEM solid (neo-Hookean, lumped mass, backward Euler, Newton) that runs
partial_cut's partial-cut enrichment on check 4's nodal ring and shows it in Polyscope.
It is part of the partial_cut solution: the example calls partial_cut's code directly
(`partial_cut_core`, `buildXfemRing` in `src/coupling_gate.h`), so the ring and the glued
harmonic field h are computed in the same process, about a second per crack. The XFEM is
ported from the motus library (not a build dependency) and simplified: one CPU, no buffers
or devices, no VBD. Every ported piece names the motus file it comes from.

A second example, `xfem_beam`, runs the same method on a real mesh: the cantilever of the
`cutting` project's GXFEM beam example, cracked by partial, complete, curved and several
cracks ([the beam](#the-beam-a-real-mesh) below).

## Build

With partial_cut itself, from the partial_cut folder (see its README):

```
"C:/Program Files/CMake/bin/cmake.exe" --build build --config Release --parallel
```

This builds the `xfem` library, `xfem_tests`, `xfem_nodal_ring` and `xfem_beam` (folder
"XFEM" in `build/partial_cut.sln`) into `build/Release/`, next to `partial_cut.exe`. The
`xfem` library depends on Eigen only; the examples also link `partial_cut_core` and Polyscope.

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

**Tests**: `xfem_tests` passes 55 of 55. It covers:

- the export format;
- plain FEM with no enrichment;
- finite differences of every path and of the incremental potential;
- the exact Heaviside on the quadrature path, which equals the two-piece path;
- affine and rigid states, the masses, and a free run that only loses energy;
- (8) no seams along tet faces: at most 7.4e-14 of the jump;
- (9) mode I: K1's mean opening is 0.006 within 0.05 of the front and 0.081 beyond 0.4;
- (10) the front's constraint holds along the whole front to 2e-15 of the opening;
- (11) quadrature-level convergence on the synthetic ring: 0.14% held, 0.94% free;
- (12) the beam's mesh loads as the beam example meshes it: 85 nodes, 203 tets, volume 0.04;
- (13) the uncut beam is plain FEM, and its gravity is the lumped mass times g;
- (14) a cut through the beam: the per-node phi gives the plane's enrichment, the + part
  moves rigidly at no energy, and gravity does rho V_+ g.t of work on it (the signs of psi).

## The beam: a real mesh

`xfem_beam` takes the cantilever of the `cutting` project's GXFEM beam example: its mesh,
reproduced exactly (`data/beam/`, its README shows every start-up statistic matching), its
wall (every node at x = -0.5 pinned, u_i and a_i), its gravity (9.81 m/s^2 along -y),
material (E = 1e6 Pa, nu = 0.4, rho = 1000 kg/m^3) and dt = 1/60 s. The beam sags uncut
until `--cut-frame` (40, the example's); then the case's cracks go in. Their enrichment is
computed by partial_cut in the same process, once, at start-up: `buildXfemMesh` in
`src/xfem_mesh.h` (`examples/beam/live_mesh.*` hands it over). The cracked system starts from
the nodes' positions and velocities with every new jump at 0, so the field and its energy do
not change at the cut. Uncut, the tip drop is the example's to 1% (0.3588 against 0.3593 at
frame 19; the mass is lumped here, consistent in the example's newton).

```
build/Release/xfem_beam                                        # the window, the notch
build/Release/xfem_beam --case multi                           # notch | through | curved_front | curved_surface | multi | mixed
build/Release/xfem_beam --no-window --case through --frames 60 --cut-frame 20
build/Release/xfem_beam --no-window --depth 0.9                # the notch 0.9 of the way down (the example's --depth)
build/Release/xfem_beam --case notch --frames 120 --screenshot notch.png
build/Release/xfem_beam --case mixed --record out/mixed        # PNG frames, and an mp4 if ffmpeg is on the PATH
build/Release/xfem_beam --csv notch.csv --quad-level 1         # one row per frame; the export quadrature once finer
build/Release/xfem_beam --pull 100 --pull-direction 0 -1 0     # no gravity: 100 N pull the free end's bottom edge down
build/Release/xfem_beam --home-view --record out/notch         # the example's own recording camera, to set beside its frames
```

`--pull F` is the example's: no gravity, and F newtons in all pull the free end along
`--pull-direction` (+x unless given). The force is shared equally over the end's nodes that
lead the way (all nine along +x, the bottom edge's three when pulled down). It grows from
nothing over the first second and is held after.

What is drawn is styled after the example: the body and the crack sheets in its beam colour,
every tet's edges in black, its wall, the camera of its clips and 2x supersampling. Once the
cut is in, a red ball marks every enriched node, and a crack's line across a face is drawn
only where the face is cracked, not past the front. With `--pull`, a green arrow on each
pulled node points the way.

The cases, in rest space (the beam is x in [-0.5, 0.5], y and z in [-0.1, 0.1]). Every
blade comes from the top, y = +0.1, where the sag opens it: there is no contact.

| case | cracks |
|---|---|
| notch | x = -0.245 from the top to 0.7 of the height: the example's default cut, partial |
| through | x = 0.2 all the way down (the example's `--depth 1 --cut-x 0.2`): the free end falls away |
| curved_front | x = -0.245, 0.7 deep at z = 0 and 0.4 at the sides (y_front = -0.04 + 6 z^2) |
| curved_surface | a cylinder, axis along z, radius 0.15, through x = -0.245 at the top and bending to +x, 0.6 deep |
| multi | notches at x = -0.245 (0.6 deep) and x = +0.15 (0.5 deep) |
| mixed | a notch at x = -0.245 (0.6 deep) and a cut through at x = +0.245 |

TetGen put nodes on x = -0.25 and x = +0.25, so the example's -0.25 and the mixed cut's
+0.25 are moved 0.005 to the side: every crack's smallest |phi| at the nodes of the tets it
cuts is then 0.005 (4.0% of the mean edge, 0.124) or more (x = 0.2: 10.4%, x = 0.15: 20.1%).

**The method on a mesh** (`src/xfem_mesh.h`). A crack is two level sets at the vertices:
phi, the signed distance to its surface, and psi, negative where the surface is cracked.
Both are linear in each tet, so in every tet the crack is a flat polygon with a straight
front, the per-tet cut check 4 already solves, and a curved crack is piecewise flat and
continuous across faces. The tet classes and the per-node rule are the export's. h is
solved only in each crack's region, the tets with a harmonic corner; a crack cut through
the mesh everywhere has none and needs no walks. In every region tet, in its own xi frame
with all faces Dirichlet: H0 = J_xi + R and P_k = lambda_k + Q_k (the fast path), with the
face function J = atan2(phi_h, psi_h)/pi of the interpolated level sets (continuous across
faces, +-1 on the crack) and Level C's half-plane of the tet's own plane and front as the
control variate J_xi. The glue minimises the region's energy over its vertex values, one
sparse LDLT (check 4's 7 x 7 system), with the region's outer vertices free as O's are. The
export pass then solves each region tet's quadrature (split by its own phi, refined around
the front polyline; bands scaled by the tet's longest edge) at 256 walks per point. The XFEM
takes the per-node phi (`build_enrichment`'s second form), holds the jump at zero at both
ends of every front tet's piece of the front, and loads gravity as rho g int N_a psi_a on an
enrichment (`xfem_system::body_force`).

**Several cracks.** The XFEM has one enrichment per node, and each tet's quadrature carries
one h. So partial_cut refuses (with the cracks' names) two cracks in one tet, a node
enriched for two cracks, and a tet in the supports of nodes enriched for two cracks. multi's
and mixed's cracks are 0.4 and 0.49 apart, and none of the three happens.

**Cost.** The enrichment takes about 1.3 s per partial crack (glue 0.9 s, about 230 000
points in 56 to 64 region tets, 7.6 M walks; export 0.3 s, about 3400 quadrature points at
level 0) and nothing for a cut through. A step takes 5 to 10 ms uncut and, at level 0, about
100 ms with one partial crack and 250 ms with two.

## Limits

- **Flat in each tet.** A crack is flat with a straight front in every tet (phi, psi
  linear), so a curved crack is piecewise flat. The ring takes one plane and one front
  segment (`zero_jump_on_front`); the beam holds the jump tet by tet.
- **Cracks kept apart.** One enrichment per node: see "Several cracks" above.
- **Meshes.** The 8-tet ring and the 203-tet beam have been run.
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
  says so when the mean gap is negative (the beam: the mouth's gap along the crack normal).
- **The beam's mass is lumped**, the beam example's newton takes the consistent one: the
  same static sag, a slightly different swing.
- **K3's free half turns slowly under the sine**, because the loads keep their direction
  while it rocks. One edge dips into the - half by 0.075 within 8 s.
- **Cost.** Level 1 costs about 110 ms per step for K1. Level 0 gives the held answers
  within 0.6% (K1) and 1.8% (K2) of level 2, with a fifth of level 1's points.

## Layout

| path | what |
|---|---|
| `xfem/types.hpp`, `tet_mesh.hpp`, `tet4.hpp`, `neo_hookean.hpp` | scalars, the mesh, the linear tet, the material |
| `xfem/tetgen_mesh.*` | a tet mesh from TetGen's .node / .ele (the beam's) |
| `xfem/xfem_cutter.*` | the +3/-1 and 2+/2- pieces, piece volumes, shape integrals |
| `xfem/xfem_enrichment.*` | tet classes, the per-node rule, psi, the element paths (from a plane or from phi at the vertices) |
| `xfem/xfem_system.*` | the element paths, the mass, the body force, energy / gradient / hessian, the inversion filter, the front constraint |
| `xfem/backward_euler.*` | the incremental potential, Newton (in the constraints' null space), the line search |
| `xfem/ring_export.*`, `synthetic_ring.*` | the export reader, writer and checks; a stand-in export for the tests |
| `examples/nodal_ring/` | the Polyscope example; `live_ring.*` takes the ring from partial_cut_core (buildXfemRing) |
| `examples/beam/` | the beam example; `live_mesh.*` takes the mesh's enrichment from partial_cut_core (buildXfemMesh) |
| `tests/` | the self-tests; `tests/data/` the synthetic rings (`xfem_tests --write-synthetic`) |
| `data/` | partial_cut's exports, `ring2_K1.txt`, `ring2_K2.txt`, `ring2_K3.txt` (for `--data` and the self-tests); `beam/` the beam's mesh and how it was made |
