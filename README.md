# partial_cut: harmonic partial-cut enrichment on one tetrahedron

This is a small gate experiment. It asks whether a Monte Carlo solver gives a usable
harmonic enrichment `H(X)` for a partial cut inside a single tetrahedron, and whether `H`
tends to the ±1 Heaviside as the cut completes. It computes `H` only: there is no
simulation, no caching and only one cut.

- **Domain**: one tetrahedron, by default with vertices (0,0,0), (1,0,0), (0,1,0) and (0,0,1).
- **Cut**: one planar, zero-thickness, double-sided cut with a straight front. The front
  position is a parameter `t` in [0,1]: `t = 0` means no cut and `t = 1` means the plane
  crosses the whole tet.
- **PDE**: `ΔH = 0` in the tet, with zero Neumann (reflecting) conditions on all four faces
  and Dirichlet `H = ±g` on the cut (absorbing): `+g` on the side `n·X > d` and `-g` on
  the other.
- **Solver**: Walk on Stars from the [Zombie](https://github.com/rohan-sawhney/zombie)
  library, which uses [FCPW](https://github.com/rohan-sawhney/fcpw) for geometric queries.
  This is Zombie itself, not a reimplementation of WoSt.
- **Viewer**: [Polyscope](https://polyscope.run).

**Result in one paragraph.** With g = 1, Zombie's WoSt gives a usable H. The error falls
like 1/√N. The jump across the cut is 2g to within 1.1% at ε = 0.01·L, and the
shortfall is the expected O(ε) offset. At t = 1, H is exactly the ±1 Heaviside. A
separate reflected-Brownian-motion simulator agrees with Zombie at all 8 probe values to
within 1 standard error. With g_b = min(1, sqrt(r/L)), H has the right jump but can never
reach ±1 in this tet, because of how g_b is defined. A default solve takes 1.05 s. See
[Validation results](#validation-results).

## Build

Requirements:

- CMake 3.22 or newer.
- A C++20 compiler: GCC 10 or newer, Clang 12 or newer / AppleClang, or MSVC 2022.
- git, and internet access for the first configure.

All dependencies are fetched by CMake `FetchContent` at pinned commits (see
`CMakeLists.txt`):

- Eigen 5.0.1
- fcpw `e498ca8`
- zombie `b232b1c`, with its pcg32, nanoflann and oneTBB submodules
- Polyscope v2.6.1, with its own submodules

The first configure clones everything and takes about 2 minutes (134 s on the test
machine). The first Release build takes about 40 s.

**Linux** (not tested yet). GLFW needs the X11 and OpenGL headers, e.g. on Debian/Ubuntu:

    sudo apt install build-essential cmake git xorg-dev libgl1-mesa-dev
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build --parallel
    ./build/partial_cut

**macOS** (not tested yet). Install the Xcode command line tools, then:

    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build --parallel
    ./build/partial_cut

**Windows** (tested with Visual Studio 2022, MSVC 19.44, CMake 4.4):

    cmake -S . -B build -G "Visual Studio 17 2022" -A x64
    cmake --build build --config Release --parallel
    build\Release\partial_cut.exe

On Windows `tbb12.dll` is copied next to the executable after each build.

From the Visual Studio IDE, open `build\partial_cut.sln`, set the configuration dropdown
to **Release** and press Run. `partial_cut` is the startup project and runs from the
project folder. A Debug build works but is about 60x slower (61 s per solve), and the
window stays blank until the first solve finishes.

## Run

| command | what it does |
|---|---|
| `partial_cut` | Opens the Polyscope GUI. |
| `partial_cut --validate` | Runs the four validation checks headless, prints them and writes `results.txt` into the current working directory. |
| `partial_cut --bench` | Runs one full default solve (lattice with gradient, plus slice) and prints wall time and walk statistics. |
| `partial_cut --screenshot out.png [--t T] [--with-ui]` | Builds the scene with the front at `T` (default 0.6), solves, saves a 1280×960 Polyscope screenshot and exits. `--with-ui` also draws the control panels. |
| `partial_cut --help` | Prints the usage. |

A full default solve takes 1.05 s (median of 3 `--bench` runs: 1.04, 1.05 and 1.07 s) on
the 16-thread Windows test machine. It covers 433 lattice points, with ∇H at 156 of them,
and 1481 slice points, at N = 256. The slice takes about 0.8 s of that. `--validate` takes
4.6 s.

The GUI writes `imgui.ini` and `.polyscope.ini` into the working directory. Both are
gitignored. `--screenshot` writes neither.

## GUI controls

Every change re-solves the whole scene from scratch with a fresh seed; nothing is reused.
A slider triggers the solve when you release it, not while you drag it.

| control | effect |
|---|---|
| **t (front position)** | Front position in [0,1]. The cut is the part of the plane∩tet polygon with `X·u <= s(t)`, where `s(t) = s_min + t (s_max - s_min)` over the polygon's extent along the front direction `u`. At `t = 0` there is no cut. |
| **N walks per point** | Walks per query point, 16 to 4096 (logarithmic slider). |
| **epsilon / L** | The user epsilon, in units of the longest tet edge `L`. It is the lattice exclusion radius (lattice points closer than epsilon to the cut are dropped) and the half separation of the jump-check pairs. The range is 1e-3 to 0.1 (logarithmic), so epsilon always stays at least 10x the walk termination shell and 7x the distance within which walks stop (see below). |
| **(a) g = 1 / (b) g = min(1, sqrt(r/L))** | Dirichlet data type. `r` is the distance to the front line: the infinite line in the cut plane at `X·u = s(t)`, with direction `n×u`. |
| **k (lattice subdivisions)** | Barycentric lattice resolution, 1 to 30. There are `(k+1)(k+2)(k+3)/6` points before exclusion. |
| **offset from centroid** | Moves the slice plane along its normal. |
| **yaw (deg, about n)** | Rotates the slice normal about the cut normal `n`, towards `u`. |
| **pitch (deg, towards n)** | Tilts the slice normal towards `n`. At yaw = pitch = 0 the slice normal is `n×u` and the plane passes through the centroid, so the cut appears on the slice as a segment that ends at the front. |
| **Re-solve (new seed)** | Same settings, new independent sample. |
| **Run validation** | Runs the same four checks as `--validate`, with the current epsilon and k. The current N is used for checks 2 and 4; check 1 has its own N list and check 3 always uses N = 4096. Checks 1 and 3 always run at t = 0.6, both g types are always run, and the base seed is the same as for `--validate`. It runs on the GUI thread, so the window freezes for about 5 s. It writes `results.txt` into the working directory. |

The status lines show:

- epsilon in absolute units, and the distance from the cut within which walks stop;
- the cut polygon (vertex count, area, front position `s`) and the plane∩tet quad check;
- lattice point count (and how many were excluded) and slice point count;
- the time of the last solve, split into lattice and slice;
- walk counts (counted, escaped and discarded, capped) and the mean walk length;
- how many lattice points have a ∇H estimate;
- after **Run validation**, where `results.txt` was written.

At `t = 0` the panel says that no walks are run and that H is set to 0.

What Polyscope shows:

- **tet**: a transparent volume mesh.
- **cut**: a surface mesh of the cut polygon.
- **lattice**: a point cloud of all lattice points with these quantities:
  - `H`, with the coolwarm diverging map fixed at [-1,1];
  - `H std. error`.
- **lattice grad H**: a second point cloud that holds only the interior lattice points,
  the ones where ∇H is estimated (156 of 433 at the defaults; none for k <= 3). It has:
  - `grad H`, as vectors;
  - `|std. error of grad H|`.

  It is hidden at startup; tick it in Polyscope's structure list to show the arrows.
  Points on the tet faces have no estimate, so they are left out instead of showing a
  zero arrow. The cloud's points are drawn smaller than the lattice points, so `H` stays
  visible. To see the `|std. error of grad H|` colours, hide `lattice`.

The control panel and the console also show H ± 1 s.e. at the four tet vertices (the
lattice points with barycentric index k), which a shifted enrichment N_i (H − H(x_i))
needs.
- **slice**: a surface mesh with `H` (coolwarm, [-1,1]).

The arrows are scaled for display only: the arrow length maps to the 90th percentile of
|∇H|, never below 1/L, so the few very noisy estimates next to the cut do not shrink every
other arrow.

## How it works

**Geometry.** Everything is set in one place: `GeometryConfig` in `src/cut_problem.h`.
It holds the 4 tet vertices, the plane normal `n` and offset `d` (with `n̂·X = d`), and the
front direction `u`.

- The default plane is `x + y - z = 0.5`, i.e. `n = (1,1,-1)/√3`, `d = 0.5/√3`. It
  separates {v0, v3} from {v1, v2}, so it cuts the tet in a quad. This is checked and
  printed at startup.
- The default `u = (1,-1,0)/√2`. The quad's vertices, and the `t` at which the front
  passes each one, are:
  - D = (0, 0.75, 0.25) on edge v2-v3, at t = 0;
  - B = (0, 0.5, 0) on edge v0-v2, at t = 1/6;
  - A = (0.5, 0, 0) on edge v0-v1, at t = 5/6;
  - C = (0.75, 0, 0.25) on edge v1-v3, at t = 1.
- The cut starts at the single vertex D and grows towards the vertex C. One end of the
  front line is always on the face x+y+z = 1. The other end is on the face x = 0 for
  t < 1/6, on z = 0 for 1/6 < t < 5/6 (this includes the default t = 0.6), and on y = 0
  for t > 5/6.
- At t = 1 the front line only touches the polygon, at C. This is why g type (b) stays
  below 1 at t = 1 (see the caveats).
- The plane∩tet polygon is clipped by the half-plane `X·u <= s(t)` and fan-triangulated.
- Zero-area triangles are skipped. Every triangle is wound so that `((b-a)×(c-a))·n > 0`.

**Solver (Zombie WoSt).**

- The four tet faces, oriented outward, are the reflecting boundary. The cut triangles
  are the absorbing boundary.
- The PDE is Laplace with zero source and zero Neumann data.
- `solveDoubleSided = true`, because the cut is an open two-sided surface.
- Each solve runs in parallel over query points with Zombie's TBB `parallel_for`. All
  walks of one point run on one thread.
- Settings, all scaled by `L`:
  - walk termination shell 1e-4·L. Zombie is given distances to the cut that are smaller
    than the true ones by a face-offset margin of 4.3e-5·L, so walks stop within
    1.43e-4·L of the cut (see "Epsilons and the face-offset margin");
  - minimum star radius 1e-5·L;
  - walk-length cap 4096.
- Zombie's default escape test only checks the bounding box. It is replaced by exact
  tet containment, because walks leaking through the slanted face otherwise kept walking
  outside the tet.
- Walks that leave the tet through float round-off are discarded by Zombie and reported
  as "escaped". This is about 0.13% of walks. The worst point is the lattice point at the
  apex (0,0,1), at 5.5%.
- Walks that hit the cap are reported as "capped". They still count, with value 0. None
  hit the cap at the default settings.

**Which side, +g or -g.**

- Zombie calls the Dirichlet callback with the walk's end point *projected onto the cut*,
  so `sign((x-p)·n)` at that point would be meaningless.
- It also passes a flag. The flag is true when the unprojected end point lies on the side
  that the nearest cut triangle's normal points to.
- Since every cut triangle is wound with normal `+n`, the flag means "the walk ended on
  the `n·X > d` side". The callback returns `+g(x)` if the flag is set and `-g(x)`
  otherwise, with `g` evaluated at the projected point.

**Epsilons and the face-offset margin.**

- The *user epsilon* (UI) excludes lattice points near the cut and sets the jump-check
  offset.
- The *walk termination shell* is a fixed 1e-4·L, at least 10x smaller. If it were not,
  the jump-check points would terminate immediately and the check would be trivially exact.
- **Face-offset margin.** A step that starts on a tet face starts from a point that Zombie
  pushes slightly off the face (`offsetPointAlongDirection`: 2^-16, or 256 float ulps, per
  coordinate, about 1.5e-5 here). The step's ball radius, however, is 0.99x the distance
  to the cut measured at the un-pushed point. Close to the cut the 1% safety margin (about
  1.4e-6 at the shell) is smaller than the push. So a walk that is on a face within about
  6e-4·L of the cut can step across the cut and end on the wrong side.
- To prevent this, the solver gives Zombie distances to the cut that are smaller by twice
  the largest possible push for the tet's coordinates (6.1e-5 = 4.3e-5·L for the default
  tet). Every ball then stays on its own side of the cut. Walks stop within shell + margin
  = 1.43e-4·L. The measurements that led to this are under
  [Fixes made during validation](#fixes-made-during-validation).
- The minimum star radius is kept below the shell. Zombie shrinks the star radius to 0.99x
  the distance to the cut only while that distance is at least the minimum star radius;
  below that the radius equals the distance.

**t = 0.** With no Dirichlet boundary, Laplace with zero Neumann data on every face fixes
H only up to a constant. No walks are run and H is set to 0. The UI says so.

**Gradient.** ∇H comes from Zombie's gradient estimator (control variates plus antithetic
pairs).

- The estimator starts with a ball of radius min(distance to cut, distance to faces), and
  its noise grows like 1/radius.
- ∇H is therefore estimated only at *interior* lattice points, those not on a tet face,
  edge or vertex. On the faces it would be useless (s.e. about 600 at N = 256). This is a
  deviation from the spec, which asks for ∇H at the lattice points; 277 of the 433 default
  lattice points are on a face.
- Even at interior points it is noisy at N = 256: in `--bench`, median |∇H| is 1.66 and
  the median standard error is 1.59. Raise N to see a clean field.

**Query points.**

- Every lattice and slice point is pulled 0.1% of the way towards the centroid, a uniform
  shrink, so that walks start strictly inside. Points exactly on a slanted face otherwise
  lose all their walks to round-off.
- The slice is a square grid (64 cells across) over the plane∩tet section. Cells with a
  corner inside are kept, and their outside corners are snapped onto the section boundary.

**Random numbers.** Point `i` of a solve uses its own pcg32 stream `i`, with a state
derived from the solve's seed. This replaces Zombie's clock-based seeding. The validation
is therefore reproducible for a given build and base seed, and every point, repeat and
solve is independent. The seeds are recorded in `results.txt`.

## Validation

`--validate` (or the GUI button) prints the following and writes the same text to
`results.txt`, together with all settings and seeds.

1. **Convergence.** Four fixed probes near or ahead of the front at the default t = 0.6.
   Walks there end on both sides of the cut, so the estimate is genuinely random. For
   N = 16, 64, 256, 1024 and 4096, with 8 independent repeats each, it prints:
   - mean and std of the 8 repeat means;
   - std·√N;
   - the std predicted from the per-walk variance;
   - the number of distinct repeat values;
   - the log-log slope of std against N (1/√N gives −0.5).

   This is done for both g types.
2. **Complete-cut limit, t = 1.** For each side of the cut, over all lattice points, both
   g types:
   - max and mean |H − sign((X−p)·n)|;
   - the number of points above 0.1;
   - max |H|;
   - for g = 1 (when no walk hit the cap), the exact number of walks that ended on the
     other side of the cut.
3. **Jump.** 20 points X sampled uniformly on the cut (t = 0.6) at least max(0.05 L, 2ε)
   from its edges, evaluated at X ± ε·n with N = 4096. Per pair it prints H+, H−, the jump
   and 2g(X), then the mean jump, mean 2g and mean jump/2g. Both g types.
   - The deepest point of the default t = 0.6 cut is only 0.108·L from its edges. The
     margin is capped at 0.8x that depth (0.087·L), so for ε/L above about 0.043 the cap
     applies, the output says so, and all X lie near the cut's deepest point.
   - If no admissible X exists, the check prints that it was not run instead of printing
     means.
4. **Progressive sweep.** For t = 0.2, 0.4, 0.6, 0.8 and 1.0, the fraction of lattice
   points with |H| > 0.9, for both g types, plus the maximum of g_b on the cut.

## Validation results

All numbers below come from `results.txt`, which is the source of truth. It was written by
`--validate` on 2026-10-31 at 15:54. The binary was an up-to-date Release build (MSVC
19.44, 16 threads). The run used the defaults: t = 0.6, N = 256, ε = 0.01·L = 0.01414,
k = 12 and base seed 20261031. Over all solves, 2,001,409 of 2,003,968 walks were
counted. 2,559 (0.128%) escaped and were discarded, and none hit the cap. Wall time was
4.64 s. A second run gave identical H, std, counts and walk tallies; only the date and
the timings changed.

| check | key number | agrees with theory? |
|---|---|---|
| 1. convergence, both g | mean log-log slope −0.488 over 8 fits (range −0.64 to −0.39) | yes |
| 2. t = 1, g = 1 | max \|H ∓ 1\| = 0 on both sides; 0 of 107,651 walks on the wrong side | yes, exactly ±1 |
| 2. t = 1, g_b | max \|H\| = 0.811, below the bound 0.866 | yes, but no ±1 limit (by construction) |
| 3. jump, g = 1 / g_b | mean jump/2g = 0.9887 / 0.9844 | yes, the shortfall is the O(ε) offset |
| 4. sweep, g = 1 | fraction with \|H\| > 0.9 rises from 0.002 to 1.000 | yes |
| 4. sweep, g_b | 0 at every t | yes, but the test cannot detect a front here |
| extra: independent simulator | 8 of 8 probe values within 1 s.e. of Zombie | yes |
| solve time at defaults | 1.05 s | under the ~2 s target |

No check disagreed with theory. The one stated expectation that does not hold is the
spec's ±1 limit for g type (b). That is a property of the spec's g_b, not a solver error
(see the caveats).

### 1. Convergence

Std of the 8 repeat means, at the smallest and largest N:

| probe | g = 1, N = 16 | g = 1, N = 4096 | g_b, N = 16 | g_b, N = 4096 |
|---|---|---|---|---|
| P1 | 0.135 | 0.0175 | 0.057 | 0.0050 |
| P2 | 0.160 | 0.0080 | 0.040 | 0.0058 |
| P3 | 0.133 | 0.0137 | 0.049 | 0.0034 |
| P4 | 0.229 | 0.0140 | 0.080 | 0.0015 |

- **Slopes.** The log-log slopes are −0.387, −0.503, −0.463 and −0.532 for g = 1, and
  −0.460, −0.430, −0.493 and −0.638 for g_b. Their mean is −0.488 with a standard error of
  0.022, within 0.6 s.e. of −0.5.
- **Scatter.** The std from 8 repeats is only known to about ±27%, so a single slope has an
  s.e. of about 0.06. The 8 slopes have an sd of 0.075. The extremes, −0.39 and −0.64, are
  1.9 and 2.3 s.e. from −0.5.
- **Measured against predicted.** The ratio of measured std to the std predicted from the
  per-walk variance has a median of 0.92 over the 40 (probe, N, g) cells. 5 of the 40 fall
  outside the 90% band [0.56, 1.42], where 4 are expected. So the repeats are
  independent: correlated RNG streams would give a std well below the prediction.
- **Distinct values.** At N = 16 only 4 of 8 repeat values are distinct at some probes.
  This is quantization: with g = 1 and 16 walks, H is a multiple of 1/8. It does not mean
  duplicate streams.
- **Means.** The means are stable across N. For example, P3 with g = 1 stays between 0.797
  and 0.830.

The estimator behaves like plain Monte Carlo, with error ∝ 1/√N.

### 2. Complete-cut limit (t = 1)

The lattice has 455 points. 34 lie within ε of the cut and are excluded, which leaves 126
points on the + side and 295 on the − side.

- **g = 1:** max |H − 1| = 0 and max |H + 1| = 0. All 107,651 counted walks ended on the
  side they started on, none hit the cap, and 125 (0.116%) escaped and were discarded.
  At t = 1 the sides are disconnected, so this is exactly what theory requires. There is
  no ε bias, because the exact solution is constant on each side.
- **g_b:** H does not approach ±1. On the + side H lies between 0.281 and 0.811 (mean
  |H − 1| = 0.419). On the − side max |H| = 0.811 and mean |H + 1| = 0.404. At t = 1 the
  front line still touches the section at C, and the farthest cut point is only 0.75·L
  from it. So g_b ≤ sqrt(0.75) = 0.866 on the whole cut, and by the maximum principle
  |H| ≤ 0.866. The results satisfy this bound.

### 3. Jump across the cut (t = 0.6, X ± ε·n, N = 4096)

- **g = 1:** mean jump 1.9773 against 2g = 2, so mean jump/2g = 0.9887 (per pair 0.9687
  to 0.9980).
- **g_b:** mean jump 0.8396 against a mean 2g of 0.8526, so mean jump/2g = 0.9844 (per
  pair 0.9773 to 0.9890).
- **The shortfall is real.** Per-point standard errors are 0.0006 to 0.004, and the noise
  in the mean ratio is about 0.0007, far below the 1.1 to 1.6% shortfall.
- **The shortfall is expected.** The check measures at a distance ε from the cut, where
  H(X ± ε·n) ≈ ±g ∓ ε·|∂H/∂n|. A 1.1% shortfall for g = 1 means a mean |∂H/∂n| of about
  0.8 on each side.
- **It grows towards the front.** There ∂H/∂n is singular (∝ r^−1/2). The lowest g = 1
  ratios are at the pairs nearest the front: 0.969 at r/L = 0.083, 0.979 at 0.092 and
  0.977 at 0.096. The highest is 0.998 at r/L = 0.308.
- **Termination contributes little.** Walks stop within 1.43e-4·L of the cut, 70x closer
  than ε.

The solver reproduces the prescribed jump 2g for both data types, up to the expected O(ε)
offset.

### 4. Progressive sweep (fraction of lattice points with |H| > 0.9)

| t | 0.2 | 0.4 | 0.6 | 0.8 | 1.0 |
|---|---|---|---|---|---|
| (a) g = 1 | 0.0022 (1/449) | 0.0975 (43/441) | 0.1894 (82/433) | 0.7330 (313/427) | 1.0000 (421/421) |
| (b) g_b | 0 | 0 | 0 | 0 | 0 |
| max g_b on the cut | 0.387 | 0.548 | 0.671 | 0.775 | 0.866 |
| cut area (of 0.271) | 0.038 | 0.103 | 0.168 | 0.233 | 0.271 |

- **g = 1.** The fraction never decreases and reaches exactly 1 at t = 1. The big step
  between t = 0.6 and 0.8 comes from the geometry. At t = 0.8 the cut already covers 86%
  of the section, and only a thin uncut ligament near A and C connects the two sides.
- **g_b.** The fraction must be 0: max g_b < 0.9 at every t, and |H| ≤ max g_b. The
  |H| > 0.9 test therefore cannot show a saturation front for g_b with this g_b and
  geometry.

### Extra: independent cross-check (not part of the spec)

None of the spec's checks compares H at a partial cut with an independent method; check 2
is trivially piecewise constant. So a separate numpy reflected-Brownian-motion simulator
was written that shares no code with Zombie. It uses:

- Gaussian steps with σ = clip(0.25·dist, 5e-4, 0.02);
- mirror reflection on the tet faces;
- absorption when a step crosses the cut polygon, with a Brownian-bridge correction;
- 40,000 walkers per probe.

The Zombie values are the check-1 means pooled over all N and repeats (43,648 walks per
probe).

| probe | g = 1, simulator | g = 1, Zombie | g_b, simulator | g_b, Zombie |
|---|---|---|---|---|
| P1 | −0.1490 | −0.1527 | −0.0568 | −0.0581 |
| P2 | +0.4247 | +0.4267 | +0.0849 | +0.0856 |
| P3 | +0.8247 | +0.8229 | +0.2012 | +0.2019 |
| P4 | −0.3482 | −0.3546 | −0.1180 | −0.1164 |

All 8 differences are below 1 combined s.e. (the largest is 0.98 s.e.). With a coarser
step, all are within 1.75 s.e. So Zombie's H at the partial cut is unbiased to within
about 0.005 for g = 1 and 0.0015 for g_b. This covers the reflecting faces, the ±g side
assignment and data type (b). The script was a one-off and is not in this repository, and
these numbers are not in `results.txt`.

### Solve time

`--bench` ran one full default solve three times. Totals were 1.073, 1.042 and 1.052 s,
median 1.052 s. Of that, the lattice took 0.25 to 0.28 s and the slice 0.79 to 0.80 s. The
solve covers 433 lattice points (156 with ∇H) and 1481 slice points, at N = 256.

Walk counts were identical across the three runs:

- lattice: 110,682 counted of 110,848;
- slice: 378,595 counted of 379,136;
- none capped, and an escape rate of 0.14 to 0.15%.

The machine was the Release build, MSVC 19.44, 16 hardware threads. That is under the
spec's ~2 s target.

### Screenshots

Made with `--screenshot` with a real OpenGL 3.3 window (`screenshots/`):

| t = 0.6 | t = 1.0 |
|---|---|
| ![H at t = 0.6](screenshots/t060.png) | ![H at t = 1](screenshots/t100.png) |

- **t = 0.6.** The slice shows the red (+) / blue (−) split ending at the front, with a
  smooth transition and ∇H arrows ahead of it.
- **t = 1.** The two sides are uniformly +1 and −1, and the arrows have zero length,
  consistent with check 2.
- **`screenshots/t060_ui.png`.** The same view with the control panel and status lines.

### Caveats

#### Fixes made during validation

- **Final validation pass.** This `results.txt` needed no code change. One README number
  was corrected: the median measured/predicted std ratio, from 0.93 to 0.92.
- **Wrong-side walks.** An earlier validation run found walks ending on the wrong side of
  the cut at t = 1:
  - 10 of 6.89M walks at N = 4096;
  - at the default N, one wrong walk, which gave max |H − 1| = 0.0078.

  Instrumenting such walks showed that each started its last step on a tet face, close to
  the cut. The cause is the face push-off described above. The fix is the face-offset
  margin. After it, 0 of 13.8M walks at N = 4096 and 0 of 107,651 here ended on the wrong
  side, with the same mean walk length and run time.
- **Minimum star radius.** In the same investigation it was lowered from 1e-3·L to
  1e-5·L. Before the margin existed, the 1e-3·L value put 91 of 1.72M walks on the wrong
  side.
- **First smoke tests.** These fixes came earlier:
  - an exact tet-containment escape test (the bounding-box test gave H = 0.911 instead of
    1 at the apex);
  - `solveDoubleSided = true` (otherwise H = −g everywhere);
  - the 0.1% pull of query points (points on the slanted face lost all their walks);
  - ∇H only at interior points;
  - the t = 0 special case;
  - per-point pcg32 streams instead of clock seeding.

#### Zombie, not a plain WoSt

All results validate Zombie's `WalkOnStars` (commit `b232b1c`) as configured here, not an
independent implementation.

- Zombie works in float32, with FCPW BVH queries for distances and ray hits. The 0.13% of
  walks that escape come from that round-off.
- This code overrides four things: the escape test, the distance to the cut (shrunk by the
  margin), the seeding, and the side-dependent ±g callback.
- ∇H uses Zombie's control-variate and antithetic estimator.
- No plain textbook WoSt was written for comparison. The independent check is the
  Brownian-motion simulator, which checks the PDE solution, not WoSt as such.

#### g type (b) at t = 1

The spec's g_b = min(1, sqrt(r/L)) uses r from a front line that still exists at t = 1:
it touches the section at C, and r ≤ 0.75·L on the cut.

- So g_b ≤ 0.866 and |H| ≤ 0.866 everywhere.
- Check 2(b) can never approach ±1, and check 4(b) can never exceed 0.9.
- A ±1 limit for (b) would need a different definition, for example a length scale well
  below 0.75·L, or g := 1 once the cut is complete.

#### ε bias

- **Jump check.** Measuring at X ± ε·n understates the jump by ε·(|∂H+/∂n| + |∂H−/∂n|):
  1.1% (g = 1) and 1.6% (g_b) at ε = 0.01·L. The earlier run measured jump/2g = 0.932 at
  ε = 0.06·L and 0.895 at ε = 0.1·L. That is roughly linear in ε, but it was not
  rerun for this `results.txt`, and at those ε the margin cap moves all X near the cut's
  deepest point.
- **Lattice exclusion.** Lattice points within ε of the cut are excluded, so checks 2 and
  4 never test the near-cut region.
- **Termination.** Walks stop within 1.43e-4·L of the cut and take ±g at the projected
  point. This adds an O(1.4e-4·L) bias.

#### Other limitations

- **Escapes.** Escaped walks (0.128% overall, 5.5% at the apex lattice point) are
  discarded, so H at such points is a conditional average. The size of that bias was not
  measured separately. Check 2a stays exact despite 125 escapes, but it cannot show such a
  bias, because H is constant on each side there.
- **∇H.** It is estimated only at the 156 interior lattice points, and it is noisy at
  N = 256 (median s.e. ≈ median |∇H|).
- **Cosmetic.** When the max deviation in check 2 is exactly 0, the line prints
  "at (0.0000, 0.0000, 0.0000)". That is the initial value of the location, not a lattice
  point; the deviation itself is correct.
- **Platforms.** The Linux and macOS builds have not been tested; only Windows/MSVC has.

## Files

- `CMakeLists.txt`: the project and its pinned FetchContent dependencies.
- `src/cut_problem.h`, `src/cut_problem.cpp`:
  - the geometry config;
  - tet and cut geometry;
  - the Zombie solver wrapper;
  - the lattice and slice;
  - the full-scene solve.
- `src/validation.h`, `src/validation.cpp`: the four checks.
- `src/main.cpp`: the CLI modes and the Polyscope GUI.
- `results.txt`: output of the validation run reported above.
- `screenshots/`: `t060.png`, `t100.png` and `t060_ui.png`.
