# The beam's mesh

The tet mesh of the `cutting` project's GXFEM beam example (`examples/gxfem_beam/gxfem_beam.cpp`,
`cutting_example_gxfem_beam`), reproduced for `xfem_beam` and the self-tests.

The example meshes `geometry::make_box(-size / 2, size / 2)` with `size = (1, 0.2, 0.2)`:
`meshing::tetrahedralize` with `max_volume = regular_tet_volume(0.15)`, `min_dihedral_angle = 10` and every
other option at its default (radius-edge ratio 1.414, Steiner points allowed on the surface).

- `beam_surface.off`: `make_box`'s eight corners and twelve triangles, in its order (`geometry/primitives.cpp`).
- `beam.node`, `beam.ele` (TetGen's format, 0-based), `beam.vtk`, `beam_surface.obj`, `beam.json` (the manifest):
  written by `cutting`'s `tetmesh` tool, from this folder:

```
tetmesh beam_surface.off --keep-scale --edge 0.15 --min-dihedral 10 --out . --name beam
```

`--keep-scale` leaves the coordinates as they are, `--edge 0.15` gives `max_volume = regular_tet_volume(0.15)` (the
same function), and `--min-dihedral 10` the dihedral bound; the radius-edge ratio defaults to 1.414 in both. TetGen is
deterministic, so the same input and switches give the same mesh. Every number the example prints at start-up
(`cutting_example_gxfem_beam --no-window --frames 1`) against tetmesh's:

| | the beam example | tetmesh |
|---|---|---|
| TetGen switches | `-pzq1.414/10a0.000397748Q` | `-pzq1.414/10a0.000397748Q` |
| nodes | 85 | 85 |
| tets | 203 | 203 |
| boundary triangles | 154 | 154 |
| volume | 0.04 (box 0.04) | 0.04 (surface encloses 0.04) |
| smallest tet | 8.67826e-05 | 8.67826e-05 |
| smallest dihedral angle | 25.6609 deg | 25.6609 deg |
| nodes held at x = -0.5 | 9 | 9 |

The dynamics agree as well: uncut, under gravity with the example's material and dt, `xfem_beam`'s tip drop is
0.0029 at frame 0 (the example's 0.0029) and 0.3588 at frame 19 (0.3593); the rest is the mass, lumped here and
consistent in the example's newton.
