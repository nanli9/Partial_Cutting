#pragma once

/// @file
/// The ring's enrichment straight from partial_cut, in this process: the same
/// solve as partial_cut --gate coupling --export-xfem (coupling_gate.h,
/// buildXfemRing), handed over as the ring_export the export file would have
/// held. No file in between; data/ only holds exports for --data and the tests.
///
/// Kept out of the xfem library on purpose: the library knows the file format
/// and nothing of partial_cut, and only this example links the two.

#include "xfem/ring_export.hpp"

namespace xfem {

struct live_ring_options {
  int ring = 2;        ///< check 4's ring 1 (node at the centre) or 2 (node moved off it)
  int crack = 0;       ///< 0, 1, 2: K1, K2, K3
  int quad_level = 1;  ///< uniform refinements of every sub-tet (8 children each)
};

/// Solves check 4's ring with partial_cut's fast method (no GT) and returns it
/// as an export. Prints partial_cut's own lines (timings, self-checks); throws
/// std::runtime_error when the self-checks fail.
ring_export live_ring(const live_ring_options& options);

}  // namespace xfem
