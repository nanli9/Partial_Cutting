#pragma once

#include "cut_problem.h"

#include <string>

// --check-xi: H is solved in the tet's local coordinates ξ, on the reference tet, and a rest tet
// only pulls the cut back into ξ (pullBackCut). The checks, printed and written into their own
// section of resultsPath:
//   X1 pull-back geometry: for random rest tets and rest-space cuts, the ξ cut at t, mapped back by
//      X = X0 + J ξ, is the rest-space cut at the same t, with the same "+g" side;
//   X2 round trip: the default ξ cut pushed into a skewed rest tet and pulled back returns to round-off;
//   X3 the solve sees only the ξ cut: two different rest tets carrying the same ξ cut give bit-identical
//      cuts, H and ∇H (Zombie, solveScene) and bit-identical walk sums (own loop, WalkCache) over an
//      advancing front.
// N and the seed come from s (N is capped for speed). Returns false if a check fails or the file
// could not be written.
bool runXiCheck(const Settings& s, const std::string& resultsPath);
extern const char* const kResultsSectionXi;   // "local coordinates xi (--check-xi)"
