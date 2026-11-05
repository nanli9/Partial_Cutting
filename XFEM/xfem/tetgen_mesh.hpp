#pragma once

/// @file
/// A tet mesh from TetGen's two text files, prefix.node and prefix.ele: the
/// format cutting's tetmesh writes (data/beam/ holds one). Comments (#) and
/// blank lines are skipped, the index base is read off the first node, and
/// columns past the ones needed (attributes, markers) are ignored.

#include "xfem/tet_mesh.hpp"

#include <string>

namespace xfem {

/// Throws std::runtime_error naming the file and the line on anything
/// malformed, and on a tet that is not positively oriented (tet_mesh's
/// promise): none is swapped, since a mesh read here must be the mesh written.
tet_mesh read_tetgen_mesh(const std::string& prefix);

}  // namespace xfem
