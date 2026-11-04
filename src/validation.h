#pragma once

#include "cut_problem.h"

#include <string>

// Runs the four validation checks (convergence, complete-cut limit, jump, progressive
// sweep), prints them to stdout and writes the same text into the first-gate section of
// resultsPath (writeResultsSection below). N walks, epsilon, k, the walk parameters and the
// base seed come from s; the front positions are fixed by the checks themselves (default
// t = 0.6 for checks 1 and 3, t = 1 for check 2, the sweep for check 4). Returns false if
// the file could not be written.
bool runValidation(const GeometryConfig& cfg, const Settings& s, const std::string& resultsPath);

// results.txt holds one marked section per validation mode, and each mode rewrites only its own:
//   ======== BEGIN <name> ========
//   ...
//   ======== END <name> ========
// The known sections are written in a fixed order (first gate, then reuse); other sections are kept
// after them. A file without any marker (written before sections existed) is taken to be the
// first-gate section. Returns false if the file could not be written.
extern const char* const kResultsSectionGate1;   // "first gate (--validate)"
extern const char* const kResultsSectionReuse;   // "incremental walk reuse (--validate-reuse)"
bool writeResultsSection(const std::string& path, const std::string& section, const std::string& text);
