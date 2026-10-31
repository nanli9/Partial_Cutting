#pragma once

#include "cut_problem.h"

#include <string>

// Runs the four validation checks (convergence, complete-cut limit, jump, progressive
// sweep), prints them to stdout and writes the same text to resultsPath. N walks, epsilon,
// k, the walk parameters and the base seed come from s; the front positions are fixed by
// the checks themselves (default t = 0.6 for checks 1 and 3, t = 1 for check 2, the sweep
// for check 4). Returns false if the file could not be written.
bool runValidation(const GeometryConfig& cfg, const Settings& s, const std::string& resultsPath);
