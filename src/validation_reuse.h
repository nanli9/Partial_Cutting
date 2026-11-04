#pragma once

#include "cut_problem.h"

#include <string>

// Validation of incremental walk reuse (walk_cache.h): V0 (own WoSt loop vs Zombie), the exact
// checks R1 (replay), R2 (resume = same-seed fresh solve), R3 (resume statistics and resume state),
// R4 (bbox prefilter soundness), and the statistical checks 1-4 of the spec (unbiasedness,
// convergence under reuse, reuse statistics vs strip size, temporal coherence). Prints everything
// and writes it into the reuse section of resultsPath (writeResultsSection). N, epsilon, k, the walk
// parameters and the base seed come from s; t is fixed by the checks. quick: a reduced-size R1-R4
// run only (k = 6, N = 32, 6 steps, both g types), printed, nothing written. Returns false if an
// exact check failed or the file could not be written.
bool runReuseValidation(const GeometryConfig& cfg, const Settings& s, const std::string& resultsPath, bool quick);

// --bench-reuse: the default scene (lattice + slice at the defaults of s), full solve at t = 0.2,
// then 16 incremental steps to t = 1, both g types; per step the update statistics and the time of
// a same-seed full re-solve on the same geometry.
int runReuseBench(const GeometryConfig& cfg, const Settings& s);

// --compare-options: the exact reuse options 1 (frozen seeds, full recompute), 2a (every step stored, bbox
// test), 2b (every step stored, queue keyed by the entry front) and 3 (seed replay, WalkCache::advance), plus
// the control 3q (seed replay with 2b's queue, no stored steps), on (i) the lattice sweep t = 0.2 -> 1 in 16
// steps and (ii) check 3's schedule on the default scene (full solve at t = 0.5, then 4 steps of 1/8, 1/16,
// 1/32, 1/64), both g types, each run from each option's own initial full solve. At every step: per-walk
// bitwise equality with option 1's same-seed walks, stored paths against a fresh recorded solve, per-point H
// and grad H within a worst-case round-off bound of option 1's, 2b soundness against a brute-force test of
// every stored ball, 3q's keys against 2b's; per-pass timings (same clock, same threads, a random order per
// step, validation work outside the timed regions), totals with and without setup, by re-walked share, and
// memory (with the updates' temporaries and OS cross-checks). Prints everything and writes it to outPath
// (never results.txt). quick: lattice k = 6, slice res 16, N = 32. threads > 0 limits TBB
// (tbb::global_control). repeats: complete samples per run (medians and ranges). Returns false if a check
// failed or the file could not be written.
bool runCompareOptions(const GeometryConfig& cfg, const Settings& s, const std::string& outPath, bool quick, int threads,
                       int repeats);
