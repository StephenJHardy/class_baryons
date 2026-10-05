#!/usr/bin/env bash
# 2D profile grid in (f_cl, u) around the small-fraction valley (prior-free map of its shape).
# The fixed-fraction profiles have their minimum near f_cl * u ~ 6e-4 for f_cl = 0.03 and 0.01, so
# each new fraction gets couplings at f_cl * u = 2e-4, 4e-4, 6e-4, 9e-4, 1.35e-3 (plus u = 0),
# and the existing f_cl = 0.03 / 0.01 grids are reused. Isotropic kernel first (the physically
# motivated law for opaque clumps), then Thomson on the same starts (--start mirror).
# Resumable; launch with scripts/run_detached.sh results/mix/logs/grid2d.log scripts/run_mix_grid2d.sh
set -uo pipefail
cd "$(dirname "$0")/../src"
run() { uv run python profile_mix.py run "$@" 2>&1 | grep --line-buffered -v '^\[\|^$'; }

# isotropic
run 0.05  isotropic 0 0.004 0.008 0.012 0.018 0.027
run 0.02  isotropic 0 0.01 0.02 0.03 0.045 0.0675
run 0.01  isotropic 0.3 0.45                      # upper wall not reached by u = 0.2 in the profile
run 0.005 isotropic 0 0.04 0.08 0.12 0.18 0.27 0.4 0.6
run 0.003 isotropic 0 0.0667 0.133 0.2 0.3 0.45 0.67 1.0
run 0.001 isotropic 0 0.2 0.4 0.6 0.9 1.35 2.0 3.0
# Thomson on the same grid
run 0.05  thomson 0.004 0.008 0.012 0.018 0.027
run 0.02  thomson 0.01 0.02 0.03 0.045 0.0675
run 0.005 thomson 0.04 0.08 0.12 0.18 0.27
run 0.003 thomson 0.0667 0.133 0.2 0.3 0.45
run 0.001 thomson 0.2 0.4 0.6 0.9 1.35
