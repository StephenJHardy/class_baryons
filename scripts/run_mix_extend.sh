#!/usr/bin/env bash
# Stage 4 of the mixed CDM + clump profiles: extend the coupling grids after the coarse pass.
# The coarse profiles for f_cl = 0.03 and 0.01 were still falling at the largest coupling
# (Delta chi^2 = -0.85 and -1.67 relative to u = 0), and f_cl = 0.1 had not reached a threshold,
# so the grids continue upward in steps of 1.5 until the profile turns and crosses, or the
# calculation fails. f_cl = 1 and 0.5 get two larger couplings to check for a second minimum.
# Resumable; launch with scripts/run_detached.sh results/mix/logs/extend.log scripts/run_mix_extend.sh
set -uo pipefail
cd "$(dirname "$0")/../src"
run() { uv run python profile_mix.py run "$@" 2>&1 | grep --line-buffered -v '^\[\|^$'; }

run 0.03 isotropic 3e-3 --start mirror --seed 7             # rerun of a point that did not converge: new seed, same start
run 0.01 thomson 0.027 0.04 0.06 0.09 0.135 0.2
run 0.03 thomson 0.009 0.0135 0.02 0.03 0.045
run 0.1  thomson 0.0027 0.004 0.006 0.009
run 1    thomson 3e-4 5e-4
run 0.5  thomson 5e-4 8e-4
run 0.01 isotropic 0.027 0.04 0.06 0.09 0.135 0.2 --start mirror
run 0.03 isotropic 0.009 0.0135 0.02 0.03 0.045 --start mirror
run 0.1  isotropic 0.0027 0.004 0.006 0.009 --start mirror
run 1    isotropic 3e-4 5e-4 --start mirror
run 0.5  isotropic 5e-4 8e-4 --start mirror
