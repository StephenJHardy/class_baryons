#!/usr/bin/env bash
# Production queue for the mixed CDM + clump profiles (src/profile_mix.py), run sequentially.
# Every point is resumable, so the whole script can simply be rerun after an interruption.
# Launch detached:  scripts/run_detached.sh results/mix/logs/production.log scripts/run_mix_production.sh
# Stages: 1 repeatability at f_cl = 1; 2 isotropic at f_cl = 1; 3 coarse grids for the smaller
# fractions (Thomson, then isotropic mirrored on the Thomson starts). Refinement points are
# added by hand after looking at the coarse profiles (see the experiment log).
set -uo pipefail
cd "$(dirname "$0")/../src"
run() { uv run python profile_mix.py run "$@" 2>&1 | grep --line-buffered -v '^\[\|^$'; }

# 1. repeatability: other seeds, and independent starts from the LCDM solution
run 1 thomson 1e-4 --seed 1 --tag seed1
run 1 thomson 1e-4 --seed 2 --tag lcdmstart --start lcdm
run 1 thomson 0 --seed 1 --tag seed1
run 1 thomson 1.5e-4 --seed 1 --tag seed1

# 2. isotropic kernel at f_cl = 1, on the same starts and sample points as Thomson
run 1 isotropic 5e-5 1e-4 1.2e-4 1.5e-4 2e-4 --start mirror

# 3. smaller fractions: coarse grids chosen from the fixed-parameter scan (not assuming 1/f_cl)
run 0.5  thomson 0 1e-4 2e-4 2.7e-4 3.6e-4
run 0.1  thomson 0 4.5e-4 9e-4 1.35e-3 1.8e-3
run 0.03 thomson 0 1.5e-3 3e-3 4.5e-3 6e-3
run 0.01 thomson 0 4.5e-3 9e-3 1.35e-2 1.8e-2
run 0.5  isotropic 1e-4 2e-4 2.7e-4 3.6e-4 --start mirror
run 0.1  isotropic 4.5e-4 9e-4 1.35e-3 1.8e-3 --start mirror
run 0.03 isotropic 1.5e-3 3e-3 4.5e-3 6e-3 --start mirror
run 0.01 isotropic 4.5e-3 9e-3 1.35e-2 1.8e-2 --start mirror
