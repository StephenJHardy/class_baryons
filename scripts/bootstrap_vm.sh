#!/usr/bin/env bash
# Bootstrap a fresh Ubuntu 24.04 VM (e.g. GCP spot instance) for the likelihood runs.
# Usage on the VM:  bash bootstrap_vm.sh   (then run jobs with QF_WORKERS/QF_THREADS set)
set -euo pipefail
sudo apt-get update -q
sudo apt-get install -y -q build-essential gfortran git curl rsync
command -v uv >/dev/null || curl -LsSf https://astral.sh/uv/install.sh | sh
export PATH="$HOME/.local/bin:$PATH"
[ -d class_baryons ] || git clone https://github.com/StephenJHardy/class_baryons.git
cd class_baryons
scripts/setup_class.sh
export COBAYA_PACKAGES_PATH="$HOME/cobaya_packages"
uv run cobaya-install planck_2018_lowl.TT planck_2018_lowl.EE planck_2018_highl_plik.TTTEEE_lite_native \
    planck_2018_highl_plik.TTTEEE -p "$COBAYA_PACKAGES_PATH" --skip-global --no-progress-bars
echo "export COBAYA_PACKAGES_PATH=$COBAYA_PACKAGES_PATH" >> ~/.bashrc
echo "Bootstrap complete. Suggested: export QF_WORKERS=45 QF_THREADS=4   (180 vCPUs)"
