#!/usr/bin/env bash
# Launch (or relaunch after preemption) a resumable job on the VM, fully detached from SSH
# (the job runs in its own session with all file descriptors redirected, so ssh returns at once).
# Usage: scripts/vm_run.sh <host> <log name> <command...>
set -euo pipefail
HOST=$1; LOG=$2; shift 2
ssh -n "$HOST" "cd ~/class_baryons && export PATH=\$HOME/.local/bin:\$PATH COBAYA_PACKAGES_PATH=\$HOME/cobaya_packages \
  QF_WORKERS=45 QF_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1; \
  setsid nohup $* > results/likelihood/$LOG 2>&1 < /dev/null > /dev/null 2>&1 & disown; exit 0" > /dev/null 2>&1
