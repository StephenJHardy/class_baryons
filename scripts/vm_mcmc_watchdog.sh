#!/usr/bin/env bash
# Keep N resumable Cobaya MCMC chains running on a spot VM: restart the VM after preemption,
# relaunch any chain that is not running (Cobaya resumes from its output), rsync results back.
# Exits when every chain's log reports convergence.
# The running-chain check skips shell processes: the remote bash -c command line itself
# contains the chain's command text and would otherwise always match.
# Usage: scripts/vm_mcmc_watchdog.sh <instance> <zone> <likelihoods> <n_chains> <threads_per_chain>
set -uo pipefail
INST=$1; ZONE=$2; LIKE=$3; N=$4; THREADS=$5
PROJECT=astro-compute
HOST="$INST.$ZONE.$PROJECT"
cd "$(dirname "$0")/.."
launch_missing() {
  for i in $(seq 0 $((N - 1))); do
    ssh -n -o ConnectTimeout=15 "$HOST" "cd ~/class_baryons && \
      if ! ps -eo comm=,args= | grep -Ev '^(bash|sh|ssh) ' | grep -q '[p]ython src/mcmc_run.py $i $LIKE\$'; then \
        if ! grep -q 'has converged' results/mcmc/$LIKE/log_chain_$i.txt 2>/dev/null; then \
          export PATH=\$HOME/.local/bin:\$PATH COBAYA_PACKAGES_PATH=\$HOME/cobaya_packages OMP_NUM_THREADS=$THREADS \
                 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1; mkdir -p results/mcmc/$LIKE; \
          setsid nohup uv run python src/mcmc_run.py $i $LIKE >> results/mcmc/$LIKE/log_chain_$i.txt 2>&1 < /dev/null & disown; \
          echo launched $i; fi; fi; exit 0" 2>/dev/null
  done
}
all_converged() {
  for i in $(seq 0 $((N - 1))); do
    grep -q "has converged" "results/mcmc/$LIKE/log_chain_$i.txt" 2>/dev/null || return 1
  done
}
while true; do
  status=$(gcloud compute instances describe "$INST" --zone="$ZONE" --project=$PROJECT --format="value(status)" 2>/dev/null)
  if [ "$status" = "TERMINATED" ]; then
    echo "$(date +%T) VM preempted/stopped: restarting"
    if gcloud compute instances start "$INST" --zone="$ZONE" --project=$PROJECT >/dev/null 2>&1; then
      gcloud compute config-ssh --quiet --project=$PROJECT >/dev/null 2>&1
      for i in $(seq 1 30); do ssh -o ConnectTimeout=10 "$HOST" true 2>/dev/null && break; sleep 10; done
    else
      echo "$(date +%T) restart failed (no capacity?); retrying in 5 min"; sleep 300; continue
    fi
  fi
  if [ "$status" = "RUNNING" ] || [ "$status" = "TERMINATED" ]; then
    out=$(launch_missing); [ -n "$out" ] && echo "$(date +%T) $out" | tr '\n' ' ' && echo
    rsync -a -e ssh "$HOST:class_baryons/results/mcmc/" results/mcmc/ 2>/dev/null
    if all_converged; then echo "$(date +%T) all chains converged; watchdog exiting"; exit 0; fi
  fi
  sleep 300
done
