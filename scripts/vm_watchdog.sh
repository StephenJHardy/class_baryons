#!/usr/bin/env bash
# Keep a resumable job running on a spot VM: restart the VM after preemption, relaunch the
# job (it resumes from checkpoints / skips finished points), and copy results back here.
# Stops when every expected result file exists locally.
# Usage: scripts/vm_watchdog.sh <instance> <zone> <log name> "<expected result dirs, space-separated>" <command...>
set -uo pipefail
INST=$1; ZONE=$2; LOG=$3; EXPECT=$4; shift 4
PROJECT=astro-compute
HOST="$INST.$ZONE.$PROJECT"
cd "$(dirname "$0")/.."
done_all() { for d in $EXPECT; do [ -f "results/likelihood/$d/result.json" ] || return 1; done; }
while true; do
  status=$(gcloud compute instances describe "$INST" --zone="$ZONE" --project=$PROJECT --format="value(status)" 2>/dev/null)
  if [ "$status" = "TERMINATED" ]; then
    echo "$(date +%T) VM preempted/stopped: restarting"
    if gcloud compute instances start "$INST" --zone="$ZONE" --project=$PROJECT >/dev/null 2>&1; then
      gcloud compute config-ssh --quiet --project=$PROJECT >/dev/null 2>&1
      for i in $(seq 1 30); do ssh -o ConnectTimeout=10 "$HOST" true 2>/dev/null && break; sleep 10; done
      scripts/vm_run.sh "$HOST" "$LOG" "$@" && echo "$(date +%T) job relaunched"
    else
      echo "$(date +%T) restart failed (no capacity?); retrying in 5 min"; sleep 300; continue
    fi
  elif [ "$status" = "RUNNING" ]; then
    if ! ssh -o ConnectTimeout=15 "$HOST" 'ps -eo cmd | grep -q "[p]ython src/profile_u.py"' 2>/dev/null; then
      # not running: either finished or never started after a reboot
      rsync -a -e ssh "$HOST:class_baryons/results/likelihood/" results/likelihood/ 2>/dev/null
      if done_all; then echo "$(date +%T) all results present; watchdog exiting"; exit 0; fi
      scripts/vm_run.sh "$HOST" "$LOG" "$@" && echo "$(date +%T) job (re)launched"
    fi
  fi
  rsync -a -e ssh "$HOST:class_baryons/results/likelihood/" results/likelihood/ 2>/dev/null
  if done_all; then echo "$(date +%T) all results present; watchdog exiting"; exit 0; fi
  sleep 180
done
