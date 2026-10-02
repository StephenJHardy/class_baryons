#!/usr/bin/env bash
# Run a command detached from the terminal and from any parent session, so that it
# survives the launching shell (or Claude Code session) exiting.
# Usage: scripts/run_detached.sh <log file> <command...>
# The PID of the detached process group leader is written to <log file>.pid.
set -euo pipefail
LOG=$1; shift
mkdir -p "$(dirname "$LOG")"
setsid nohup bash -c "$*" >> "$LOG" 2>&1 < /dev/null &
echo $! > "$LOG.pid"
disown
echo "started (pid $(cat "$LOG.pid")), log: $LOG"
