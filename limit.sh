#!/bin/bash

MEM_LIMIT="$1"
shift
CMD=("$@")

if [[ -z "$MEM_LIMIT" || -z "$CMD" ]]; then
  echo "Usage: $0 <mem_limit> <command...>"
  exit 1
fi

MEM_BYTES=$(numfmt --from=iec "$MEM_LIMIT")

CGROUP_ID="memlimit-$$"
CGROUP_PATH="/sys/fs/cgroup/$CGROUP_ID"

sudo mkdir -p "$CGROUP_PATH"

echo "$MEM_BYTES" | sudo tee "$CGROUP_PATH/memory.max" > /dev/null
echo "max"        | sudo tee "$CGROUP_PATH/memory.swap.max" > /dev/null

"${CMD[@]}" &
PID=$!

echo "$PID" | sudo tee "$CGROUP_PATH/cgroup.procs" > /dev/null

wait "$PID"
STATUS=$?

sudo rmdir "$CGROUP_PATH"

echo "[+] Command exited with status $STATUS"
exit "$STATUS"
