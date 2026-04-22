#!/bin/sh

if [ "$(id -u)" -ne 0 ]; then
    exec pkexec /bin/sh "$(readlink -f "$0")" "$@"
fi

SCRIPT_PATH="$(readlink -f "$0")"
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"

cd "$SCRIPT_DIR/systems" || {
    echo "Cannot cd to $SCRIPT_DIR/systems"
    exit 1
}

mkdir -p results || exit 1

for f in *.txt; do
    [ -e "$f" ] || continue
    echo "$f"
    echo "Running $f..."
    "$SCRIPT_DIR/limit.sh" 4G magma "$f" > "results/output-${f%.txt}" 2>&1
done

if [ -n "$PKEXEC_UID" ]; then
    chown -R "$PKEXEC_UID":"$PKEXEC_UID" results
fi