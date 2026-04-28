#!/usr/bin/env bash
# Run ddis_skim_final on the input file list in parallel chunks, then hadd
# the per-chunk outputs into a single ROOT file.
#
# Usage:
#   ./run_parallel.sh <fileList.txt> <bins.yaml> [N_workers] [extra_skim_args...]
#
# Examples:
#   ./run_parallel.sh data/files.txt data/bins.yaml          # uses nproc
#   ./run_parallel.sh data/files.txt data/bins.yaml 4        # 4 workers
#   ./run_parallel.sh data/files.txt data/bins.yaml 8 --no-eid
#
# Output: DDIS_Combined_output.root in the current directory.
# Per-chunk outputs are written to a tmpdir and merged with hadd, then deleted.

set -euo pipefail

# Same cleanup + build prologue as run.sh: wipe stale build artifacts and
# figures so we always run with a fresh binary.
clear
rm -rf ./build/*
rm -rf ./figs/*

mkdir -p build
make -j2

if [ $# -lt 2 ]; then
    echo "Usage: $0 <fileList.txt> <bins.yaml> [N_workers] [extra_skim_args...]" >&2
    exit 1
fi

FILELIST=$1
BINS=$2
shift 2

if [ $# -ge 1 ] && [[ "$1" =~ ^[0-9]+$ ]]; then
    N=$1
    shift
else
    N=$(nproc)
fi

EXTRA_ARGS=("$@")

if [ ! -f "$FILELIST" ]; then
    echo "ERROR: file list not found: $FILELIST" >&2
    exit 1
fi
if [ ! -f "$BINS" ]; then
    echo "ERROR: bins YAML not found: $BINS" >&2
    exit 1
fi

# Resolve the skim binary. Prefer ./build/ddis_skim_final, fall back to
# whatever is on PATH (eic-shell may have an install elsewhere).
if [ -x "./build/ddis_skim_final" ]; then
    SKIM="$(pwd)/build/ddis_skim_final"
elif command -v ddis_skim_final >/dev/null 2>&1; then
    SKIM="$(command -v ddis_skim_final)"
else
    echo "ERROR: ddis_skim_final not found (looked in ./build/ and PATH)." >&2
    exit 1
fi

# Convert to absolute paths so workers can run from a tmpdir.
FILELIST=$(readlink -f "$FILELIST")
BINS=$(readlink -f "$BINS")
OUT_FINAL="$(pwd)/DDIS_Combined_output.root"

TMPDIR=$(mktemp -d -t ddis_skim_parallel.XXXXXX)
trap 'rm -rf "$TMPDIR"' EXIT

# Drop blank/comment lines from the file list, then split into N roughly-
# equal chunks.
CLEAN_LIST="$TMPDIR/files.clean.txt"
grep -v -E '^\s*(#|$)' "$FILELIST" > "$CLEAN_LIST" || true
TOTAL=$(wc -l < "$CLEAN_LIST")
if [ "$TOTAL" -eq 0 ]; then
    echo "ERROR: file list $FILELIST has no usable entries." >&2
    exit 1
fi
if [ "$N" -gt "$TOTAL" ]; then
    echo "INFO: requested $N workers but only $TOTAL files; clamping N=$TOTAL." >&2
    N=$TOTAL
fi

split -d -n "l/$N" -a 3 "$CLEAN_LIST" "$TMPDIR/chunk_"

echo "INFO: running $N parallel workers over $TOTAL files."
echo "INFO: skim binary: $SKIM"
echo "INFO: extra skim args: ${EXTRA_ARGS[*]:-(none)}"

declare -a pids=()
declare -a outs=()
for chunk in "$TMPDIR"/chunk_*; do
    name=$(basename "$chunk")
    workdir="$TMPDIR/work_$name"
    out="$TMPDIR/out_$name.root"
    mkdir -p "$workdir"
    (
        cd "$workdir"
        "$SKIM" "$chunk" "$BINS" --output "$out" "${EXTRA_ARGS[@]}"
    ) > "$TMPDIR/log_$name.txt" 2>&1 &
    pids+=($!)
    outs+=("$out")
done

# Wait, propagating any non-zero exit as a fatal error after collecting logs.
fail=0
for i in "${!pids[@]}"; do
    if ! wait "${pids[$i]}"; then
        echo "ERROR: worker $i failed; tail of its log:" >&2
        tail -n 20 "$TMPDIR/log_chunk_$(printf '%03d' "$i").txt" >&2 || true
        fail=1
    fi
done
if [ "$fail" -ne 0 ]; then
    echo "ERROR: at least one worker failed; aborting hadd." >&2
    exit 1
fi

# Merge.
echo "INFO: hadd-ing $N chunk outputs -> $OUT_FINAL"
hadd -f "$OUT_FINAL" "${outs[@]}"

echo "INFO: parallel skim done. Output: $OUT_FINAL"

# Run the plotter on the merged output. Same pattern as run.sh.
if [ -x "./build/ddis_plot_final" ]; then
    PLOT="$(pwd)/build/ddis_plot_final"
elif command -v ddis_plot_final >/dev/null 2>&1; then
    PLOT="$(command -v ddis_plot_final)"
else
    echo "WARNING: ddis_plot_final not found (looked in ./build/ and PATH); skipping plotting." >&2
    PLOT=""
fi

if [ -n "$PLOT" ]; then
    echo "INFO: running plotter: $PLOT $OUT_FINAL $BINS"
    "$PLOT" "$OUT_FINAL" "$BINS"
fi
