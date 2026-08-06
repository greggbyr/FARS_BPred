#!/bin/bash
# summarize_results.sh — print benchmark-mean metrics from a sweep directory.
#
# For each predictor and table budget, prints the mean over the five
# benchmarks (gcc, go, ijpeg, li, perl) of:
#   IPC      sim_IPC                    (forward simulated IPC)
#   FWD      bpred_dir_rate             (forward direction-prediction rate)
#   REV      reverse_bpred_dir_rate     (reverse direction-prediction rate)
#
# Usage:
#   scripts/summarize_results.sh                 # twin-logic cells, results/
#   scripts/summarize_results.sh frmt            # FRMT cells, results/
#   scripts/summarize_results.sh [frmt] results_perstage
#                                                # ablation arms instead of
#                                                # predictors
# Use this to spot-check a replication run: values should match the paper's
# figure data to within ~1e-4 (the simulated program's stack layout absorbs
# the host environment, which perturbs IPC in the fourth decimal).
#
# The means below cover the five benchmarks the paper's cross-benchmark
# figures average over; compress cells (also swept) are reported per-workload
# in the paper rather than folded into those means.
set -u
cd "$(dirname "$0")/../simulator"

MODE=""
DIR=results
for a in "$@"; do
    case $a in
        frmt) MODE="_frmt" ;;
        *) DIR=$a ;;
    esac
done

if [[ "$DIR" == *perstage* ]]; then
    PREDS=(tage tagesc tscl)
else
    PREDS=(bimod gshare ob oht mbp tscl)
fi
BENCHES=(gcc go ijpeg li perl)
SIZES=(16KB 32KB 64KB 128KB 256KB 512KB 1MB)

get() { grep -m1 "$2" "$1" 2>/dev/null | awk '{print $2}'; }

printf '%-8s %-8s %8s %9s %9s\n' SIZE PRED IPC FWD REV
for size in "${SIZES[@]}"; do
    for pred in "${PREDS[@]}"; do
        ipc=0; fwd=0; rev=0; n=0
        for b in "${BENCHES[@]}"; do
            f="$DIR/${b}_${pred}_${size}${MODE}.txt"
            v=$(get "$f" 'sim_IPC ')
            [[ -n "$v" ]] || continue
            ipc=$(echo "$ipc + $v" | bc -l)
            fwd=$(echo "$fwd + $(get "$f" 'bpred_dir_rate ')" | bc -l)
            rev=$(echo "$rev + $(get "$f" 'reverse_bpred_dir_rate ')" | bc -l)
            n=$((n+1))
        done
        [[ $n -gt 0 ]] || continue
        printf '%-8s %-8s %8.4f %9.5f %9.5f\n' "$size" "$pred" \
            "$(echo "$ipc/$n" | bc -l)" "$(echo "$fwd/$n" | bc -l)" "$(echo "$rev/$n" | bc -l)"
    done
done
