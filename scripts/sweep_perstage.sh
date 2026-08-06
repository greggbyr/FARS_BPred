#!/bin/bash
# sweep_perstage.sh — TAGE-SC-L component ablation sweep.
#
# Measures each predictor stage's contribution by re-running the full sweep
# grid with later stages disabled:
#   tage    BPRED_TSCL_NOSC=1 BPRED_TSCL_NOLOOP=1   base bimod + TAGE only
#   tagesc  BPRED_TSCL_NOSC=0 BPRED_TSCL_NOLOOP=1   + statistical corrector
#   tscl    BPRED_TSCL_NOSC=0 BPRED_TSCL_NOLOOP=0   + loop predictor (full)
#
# The disable knobs are selection-only (see bpred.c): the disabled component's
# tables still train identically, but its output never influences the used
# prediction.  All three arms set BOTH environment variables with equal-length
# values: SimpleScalar copies the host environment onto the simulated stack,
# so identical variable layout keeps the arms bit-comparable.  For the same
# reason the tscl arm here can differ from an unknobbed run in the fourth
# decimal of IPC; arm-vs-arm deltas within this sweep are the meaningful
# comparison.
#
# Usage:    scripts/sweep_perstage.sh [parallel_lanes]     (default 10)
# Grid:     3 arms x 6 benchmarks x 7 budgets x 2 reverse modes = 252 cells,
#           10M instructions each.
# Output:   simulator/results_perstage/<bench>_<arm>_<size>[_frmt].txt
#           (complete cells are skipped, so the sweep is resumable)
# Runtime:  roughly 30-45 minutes at 10 lanes on a modern machine.
set -u
cd "$(dirname "$0")/../simulator"

OUTDIR=./results_perstage
mkdir -p "$OUTDIR"
LANES=${1:-10}

BENCHES=(gcc go ijpeg li perl compress)
SIZE_ORDER=(16KB 32KB 64KB 128KB 256KB 512KB 1MB)
declare -A SIZES=(
    [16KB]=16384 [32KB]=32768 [64KB]=65536 [128KB]=131072
    [256KB]=262144 [512KB]=524288 [1MB]=1048576
)
declare -A NOSC=(   [tage]=1 [tagesc]=0 [tscl]=0 )
declare -A NOLOOP=( [tage]=1 [tagesc]=1 [tscl]=0 )

LOG="$OUTDIR/sweep_perstage.log"
echo "Sweep started: $(date)  lanes=$LANES" > "$LOG"
joblist=$(mktemp)

for arm in tage tagesc tscl; do
    for size_key in "${SIZE_ORDER[@]}"; do
        for frmt in 0 1; do
            for bench in "${BENCHES[@]}"; do
                suffix=""
                [[ "$frmt" == "1" ]] && suffix="_frmt"
                tag="${bench}_${arm}_${size_key}${suffix}"
                txt="$OUTDIR/${tag}.txt"
                if [[ -s "$txt" ]] && grep -q "sim_IPC" "$txt"; then
                    echo "[skip] $tag (already complete)" >> "$LOG"
                    continue
                fi
                echo "BPRED_TSCL_NOSC=${NOSC[$arm]} BPRED_TSCL_NOLOOP=${NOLOOP[$arm]} ./Run.pl -db ./bench.db -dir '$OUTDIR/$tag' -benchmark $bench -sim '$PWD/ss3/sim-outorder' -args '-max:inst 10000000 -bpred tscl -bpred:tscl 2 ${SIZES[$size_key]} 5 2 $frmt' >& '$txt'" >> "$joblist"
            done
        done
    done
done

total=$(wc -l < "$joblist")
echo "$total cells queued" | tee -a "$LOG"
xargs -a "$joblist" -d '\n' -P "$LANES" -I{} bash -c '{}'
rm -f "$joblist"

done_n=$(grep -l "sim_IPC" "$OUTDIR"/*.txt 2>/dev/null | wc -l)
echo "Sweep finished: $(date)  complete cells with stats: $done_n/210" | tee -a "$LOG"
