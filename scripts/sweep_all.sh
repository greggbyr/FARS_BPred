#!/bin/bash
# sweep_all.sh — full replication sweep of every predictor configuration.
#
# Runs all 504 cells: 6 predictor configurations (bimod, gshare, gshare+OB,
# gshare+OHT, MBP, TAGE-SC-L) x 6 benchmarks (gcc, go, ijpeg, li, perl,
# compress) x 7 table budgets (16K..1M entries) x 2 reverse-mode
# implementations (twin logic, FRMT), 10M instructions each.
#
# Usage:    scripts/sweep_all.sh [parallel_lanes]     (default 10)
# Output:   simulator/results/<bench>_<pred>_<size>[_frmt].txt
#           (each file is the full simulator output including the stats dump;
#            cells that already contain a stats dump are skipped, so the
#            sweep can be interrupted and resumed)
# Runtime:  roughly 1-2 minutes per cell; ~1 hour at 10 lanes on a modern
#           machine for a fresh run.
#
# The gshare-family history width is log2(budget); TAGE-SC-L uses the
# geometry (l1=2, shift=5, xor=2) evaluated in the paper.
set -u
cd "$(dirname "$0")/../simulator"

OUTDIR=./results
mkdir -p "$OUTDIR"
LANES=${1:-10}

BENCHES=(gcc go ijpeg li perl compress)
SIZE_ORDER=(16KB 32KB 64KB 128KB 256KB 512KB 1MB)
declare -A SIZES=(
    [16KB]=16384 [32KB]=32768 [64KB]=65536 [128KB]=131072
    [256KB]=262144 [512KB]=524288 [1MB]=1048576
)
declare -A HIST=(
    [16KB]=14 [32KB]=15 [64KB]=16 [128KB]=17 [256KB]=18 [512KB]=19 [1MB]=20
)

args_for() {  # pred size_key frmt -> simulator arguments
    local p=$1 sv=${SIZES[$2]} h=${HIST[$2]} f=$3
    case $p in
        bimod)  echo "-max:inst 10000000 -bpred bimod -bpred:bimod $sv $f" ;;
        gshare) echo "-max:inst 10000000 -bpred 2lev -bpred:2lev 1 $sv $h 1 $f" ;;
        ob)     echo "-max:inst 10000000 -bpred ob -bpred:ob 1 $sv $h 1 $f" ;;
        oht)    echo "-max:inst 10000000 -bpred oht -bpred:oht 1 $sv $h 1 $f" ;;
        mbp)    echo "-max:inst 10000000 -bpred mbp -bpred:mbp 1 $sv $h 1 $f" ;;
        tscl)   echo "-max:inst 10000000 -bpred tscl -bpred:tscl 2 $sv 5 2 $f" ;;
    esac
}

LOG="$OUTDIR/sweep_all.log"
echo "Sweep started: $(date)  lanes=$LANES" > "$LOG"
joblist=$(mktemp)

for size_key in "${SIZE_ORDER[@]}"; do
    for frmt in 0 1; do
        for bench in "${BENCHES[@]}"; do
            for pred in bimod gshare ob oht mbp tscl; do
                suffix=""
                [[ "$frmt" == "1" ]] && suffix="_frmt"
                tag="${bench}_${pred}_${size_key}${suffix}"
                txt="$OUTDIR/${tag}.txt"
                if [[ -s "$txt" ]] && grep -q "sim_IPC" "$txt"; then
                    echo "[skip] $tag (already complete)" >> "$LOG"
                    continue
                fi
                echo "./Run.pl -db ./bench.db -dir '$OUTDIR/$tag' -benchmark $bench -sim '$PWD/ss3/sim-outorder' -args '$(args_for $pred $size_key $frmt)' >& '$txt'" >> "$joblist"
            done
        done
    done
done

total=$(wc -l < "$joblist")
echo "$total cells queued" | tee -a "$LOG"
xargs -a "$joblist" -d '\n' -P "$LANES" -I{} bash -c '{}'
rm -f "$joblist"

done_n=$(grep -l "sim_IPC" "$OUTDIR"/*.txt 2>/dev/null | wc -l)
echo "Sweep finished: $(date)  complete cells with stats: $done_n" | tee -a "$LOG"
