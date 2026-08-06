# FARS_BPred

Code and benchmarks for Fully Adiabatic, Reversible, and Superscalar (FARS)
branch predictor implementations in the SimpleScalar simulator.

Six reversible predictor configurations are implemented in
`simulator/ss3/bpred.c`: bimodal, gshare, gshare+OB (Outcome Buffer),
gshare+OHT (Outcome History Table), MBP (Mirrored Branch Predictor), and
TAGE-SC-L. Each can run its reverse mode either with twin prediction logic
(a mirrored set of reverse-side tables) or with a Forward-to-Reverse Mapping
Table (FRMT, sized at one quarter of the configured predictor budget) in
place of the twin structures.

## Requirements

A Linux environment (native or WSL2) with `gcc`, `make`, and `perl`.
The simulator is 32-bit-era C; the build flags in `scripts/build.sh` handle
modern-compiler strictness, and `simulator/ss3/compat/` carries a `termio.h`
shim for glibc versions that no longer ship one (applied automatically when
needed).

## Build

```sh
scripts/build.sh
```

This produces `simulator/ss3/sim-outorder`. The script builds `sysprobe`
first (the Makefile uses it to derive endianness defines) and avoids a flex
dependency when building libexo.

## Running a single simulation

From `simulator/`:

```sh
./Run.pl -db ./bench.db -dir out_dir -benchmark gcc \
    -sim "$(pwd)/ss3/sim-outorder" \
    -args "-max:inst 10000000 -bpred tscl -bpred:tscl 2 131072 5 2 0"
```

The `-bpred:tscl` arguments are `<l1size> <budget> <shift_width> <xor> <frmt>`;
the trailing flag selects twin logic (`0`) or FRMT (`1`). The gshare-family
predictors (`2lev`, `ob`, `oht`, `mbp`) take
`<l1size> <budget> <hist_width> <xor> <frmt>` with `hist_width = log2(budget)`.
The stats dump on stderr includes, per predictor, `bpred_dir_rate` (forward
direction-prediction rate) and `reverse_bpred_dir_rate` (reverse rate).

## Replicating the paper's results

```sh
scripts/sweep_all.sh [lanes]        # 420 cells: 6 predictors x 5 benches x
                                    # 7 budgets x {twin, FRMT}; ~1h at 10 lanes
scripts/sweep_perstage.sh [lanes]   # 210 cells: TAGE-SC-L component ablation
                                    # (TAGE / TAGE+SC / full)
scripts/summarize_results.sh        # benchmark-mean IPC and direction rates
scripts/summarize_results.sh frmt   # same for the FRMT cells
scripts/summarize_results.sh [frmt] results_perstage   # ablation summary
```

Sweep outputs are one text file per cell in `simulator/results/` (and
`simulator/results_perstage/` for the ablation), named
`<bench>_<pred>_<size>[_frmt].txt`; completed cells are skipped on re-run, so
sweeps are resumable. The repository ships with the complete result set — to
replicate from scratch, delete (or move aside) the results directories first,
then run the sweeps and compare against the shipped values with
`scripts/summarize_results.sh`.

Reproduction notes:

- Simulation is deterministic for a fixed binary, arguments, and environment.
  SimpleScalar copies the host environment onto the simulated program's
  stack, so differing environment variables can shift IPC in the fourth
  decimal place. When comparing configurations, keep the environment
  identical across arms (the sweep scripts do this).
- The benchmarks known to run correctly are gcc, go, ijpeg, li, perl, and
  compress; the sweeps use these six. The remaining SPEC95 binaries in
  `simulator/bench/` cannot currently run: the floating-point suite
  (tomcatv, swim, su2cor, hydro2d, mgrid, turb3d, fpppp, wave5) and vortex
  are missing their SPEC input files, and m88ksim's simulated memory map
  cannot hold the SPEC workload programs (it needs a rebuild from source
  with a larger memory configuration).

## Component ablation knobs

Two environment variables gate TAGE-SC-L stages at prediction-selection time
only (all tables continue to train identically):

- `BPRED_TSCL_NOSC=1` — disable the statistical-corrector inversion decision
- `BPRED_TSCL_NOLOOP=1` — disable loop-predictor commit selection

Both unset (or `0`) gives the full predictor. `scripts/sweep_perstage.sh`
runs the three ablation arms; see the header comments in that script and in
`simulator/ss3/bpred.c` for the measurement methodology.

## Repository layout

```
scripts/                  build, sweep, and summary scripts (start here)
simulator/Run.pl          per-benchmark run harness (used by the sweeps)
simulator/bench.db        benchmark command database for Run.pl
simulator/bench/          SPEC95 little-endian benchmark binaries
simulator/input/          benchmark inputs
simulator/output/ref/     expected benchmark outputs (Run.pl validation)
simulator/ss3/            SimpleScalar 3.0 source; predictors in bpred.c
simulator/results/        main sweep results (one .txt per cell)
simulator/results_perstage/  component ablation results
```
