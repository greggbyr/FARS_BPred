#!/usr/bin/env python3
"""compile_results.py — flatten the sweep outputs into tidy CSV files.

Reads every per-cell output file in simulator/results/ (main sweep) and
simulator/results_perstage/ (component ablation), named
<bench>_<pred>_<size>[_frmt].txt, and writes:

    simulator/compiled_results/main_sweep.csv
    simulator/compiled_results/perstage_sweep.csv

One row per simulation cell, long/tidy format (spreadsheet- and
pandas-friendly).  Columns:

    benchmark, predictor, budget_label, budget_entries, mode,
    instructions, ipc,
    fwd_dir_rate, rev_dir_rate, fwd_addr_rate, rev_addr_rate,
    tage_fwd_tag_matches, sc_fwd_inversions, loop_fwd_hits, loop_rev_hits

Component columns are blank for predictors without those structures, and
rate columns are blank when the simulator reported a divide-by-zero (no
updates of that kind occurred).  Rows are sorted, so regenerated files are
byte-stable and diff-able.

Usage: python3 scripts/compile_results.py       (from anywhere in the repo)
"""
import csv
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SIM = os.path.join(HERE, "..", "simulator")
OUT_DIR = os.path.join(SIM, "compiled_results")

CELL_RE = re.compile(
    r"^(?P<bench>[A-Za-z0-9]+)_(?P<pred>[A-Za-z0-9]+)_"
    r"(?P<size>\d+(?:KB|MB))(?P<frmt>_frmt)?\.txt$")

SIZE_ENTRIES = {"16KB": 16384, "32KB": 32768, "64KB": 65536, "128KB": 131072,
                "256KB": 262144, "512KB": 524288, "1MB": 1048576, "0KB": 0}

# stat name -> csv column; names are matched with word boundaries so
# e.g. bpred_dir_rate never matches inside reverse_bpred_dir_rate.
STATS = [
    ("sim_num_insn",            "instructions"),
    ("sim_IPC",                 "ipc"),
    ("bpred_dir_rate",          "fwd_dir_rate"),
    ("reverse_bpred_dir_rate",  "rev_dir_rate"),
    ("bpred_addr_rate",         "fwd_addr_rate"),
    ("reverse_bpred_addr_rate", "rev_addr_rate"),
    ("tage_fwd_tag_matches",    "tage_fwd_tag_matches"),
    ("sc_fwd_inversions",       "sc_fwd_inversions"),
    ("loop_fwd_hits",           "loop_fwd_hits"),
    ("loop_rev_hits",           "loop_rev_hits"),
]

COLUMNS = (["benchmark", "predictor", "budget_label", "budget_entries", "mode"]
           + [col for _, col in STATS])


def parse_cell(path):
    """Extract the stat values from one simulator output file."""
    pats = {stat: re.compile(r"\b" + re.escape(stat) + r"\b\s+(\S+)")
            for stat, _ in STATS}
    vals = {}
    with open(path, errors="replace") as f:
        for line in f:
            line = line.split("#")[0]
            for stat, col in STATS:
                if col in vals:
                    continue
                m = pats[stat].search(line)
                if m:
                    tok = m.group(1)
                    # divide-by-zero stats print "<error:..." -> leave blank
                    vals[col] = tok if not tok.startswith("<") else ""
    return vals


def compile_dir(results_dir, out_csv):
    rows, skipped = [], []
    for name in sorted(os.listdir(results_dir)):
        if not name.endswith(".txt"):
            continue
        m = CELL_RE.match(name)
        if not m:
            skipped.append(name)
            continue
        vals = parse_cell(os.path.join(results_dir, name))
        if "ipc" not in vals:
            skipped.append(name)     # incomplete cell (no stats dump)
            continue
        size = m.group("size")
        rows.append({
            "benchmark": m.group("bench"),
            "predictor": m.group("pred"),
            "budget_label": size,
            "budget_entries": SIZE_ENTRIES.get(size, ""),
            "mode": "frmt" if m.group("frmt") else "twin",
            **{col: vals.get(col, "") for _, col in STATS},
        })
    rows.sort(key=lambda r: (r["benchmark"], r["predictor"],
                             int(r["budget_entries"] or 0), r["mode"]))
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    return len(rows), skipped


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    total = 0
    for sub, out in (("results", "main_sweep.csv"),
                     ("results_perstage", "perstage_sweep.csv")):
        d = os.path.join(SIM, sub)
        if not os.path.isdir(d):
            print(f"skip {sub}/ (not present)")
            continue
        n, skipped = compile_dir(d, os.path.join(OUT_DIR, out))
        total += n
        print(f"{out}: {n} rows from {sub}/")
        for s in skipped:
            print(f"  skipped (not a completed cell): {s}", file=sys.stderr)
    if total == 0:
        sys.exit("no cells found — run the sweeps first (see README)")


if __name__ == "__main__":
    main()
