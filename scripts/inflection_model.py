#!/usr/bin/env python3
"""Energy model and dataset generator for the architectural inflection point.

Every quantity that study reports is computed here from the SimpleScalar result
files in simulator/results/, the same files that support the branch-prediction
study, so the two papers rest on one set of measurements.  Nothing is typed
into the manuscript by hand: run with --write and this script emits the
fars_data.tex that the manuscript inputs, along with two of its tables.

Usage: python3 scripts/inflection_model.py            (report to stdout)
       python3 scripts/inflection_model.py --write=DIR  (also write LaTeX)

WHAT IS MEASURED (first party, from sim-outorder statistics)
    A_l   accesses per committed instruction at level l in {L1I, L1D, L2, DRAM}
    W     wrong-path instructions executed per committed instruction, i.e.
          sim_total_insn/sim_num_insn - 1.  This is speculation waste as the
          simulator actually saw it, not an estimate.
    R_l   cache line replacements per committed instruction at level l

WHAT IS ASSUMED (stated in the manuscript, swept in the sensitivity analysis)
    word width, the share of instructions that write a register, the per-bit
    cost of destroying stored state, and the per-transition dissipation of each
    logic family.  Each assumption is a named constant below.

THE ARGUMENT the numbers are used for: B, the bits a workload destroys per
instruction, is set by the architecture and the program.  It does not fall when
the device improves.  Every adiabatic logic family dissipates a floor amount
per switching transition, and device research lowers that floor, so for any
workload with B > 0 there is a device efficiency beyond which erasure is the
larger term.  The inflection is therefore structural; the measurements only say
where it lands.
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CONTAINER = os.path.dirname(HERE)

# The per-cell result files live in simulator/results/ of this repository.
# This script is run both from inside the repository and from the manuscript
# working tree that sits beside it, so try the repository layout first and fall
# back to the sibling checkout.
_RESULT_DIRS = [
    os.path.join(CONTAINER, "simulator", "results"),
    os.path.join(os.path.dirname(CONTAINER), "FARS_BPred_container",
                 "FARS_BPred", "simulator", "results"),
]
RESULTS = next((p for p in _RESULT_DIRS if os.path.isdir(p)), _RESULT_DIRS[0])

BENCHES = ["gcc", "go", "ijpeg", "li", "perl", "compress"]
PRED, BUDGET = "tscl", "128KB"

# ---------------------------------------------------------------- assumptions
WORD_BITS = 32          # SimpleScalar PISA general register width
REG_WRITE_FRAC = None   # derived per benchmark: non-store, non-branch commits
L1_LINE_B = 32          # -cache:il1/dl1 line size, from the run configuration
L2_LINE_B = 64          # -cache:dl2 line size, from the run configuration

# Per-bit cost of destroying stored state.  Horowitz's 45 nm survey gives
# ~10 pJ for an 8 KB SRAM access and ~100 pJ for 1 MB; taken over a 32-bit
# access these are the per-bit figures below.  DRAM is his 1.3-2.6 nJ per
# access over a 64-byte line.  These are the physical costs, which are four to
# six orders of magnitude above the Landauer bound also listed.
E_BIT = {                     # joules per bit destroyed
    "RF":   0.30e-15,         # register file, ~10 fJ per 32-bit write
    "L1":   0.31e-12 / 32,    # 8-32 KB SRAM, ~10 pJ per access / 32 bits
    "L2":   3.1e-12 / 32,     # 256 KB-1 MB SRAM, ~100 pJ per access / 32 bits
    "DRAM": 1.3e-9 / (64 * 8),
}
K_B = 1.380649e-23
T_AMBIENT = 300.0
E_LANDAUER = K_B * T_AMBIENT * 0.6931471805599453   # 2.87 zJ at 300 K

# Minimum per-cycle FET dissipation for each logic family, in joules.  These
# are the figures already published in Figure 1 of the branch-prediction paper,
# compiled there from Frank and from the ITRS survey.
E_FET = {
    "CMOS 14nm":     1350e-18,
    "CMOS 10nm":      670e-18,
    "CMOS 7nm":       430e-18,
    "2LAL 180nm":      66e-18,
    "S2LAL 180nm":     95e-18,
    "2LAL 350nm":       2e-18,
}
# Transitions per committed instruction in the logic fabric.  This is the only
# term in the model with no first-party measurement behind it, and the
# crossover is directly proportional to it, so it is calibrated rather than
# assumed: N_t is chosen so that the conventional-CMOS arm reproduces the
# published energy per instruction of an out-of-order core at the same node.
# The same N_t is then held fixed while the device is changed, which is the
# comparison the paper is making -- identical logical work, different physics.
#
# Calibration anchors: reported core energy per committed instruction for a
# wide out-of-order core, excluding caches, at two nodes.  They come from the
# same published compilation as the per-transition figures above, which is
# Figure 1 of the branch-prediction paper: Frank, Tierney and Lewis (2023),
# Hoefflinger's ITRS survey (2011), and Frank et al., ICCD 2020.
#
# Two anchors rather than one, because N_t counts logical work and ought to
# come out the same whichever node it is inferred from; that it does is the
# check on the method.  The manuscript states both as reported ranges and shows
# the conclusion to be insensitive to them across a factor of three.
CAL_ANCHORS = [
    ("CMOS 7nm",   5e-12, 15e-12),
    ("CMOS 14nm", 15e-12, 40e-12),
]


def calibrate_nt():
    """(low, nominal, high) transitions per instruction over all anchors.

    Each anchor gives a range for N_t; the reported band is their union and
    the nominal is the geometric midpoint of that union.
    """
    los, his = [], []
    for fam, lo, hi in CAL_ANCHORS:
        e = E_FET[fam]
        los.append(lo / e)
        his.append(hi / e)
    lo, hi = min(los), max(his)
    return int(lo), int((lo * hi) ** 0.5), int(hi)


def anchor_table():
    """What each anchor implies on its own, for the Methods."""
    out = []
    for fam, lo, hi in CAL_ANCHORS:
        e = E_FET[fam]
        out.append((fam, lo, hi, int(lo / e), int(hi / e)))
    return out


# ---------------------------------------------------------------- recovery --
# eta_m: fraction of memory-hierarchy access energy reclaimed by an
#        energy-recovering, reversible hierarchy.  0 is a conventional
#        hierarchy, 1 is the fully adiabatic reversible memory this paper
#        assumes to be the limiting case.
# eta_r: fraction of outright information destruction removed by reversible
#        microarchitecture -- reversible rename, unwind instead of flush.
# The Landauer bound is the floor that survives even at eta_m = eta_r = 1.
ETA_GRID = [0.0, 0.25, 0.50, 0.75, 0.90, 0.99, 1.0]


def e_arch(d, eta_m=0.0, eta_r=0.0):
    """Architectural energy per instruction still dissipated at (eta_m, eta_r)."""
    return (1.0 - eta_m) * d["E_mem"] + (1.0 - eta_r) * d["E_erase"]


def read_stats(path):
    """name -> float, for every 'name value # comment' line sim-outorder emits."""
    out = {}
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            m = re.match(r"^([A-Za-z_][\w.]*)\s+(-?[\d.]+)\s*#", line)
            if m:
                try:
                    out[m.group(1)] = float(m.group(2))
                except ValueError:
                    pass
    return out


def workload(bench):
    """Everything the model needs for one benchmark, all of it measured."""
    s = read_stats(os.path.join(RESULTS, f"{bench}_{PRED}_{BUDGET}.txt"))
    n = s["sim_num_insn"]

    # register-writing commits: not stores, not branches.  Loads and ALU ops
    # each overwrite one architectural register.
    reg_writes = n - s["sim_num_stores"] - s["sim_num_branches"]

    d = {
        "bench": bench,
        "insn": n,
        "ipc": s["sim_IPC"],
        # access intensity per committed instruction
        "A_L1I": s["il1.accesses"] / n,
        "A_L1D": s["dl1.accesses"] / n,
        "A_L2": s["ul2.accesses"] / n,
        "A_DRAM": s["ul2.misses"] / n,
        # speculation waste, measured
        "W": s["sim_total_insn"] / n - 1.0,
        "W_refs": (s["sim_total_refs"] - s["sim_num_refs"]) / n,
        # line replacements per committed instruction
        "R_L1I": s["il1.replacements"] / n,
        "R_L1D": s["dl1.replacements"] / n,
        "R_L2": s["ul2.replacements"] / n,
        "reg_write_frac": reg_writes / n,
    }

    # ---- B, bits destroyed per committed instruction, by component ----------
    # 1. architectural register overwrite: the previous value is gone
    d["B_reg"] = d["reg_write_frac"] * WORD_BITS
    # 2. wrong-path results that are squashed: computed, then discarded
    d["B_spec"] = d["W"] * WORD_BITS
    # 3. cache lines overwritten on replacement.  A clean eviction still
    #    overwrites the SRAM cells that held it, which is what a reversible
    #    hierarchy would have to uncompute rather than discard.
    d["B_L1I"] = d["R_L1I"] * L1_LINE_B * 8
    d["B_L1D"] = d["R_L1D"] * L1_LINE_B * 8
    d["B_L2"] = d["R_L2"] * L2_LINE_B * 8
    d["B"] = d["B_reg"] + d["B_spec"] + d["B_L1I"] + d["B_L1D"] + d["B_L2"]

    # ---- destruction, split by whether a copy survives ---------------------
    # Register overwrite and speculative squash destroy the only copy: the
    # value is gone.  A cache replacement overwrites a copy whose original
    # survives at the next level, so it destroys redundancy, not information.
    d["B_destroyed"] = d["B_reg"] + d["B_spec"]
    d["B_redundant"] = d["B_L1I"] + d["B_L1D"] + d["B_L2"]

    # ---- energy, narrow and broad ------------------------------------------
    # Narrow: only information destruction, charged at register-file cost.
    # A cache replacement is not charged here.  The physical event that
    # overwrites a cache frame is the write that fills it, and that write is
    # already counted in E_mem, so charging the replacement as well would count
    # the same joules twice.
    d["E_erase"] = d["B_destroyed"] * E_BIT["RF"]
    d["E_landauer"] = d["B"] * E_LANDAUER
    d["E_mem"] = (d["A_L1I"] * E_BIT["L1"] * 32
                  + d["A_L1D"] * E_BIT["L1"] * 32
                  + d["A_L2"] * E_BIT["L2"] * 32
                  + d["A_DRAM"] * E_BIT["DRAM"] * 512)
    # Architectural energy with no recovery anywhere, i.e. eta_m = eta_r = 0.
    # e_arch() below gives it at any recovery level.
    d["E_arch"] = d["E_mem"] + d["E_erase"]
    return d


def crossover(d, n_t):
    """Per-transition dissipation at which logic energy equals erasure energy.

    E_logic = n_t * e_dev  and  E_erase is fixed by the workload, so
    e_dev* = E_erase / n_t.  Any family below e_dev* is erasure-dominated.
    """
    return d["E_erase"] / n_t


def main():
    rows = [workload(b) for b in BENCHES]

    def mean(k):
        return sum(r[k] for r in rows) / len(rows)

    print("=" * 78)
    print("MEASURED, per committed instruction (SPEC95, TAGE-SC-L, 128k budget)")
    print("=" * 78)
    print(f"{'bench':>9} {'IPC':>6} {'A_L1I':>7} {'A_L1D':>7} {'A_L2':>7} "
          f"{'A_DRAM':>8} {'W':>7} {'R_L1I':>7} {'R_L1D':>7} {'R_L2':>7}")
    for r in rows:
        print(f"{r['bench']:>9} {r['ipc']:6.3f} {r['A_L1I']:7.3f} {r['A_L1D']:7.3f} "
              f"{r['A_L2']:7.4f} {r['A_DRAM']:8.5f} {r['W']:7.4f} "
              f"{r['R_L1I']:7.4f} {r['R_L1D']:7.4f} {r['R_L2']:7.5f}")
    print(f"{'mean':>9} {mean('ipc'):6.3f} {mean('A_L1I'):7.3f} {mean('A_L1D'):7.3f} "
          f"{mean('A_L2'):7.4f} {mean('A_DRAM'):8.5f} {mean('W'):7.4f} "
          f"{mean('R_L1I'):7.4f} {mean('R_L1D'):7.4f} {mean('R_L2'):7.5f}")

    print()
    print("=" * 78)
    print("B, BITS DESTROYED PER COMMITTED INSTRUCTION, by component")
    print("=" * 78)
    print(f"{'bench':>9} {'B_reg':>8} {'B_spec':>8} {'B_L1I':>9} {'B_L1D':>8} "
          f"{'B_L2':>8} {'B total':>9}")
    for r in rows:
        print(f"{r['bench']:>9} {r['B_reg']:8.2f} {r['B_spec']:8.2f} "
              f"{r['B_L1I']:9.2f} {r['B_L1D']:8.2f} {r['B_L2']:8.2f} {r['B']:9.2f}")
    print(f"{'mean':>9} {mean('B_reg'):8.2f} {mean('B_spec'):8.2f} "
          f"{mean('B_L1I'):9.2f} {mean('B_L1D'):8.2f} {mean('B_L2'):8.2f} "
          f"{mean('B'):9.2f}")

    print()
    print("=" * 78)
    print("ENERGY PER INSTRUCTION")
    print("=" * 78)
    print(f"{'bench':>9} {'E_erase':>12} {'E_mem':>12} {'E_landauer':>13} "
          f"{'ratio phys/Land':>16}")
    for r in rows:
        print(f"{r['bench']:>9} {r['E_erase']*1e12:9.3f} pJ "
              f"{r['E_mem']*1e12:9.3f} pJ {r['E_landauer']*1e18:10.3f} aJ "
              f"{r['E_erase']/r['E_landauer']:16.3e}")

    print()
    print("=" * 78)
    print("CALIBRATION of N_t against published core energy per instruction")
    print("=" * 78)
    lo, nom, hi = calibrate_nt()

    print()
    print("  energy per instruction each family would spend on the SAME logical")
    print("  work, at the nominal N_t:")
    for k, v in sorted(E_FET.items(), key=lambda x: -x[1]):
        print(f"      {k:>14}: {v*nom*1e12:9.4f} pJ/inst")
    print(f"      {'destroyed':>14}: {mean('E_erase')*1e12:9.4f} pJ/inst")
    print(f"      {'memory':>14}: {mean('E_mem')*1e12:9.4f} pJ/inst")
    print()

    print()
    print("=" * 78)
    print("CROSSOVER: per-transition dissipation at which logic == erasure")
    print("=" * 78)
    print("  calibration anchors:")
    for fam, alo, ahi, nlo, nhi in anchor_table():
        print(f"      {fam:>10}: {alo*1e12:.0f}-{ahi*1e12:.0f} pJ/inst at "
              f"{E_FET[fam]*1e18:.0f} aJ/transition -> N_t {nlo:,}-{nhi:,}")
    print(f"      union: {lo:,} to {hi:,}, nominal {nom:,}")
    print()
    print("  inflection against memory recovery eta_m, at eta_r = 0:")
    print(f"      {'eta_m':>7} {'E_arch pJ':>11} {'e_dev* aJ':>11}   families "
          f"still device-dominated")
    for em in ETA_GRID:
        a = sum(e_arch(r, em) for r in rows) / len(rows)
        star = a / nom
        above = [k for k, v in sorted(E_FET.items(), key=lambda x: -x[1])
                 if v > star]
        print(f"      {em:7.2f} {a*1e12:11.4f} {star*1e18:11.2f}   "
              f"{', '.join(above) if above else 'none'}")

    print()
    print("  per-workload inflection, aJ per transition, at nominal N_t:")
    print(f"      {'bench':>9} {'eta_m=0':>10} {'eta_m=0.9':>11} {'eta_m=1':>10}")
    for r in rows:
        print(f"      {r['bench']:>9} {e_arch(r, 0)/nom*1e18:10.2f} "
              f"{e_arch(r, 0.9)/nom*1e18:11.2f} "
              f"{e_arch(r, 1.0)/nom*1e18:10.3f}")

    print()
    print("  the floor that survives full recovery everywhere "
          "(eta_m = eta_r = 1):")
    lp = sum(r["E_landauer"] for r in rows) / len(rows)
    print(f"      Landauer, {lp*1e21:.0f} zJ/inst -> e_dev* = "
          f"{lp/nom*1e21:.3f} zJ per transition")

    print()
    print("  bits per committed instruction, by what happens to the copy:")
    print(f"      information destroyed, no copy survives : "
          f"{mean('B_destroyed'):6.2f}")
    print(f"      redundant copy overwritten, original survives: "
          f"{mean('B_redundant'):6.2f}")

    print()
    print("  family dissipation for reference (min per-cycle FET, published):")
    for k, v in sorted(E_FET.items(), key=lambda x: -x[1]):
        print(f"      {k:>14}: {v*1e18:8.1f} aJ")

    print()
    print(f"  Landauer bound at {T_AMBIENT:.0f} K: {E_LANDAUER*1e21:.2f} zJ per bit")
    print(f"  mean B = {mean('B'):.1f} bits/inst -> Landauer floor "
          f"{mean('B')*E_LANDAUER*1e21:.1f} zJ/inst, "
          f"physical erasure {mean('E_erase')*1e12:.3f} pJ/inst")

    if any(a == "--write" or a.startswith("--write=") for a in sys.argv):
        paper = manuscript_dir()
        if paper is None:
            print("\n--write: no manuscript directory found, nothing written."
                  "\n  The manuscript tree is not part of this repository;"
                  "\n  name it with --write=DIR.")
        else:
            write_data(paper, rows, mean)
    return rows


def manuscript_dir():
    """Directory --write writes into, or None when it is not present.

    The LaTeX sources are not part of this repository.  --write=DIR names the
    directory explicitly; a bare --write looks for the manuscript beside this
    checkout and does nothing if it is absent.
    """
    for a in sys.argv:
        if a.startswith("--write="):
            return a.split("=", 1)[1]
    d = os.path.join(CONTAINER,
                     "Gregg_Teuscher_FARS_Inflection_Point_08_31_26")
    return d if os.path.isdir(d) else None


def write_data(paper, rows, mean):
    """Emit the pgfplots datasets the manuscript plots."""
    out = os.path.join(paper, "fars_data.tex")
    L = []
    L.append("% fars_data.tex -- every dataset the manuscript plots.\n"
             "% Generated by scripts/inflection_model.py; do not edit by hand.\n"
             "% Regenerate with: python3 scripts/inflection_model.py --write\n")

    def table(name, header, body):
        L.append("\\pgfplotstableread[row sep=" + chr(92) + chr(92)
                 + ",col sep=&]{\n" + header + chr(92) * 2 + "\n" + body
                 + "}" + chr(92) + name)

    body = "".join(
        f"{r['bench']} & {r['B_reg']:.3f} & {r['B_spec']:.3f} & {r['B_L1I']:.3f} "
        f"& {r['B_L1D']:.3f} & {r['B_L2']:.3f} & {r['B']:.3f}"
        + chr(92) * 2 + "\n" for r in rows)
    table("bdecomp", "bench & reg & spec & l1i & l1d & l2 & total", body)

    body = "".join(
        f"{r['bench']} & {r['A_L1I']:.4f} & {r['A_L1D']:.4f} & {r['A_L2']:.4f} "
        f"& {r['A_DRAM']:.5f} & {r['W']:.4f}" + chr(92) * 2 + "\n" for r in rows)
    table("access", "bench & l1i & l1d & l2 & dram & waste", body)

    lo, nom, hi = calibrate_nt()
    e_mean = sum(r["E_erase"] for r in rows) / len(rows)
    m_mean = sum(r["E_mem"] for r in rows) / len(rows)
    a_mean = sum(r["E_arch"] for r in rows) / len(rows)

    # family energies at the calibrated transition count
    body = "".join(
        f"{i} & {k.replace(' ', '~')} & {v*1e18:.1f} & {v*nom*1e12:.5f}"
        + chr(92) * 2 + "\n"
        for i, (k, v) in enumerate(sorted(E_FET.items(), key=lambda x: -x[1])))
    table("families", "idx & family & efet & einst", body)

    # per-workload erasure, memory and crossover
    body = "".join(
        f"{i} & {r['bench']} & {r['E_erase']*1e12:.5f} & {r['E_mem']*1e12:.4f} "
        f"& {crossover(r, nom)*1e18:.3f}" + chr(92) * 2 + "\n"
        for i, r in enumerate(rows))
    table("erasure", "idx & bench & eerase & emem & edevstar", body)

    # eta_m sweep: the inflection as the memory hierarchy is made recoverable.
    # Fine grid, because the interesting behavior is all above eta_m = 0.9.
    # eta_m = 1 is excluded: the figure plots against 1 - eta_m on a log
    # axis, which cannot take a zero coordinate.  The eta_m = 1 limit is the
    # destruction floor and is drawn as an asymptote from \farsEdevStar.
    fine = [i / 100.0 for i in range(0, 100)] + [0.995, 0.999]
    body = "".join(
        f"{em:.4f} & "
        + " & ".join(f"{e_arch(r, em)/nom*1e18:.5f}" for r in rows)
        + f" & {sum(e_arch(r, em) for r in rows)/len(rows)/nom*1e18:.5f}"
        + chr(92) * 2 + "\n" for em in fine)
    table("etagrid",
          "eta & " + " & ".join(r["bench"] for r in rows) + " & mean", body)

    # crossover grid: e_dev*(N_t) per workload, so Figure 4 plots from
    # data rather than from expressions typed into the picture.
    grid = [3000, 10000, 30000, 100000, 300000]
    body = "".join(
        f"{n} & " + " & ".join(f"{e_arch(r, 1.0)/n*1e18:.4f}" for r in rows)
        + chr(92) * 2 + "\n" for n in grid)
    table("crossgrid", "nt & " + " & ".join(r["bench"] for r in rows), body)
    body = "".join(
        f"{n} & " + " & ".join(f"{e_arch(r, 0.0)/n*1e18:.4f}" for r in rows)
        + chr(92) * 2 + "\n" for n in grid)
    table("crossgridbroad",
          "nt & " + " & ".join(r["bench"] for r in rows), body)

    # scalars, so a figure and the prose cannot disagree with the model
    scal = [
        ("NtNom", f"{nom:,}"), ("NtLo", f"{lo:,}"), ("NtHi", f"{hi:,}"),
        ("NtNomRaw", str(nom)), ("NtLoRaw", str(lo)), ("NtHiRaw", str(hi)),
        ("EraseMean", f"{e_mean*1e12:.3f}"),
        ("EraseMin", f"{min(r['E_erase'] for r in rows)*1e12:.3f}"),
        ("EraseMax", f"{max(r['E_erase'] for r in rows)*1e12:.3f}"),
        ("MemMean", f"{m_mean*1e12:.3f}"),
        ("EdevStar", f"{e_mean/nom*1e18:.2f}"),
        ("EdevStarLo", f"{e_mean/hi*1e18:.2f}"),
        ("EdevStarHi", f"{e_mean/lo*1e18:.2f}"),
        ("ArchMean", f"{a_mean*1e12:.3f}"),
        ("EdevStarBroad", f"{a_mean/nom*1e18:.1f}"),
        ("EdevStarBroadLo", f"{a_mean/hi*1e18:.1f}"),
        ("EdevStarBroadHi", f"{a_mean/lo*1e18:.1f}"),
        ("BDestroyed",
         f"{sum(r['B_destroyed'] for r in rows)/len(rows):.1f}"),
        ("BRedundant",
         f"{sum(r['B_redundant'] for r in rows)/len(rows):.1f}"),
        ("Bmean", f"{sum(r['B'] for r in rows)/len(rows):.1f}"),
        ("Bmin", f"{min(r['B'] for r in rows):.1f}"),
        ("Bmax", f"{max(r['B'] for r in rows):.1f}"),
        ("LandauerPerInst",
         f"{sum(r['B'] for r in rows)/len(rows)*E_LANDAUER*1e21:.0f}"),
    ]
    # the eta_m at which each family stops being device-dominated
    for name, key in (("SLAL", "S2LAL 180nm"), ("LALb", "2LAL 180nm"),
                      ("LALa", "2LAL 350nm")):
        target = E_FET[key] * nom
        span = a_mean - e_mean
        em = 1.0 if span <= 0 else max(0.0, min(1.0, (a_mean - target) / span))
        scal.append(("EtaCross" + name, f"{em:.3f}"))
        scal.append(("EtaCrossPct" + name, f"{em*100:.1f}"))

    L.append("% scalars quoted in the text, figures and captions\n"
             + "\n".join(chr(92) + "newcommand{" + chr(92) + "fars" + k
                          + "}{" + v + "}" for k, v in scal))

    with open(out, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n\n".join(L) + "\n")
    print(f"\nwrote {out}")
    write_tables(paper, rows, lo, nom, hi, e_mean, m_mean, a_mean)


GUARD_HEAD = """% {name}.tex -- {number} of the manuscript ({what}).
%
% This file is both a fragment of the manuscript and a document in its own
% right.  Compile it alone (in Overleaf: Menu -> Settings -> Main document ->
% this file, then Recompile) and the PDF it produces, {name}.pdf, is the file
% to upload for {number}.  The file NAME carries the label: the export holds
% the table body only, with no caption and no table number.  The caption for
% {number} lives with its float in fars_body.tex.
%
% Generated by scripts/inflection_model.py; do not edit by hand.
%
{bs}ifx{bs}farsmain{bs}undefined
{bs}documentclass[crop,border=8pt,varwidth=7in]{{standalone}}
{bs}IfFileExists{{tables_preamble.tex}}{{{bs}input{{tables_preamble}}}}{{{bs}input{{../tables_preamble}}}}
{bs}begin{{document}}
{bs}begin{{threeparttable}}
{bs}fi
"""
GUARD_FOOT = ("{bs}ifx{bs}farsmain{bs}undefined\n"
              "{bs}end{{threeparttable}}\n{bs}end{{document}}\n{bs}fi\n")


def write_tables(paper, rows, lo, nom, hi, e_mean, m_mean, a_mean=None):
    """Tables 2 and 3, as guarded standalone fragments."""
    bs = chr(92)
    d = os.path.join(paper, "tables")
    os.makedirs(d, exist_ok=True)

    def emit(name, number, what, body):
        head = GUARD_HEAD.format(name=name, number=number, what=what, bs=bs)
        foot = GUARD_FOOT.format(bs=bs)
        with open(os.path.join(d, name + ".tex"), "w", encoding="utf-8",
                  newline="\n") as f:
            f.write(head + body.rstrip() + "\n" + foot)

    e = bs * 2                                        # LaTeX row terminator
    rule, hl = bs + "toprule", bs + "midrule"
    bot = bs + "bottomrule"

    # ---- Table 2: measured quantities -------------------------------------
    hdr = (rf"{bs}textbf{{Workload}} & {bs}textbf{{IPC}} & "
           rf"$A_{{\mathrm{{L1I}}}}$ & $A_{{\mathrm{{L1D}}}}$ & "
           rf"$A_{{\mathrm{{L2}}}}$ & $A_{{\mathrm{{DRAM}}}}$ & $W$ & "
           rf"$R_{{\mathrm{{L1I}}}}$ & $R_{{\mathrm{{L1D}}}}$ & "
           rf"$R_{{\mathrm{{L2}}}}$ {e}")
    body = "\n".join(
        f"    {bs}texttt{{{r['bench']}}} & {r['ipc']:.3f} & {r['A_L1I']:.3f} & "
        f"{r['A_L1D']:.3f} & {r['A_L2']:.4f} & {r['A_DRAM']:.5f} & "
        f"{r['W']:.4f} & {r['R_L1I']:.4f} & {r['R_L1D']:.4f} & "
        f"{r['R_L2']:.5f} {e}" for r in rows)
    n = len(rows)
    mean = lambda k: sum(r[k] for r in rows) / n
    body += (f"\n    {hl}\n    Mean & {mean('ipc'):.3f} & {mean('A_L1I'):.3f} & "
             f"{mean('A_L1D'):.3f} & {mean('A_L2'):.4f} & {mean('A_DRAM'):.5f} & "
             f"{mean('W'):.4f} & {mean('R_L1I'):.4f} & {mean('R_L1D'):.4f} & "
             f"{mean('R_L2'):.5f} {e}")
    emit("Table_2", "Table 2", "measured per-instruction quantities",
         f"{bs}scriptsize\n"
         f"{bs}begin{{tabular}}{{lccccccccc}}\n{rule}\n{hdr}\n{hl}\n{body}\n"
         f"{bot}\n{bs}end{{tabular}}\n"
         f"{bs}begin{{tablenotes}}\n{bs}footnotesize\n"
         f"{bs}item Each workload was executed for 10 million committed "
         f"instructions. $W$ is wrong-path instructions executed per committed "
         f"instruction, taken as the ratio of executed to committed "
         f"instructions less one.\n"
         f"{bs}end{{tablenotes}}")

    # ---- Table 3: family energies vs measured architectural terms ----------
    rows3 = sorted(E_FET.items(), key=lambda x: -x[1])
    def where(x):
        if x > a_mean:
            return "device"
        if x < e_mean:
            return bs + "textbf{architecture}"
        return bs + "textbf{inside band}"

    body = "\n".join(
        f"    {k} & {v*1e18:,.0f} & {v*nom*1e12:.3f} & {where(v*nom)} {e}"
        for k, v in rows3)
    body += (f"\n    {hl}\n"
             f"    Destroyed outright, narrow & --- & {e_mean*1e12:.3f} & --- {e}\n"
             f"    Plus memory hierarchy, broad & --- & {a_mean*1e12:.3f} & --- {e}")
    emit("Table_3", "Table 3", "family energy against the measured terms",
         f"{bs}begin{{tabular}}{{lccc}}\n{rule}\n"
         f"{bs}textbf{{Logic family}} & {bs}textbf{{aJ per}} & "
         f"{bs}textbf{{pJ per}} & {bs}textbf{{Dominant}} {e}\n"
         f" & {bs}textbf{{transition}} & {bs}textbf{{instruction}} & "
         f"{bs}textbf{{term}} {e}\n{hl}\n{body}\n{bot}\n"
         f"{bs}end{{tabular}}\n"
         f"{bs}begin{{tablenotes}}\n{bs}footnotesize\n"
         f"{bs}item Energy per instruction is the calibrated transition count "
         f"of {nom:,} multiplied by the reported per-transition dissipation, so "
         f"every family is charged for the same logical work. The two rows "
         f"below the rule are measurements and bracket the inflection: a "
         f"family above the broad figure is device-dominated, one below the "
         f"narrow figure is architecture-dominated, and one between them is "
         f"inside the band.\n"
         f"{bs}end{{tablenotes}}")

    print(f"wrote {d}\\Table_2.tex and Table_3.tex")


if __name__ == "__main__":
    main()
