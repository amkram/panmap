#!/usr/bin/env python3
"""Plot panopt_spr's parsimony / CI / RI trace over optimization iterations.

Usage:
    plot_spr_trace.py <prefix>.trace.tsv [-o out.png]

The trace TSV is written by panopt_spr --apply-moves --reject-on-parsimony-increase
and has rows every --progress-every leaves (plus 'baseline' and 'final').

Columns: pass, phase, radius, cum_applied, n_variable, P_fitch, P_min, P_max, CI, RI
"""
import argparse
import sys
import csv

import matplotlib.pyplot as plt


def load_trace(path):
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            rows.append(r)
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("trace_tsv", help="path to <prefix>.trace.tsv")
    ap.add_argument("-o", "--out", default=None,
                    help="output PNG (default: <trace_tsv>.png)")
    ap.add_argument("--title", default=None)
    args = ap.parse_args()

    rows = load_trace(args.trace_tsv)
    if not rows:
        sys.exit(f"empty trace: {args.trace_tsv}")

    x = []
    p_fitch = []
    ci = []
    ri = []
    phases = []
    pass_nums = []
    radii = []
    for r in rows:
        try:
            x.append(int(r["cum_applied"]))
        except ValueError:
            x.append(x[-1] if x else 0)
        p_fitch.append(int(r["P_fitch"]))
        ci.append(float(r["CI"]))
        ri.append(float(r["RI"]))
        phases.append(r["phase"])
        try:
            pass_nums.append(int(r["pass"]))
        except ValueError:
            pass_nums.append(pass_nums[-1] if pass_nums else 0)
        try:
            radii.append(int(r["radius"]))
        except ValueError:
            radii.append(0)

    fig, axes = plt.subplots(2, 1, figsize=(10, 6.5), sharex=True)

    ax = axes[0]
    ax.plot(x, p_fitch, marker=".", linewidth=1.2, color="C0", markersize=4)
    baseline = p_fitch[0]
    ax.axhline(baseline, color="grey", linestyle=":", linewidth=0.8,
               label=f"baseline = {baseline:,}")
    ax.set_ylabel("P_fitch (Hartigan parsimony)")
    ax.set_title(args.title or f"SPR optimization — {args.trace_tsv}")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right")

    # Mark pass boundaries (where pass number changes between rows).
    for i in range(1, len(pass_nums)):
        if pass_nums[i] != pass_nums[i - 1] or phases[i] != phases[i - 1]:
            for a in axes:
                a.axvline(x[i], color="red", linestyle="--", linewidth=0.7, alpha=0.4)
            label = f"pass {pass_nums[i]} ({phases[i]}"
            if radii[i] > 0:
                label += f" R={radii[i]}"
            label += ")"
            ax.annotate(label, xy=(x[i], baseline), xytext=(3, -8),
                        textcoords="offset points", color="red",
                        fontsize=7, rotation=90, va="top")

    ax2 = axes[1]
    ax2.plot(x, ci, marker=".", linewidth=1.2, color="C1", label="CI (Kluge–Farris)",
             markersize=4)
    ax2.plot(x, ri, marker=".", linewidth=1.2, color="C2", label="RI (Farris)",
             markersize=4)
    ax2.set_xlabel("cumulative moves applied")
    ax2.set_ylabel("homoplasy index")
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="lower right")

    final_drop = baseline - p_fitch[-1]
    final_pct = 100.0 * final_drop / baseline if baseline else 0.0
    summary = (f"baseline P_fitch = {baseline:,}  →  final = {p_fitch[-1]:,}  "
               f"(Δ = −{final_drop:,}, −{final_pct:.2f}%)\n"
               f"CI: {ci[0]:.4f} → {ci[-1]:.4f}     "
               f"RI: {ri[0]:.4f} → {ri[-1]:.4f}     "
               f"trace rows: {len(rows)}")
    fig.suptitle(summary, fontsize=9, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    out = args.out or args.trace_tsv + ".png"
    fig.savefig(out, dpi=140)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
