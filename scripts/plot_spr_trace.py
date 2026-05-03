#!/usr/bin/env python3
"""Plot panopt_spr's parsimony / CI / RI trace over optimization iterations.

Usage:
    plot_spr_trace.py <prefix>.trace.tsv [-o out.png]

The trace TSV is written by panopt_spr --apply-moves --reject-on-parsimony-increase
and has one row per pass (plus a 'baseline' row at iter 0 and a 'final' row).

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

    # X = cumulative moves (more meaningful than pass#); Y = P_fitch, CI, RI.
    x = []
    p_fitch = []
    ci = []
    ri = []
    phases = []
    radii = []
    pass_labels = []
    for r in rows:
        try:
            x.append(int(r["cum_applied"]))
        except ValueError:
            # 'final' row: keep last cum value
            x.append(x[-1] if x else 0)
        p_fitch.append(int(r["P_fitch"]))
        ci.append(float(r["CI"]))
        ri.append(float(r["RI"]))
        phases.append(r["phase"])
        radii.append(int(r["radius"]) if r["radius"].isdigit() else 0)
        pass_labels.append(r["pass"])

    fig, axes = plt.subplots(2, 1, figsize=(9, 6.5), sharex=True)

    ax = axes[0]
    ax.plot(x, p_fitch, marker="o", linewidth=1.6, color="C0")
    baseline = p_fitch[0]
    ax.axhline(baseline, color="grey", linestyle=":", linewidth=0.8,
               label=f"baseline = {baseline}")
    ax.set_ylabel("P_fitch (Hartigan parsimony)")
    ax.set_title(args.title or f"SPR optimization — {args.trace_tsv}")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right")

    # Annotate phase changes
    for i in range(1, len(phases)):
        if phases[i] != phases[i - 1]:
            label = f"→ {phases[i]}"
            if radii[i] > 0:
                label += f" (R={radii[i]})"
            ax.axvline(x[i], color="red", linestyle="--", linewidth=0.8, alpha=0.5)
            ax.annotate(label, xy=(x[i], baseline), xytext=(5, -10),
                        textcoords="offset points", color="red",
                        fontsize=8)

    ax2 = axes[1]
    ax2.plot(x, ci, marker="s", linewidth=1.4, color="C1", label="CI (Kluge–Farris)")
    ax2.plot(x, ri, marker="^", linewidth=1.4, color="C2", label="RI (Farris)")
    ax2.set_xlabel("cumulative moves applied")
    ax2.set_ylabel("homoplasy index")
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="lower right")

    # Final summary text
    final_drop = baseline - p_fitch[-1]
    final_pct = 100.0 * final_drop / baseline if baseline else 0.0
    summary = (f"baseline P_fitch = {baseline:,}  →  final = {p_fitch[-1]:,}  "
               f"(Δ = −{final_drop:,}, −{final_pct:.2f}%)\n"
               f"CI: {ci[0]:.4f} → {ci[-1]:.4f}     "
               f"RI: {ri[0]:.4f} → {ri[-1]:.4f}     "
               f"passes: {len(rows)}")
    fig.suptitle(summary, fontsize=9, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    out = args.out or args.trace_tsv + ".png"
    fig.savefig(out, dpi=140)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
