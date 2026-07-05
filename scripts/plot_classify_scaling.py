#!/usr/bin/env python3
"""
Gráficos del benchmark de clasificación con MEMs (LZ77 vs sr-index).

Uso:
  python3 scripts/plot_classify_scaling.py <classify.csv> [--out-dir=<dir>]

CSV de entrada: el generado por run_classify_bench.sh / bench_classify_compare.
Columnas usadas: dataset, n_bytes, index, us_per_read, mem_us_per_read,
locate_us_per_read, lca_us_per_read, pct_correct, pct_over. Las reps se promedian.

Produce:
  classify_us_per_read.png   tiempo total por read vs tamaño (LZ77 vs sr-index)
  classify_phases.png        desglose apilado MEM/locate/LCA por dataset e índice
  classify_accuracy.png      % correcto + % menos-específico (over) por índice
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path
from statistics import mean

COLORS = {"lz77": "tomato", "sr": "seagreen"}
NAMES = {"lz77": "Índice LZ77", "sr": "sr-index"}
PHASE_COLORS = {"mem": "#8888cc", "locate": "#cc7755", "lca": "#55aa88"}
PHASE_NAMES = {"mem": "Extracción MEMs", "locate": "locate", "lca": "LCA"}


def human_bytes(x, _pos=None):
    if x >= 1e9: return f"{x/1e9:.0f} GB"
    if x >= 1e6: return f"{x/1e6:.0f} MB"
    if x >= 1e3: return f"{x/1e3:.0f} KB"
    return f"{x:.0f} B"


def f(v):
    try: return float(v)
    except (TypeError, ValueError): return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--out-dir", default=None)
    args = ap.parse_args()

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.ticker as ticker
    except ImportError:
        print("ERROR: matplotlib no disponible (usa el venv del proyecto).", file=sys.stderr)
        sys.exit(1)

    csv_path = Path(args.csv)
    with csv_path.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        print("ERROR: CSV vacío", file=sys.stderr); sys.exit(1)

    out_dir = Path(args.out_dir) if args.out_dir else csv_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    # acc[(index, n_bytes)][metric] = [valores por rep]; promediamos al final.
    acc = defaultdict(lambda: defaultdict(list))
    order, seen = [], set()
    for r in rows:
        try:
            n, idx = int(r["n_bytes"]), r["index"]
        except (KeyError, ValueError):
            continue
        for m in ("us_per_read", "mem_us_per_read", "locate_us_per_read",
                  "lca_us_per_read", "pct_correct", "pct_over"):
            v = f(r.get(m))
            if v is not None:
                acc[(idx, n)][m].append(v)
        if n not in seen:
            seen.add(n); order.append((n, r.get("dataset", "").replace("classify_", "")))

    order.sort(key=lambda x: x[0])
    ns = [x[0] for x in order]
    labels = [x[1] for x in order]
    iters = list(range(1, len(ns) + 1))

    def avg(idx, n, m):
        vals = acc.get((idx, n), {}).get(m)
        return mean(vals) if vals else None

    saved = []

    # ── 1. us/read total vs tamaño ───────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 5))
    for idx, color in COLORS.items():
        ys = [avg(idx, n, "us_per_read") for n in ns]
        xs = [i for i, y in zip(iters, ys) if y is not None]
        yv = [y for y in ys if y is not None]
        if yv:
            ax.plot(xs, yv, marker="s", linewidth=2, color=color, label=NAMES[idx])
    ax.set_xticks(iters)
    ax.set_xticklabels([f"{i}\n({human_bytes(n)})" for i, n in zip(iters, ns)], fontsize=8)
    ax.set_xlabel("Iteración (tamaño del dataset)")
    ax.set_ylabel("Tiempo (μs/read)")
    ax.set_title("Pipeline de clasificación con MEMs: tiempo por read")
    ax.legend(); ax.grid(axis="y", linestyle="--", alpha=0.4)
    fig.tight_layout(); p = out_dir / "classify_us_per_read.png"
    fig.savefig(p, dpi=150); plt.close(fig); saved.append(p)

    # ── 2. desglose de fases (barras apiladas) ───────────────────────────────
    present = [i for i in COLORS if any(avg(i, n, "us_per_read") is not None for n in ns)]
    x = list(range(len(ns)))
    width = 0.38
    offsets = [(-(len(present) - 1) / 2 + k) * width for k in range(len(present))]
    fig, ax = plt.subplots(figsize=(11, 5.5))
    for (idx, off) in zip(present, offsets):
        bottoms = [0.0] * len(ns)
        for phase in ("mem", "locate", "lca"):
            ys = [avg(idx, n, f"{phase}_us_per_read") or 0.0 for n in ns]
            xs = [xi + off for xi in x]
            ax.bar(xs, ys, width, bottom=bottoms, color=PHASE_COLORS[phase],
                   label=f"{NAMES[idx]} — {PHASE_NAMES[phase]}")
            bottoms = [b + y for b, y in zip(bottoms, ys)]
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=9)
    ax.set_xlabel("Dataset"); ax.set_ylabel("Tiempo (μs/read)")
    ax.set_title("Desglose por fase: MEM (compartida) / locate / LCA")
    ax.legend(fontsize=7, ncol=2); ax.grid(axis="y", linestyle="--", alpha=0.4)
    fig.tight_layout(); p = out_dir / "classify_phases.png"
    fig.savefig(p, dpi=150); plt.close(fig); saved.append(p)

    # ── 3. precisión (correcto + menos específico) ───────────────────────────
    fig, ax = plt.subplots(figsize=(11, 5.5))
    for (idx, off) in zip(present, offsets):
        corr = [avg(idx, n, "pct_correct") or 0.0 for n in ns]
        over = [avg(idx, n, "pct_over") or 0.0 for n in ns]
        xs = [xi + off for xi in x]
        ax.bar(xs, corr, width, color=COLORS[idx], label=f"{NAMES[idx]} — correcto")
        ax.bar(xs, over, width, bottom=corr, color=COLORS[idx], alpha=0.45,
               label=f"{NAMES[idx]} — menos específico")
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=9)
    ax.set_xlabel("Dataset"); ax.set_ylabel("% de reads")
    ax.set_ylim(0, 100)
    ax.set_title("Precisión vs ground-truth (correcto + LCA menos específico)")
    ax.legend(fontsize=8); ax.grid(axis="y", linestyle="--", alpha=0.4)
    fig.tight_layout(); p = out_dir / "classify_accuracy.png"
    fig.savefig(p, dpi=150); plt.close(fig); saved.append(p)

    print("Guardado:")
    for p in saved:
        print(f"  {p}")


if __name__ == "__main__":
    main()
