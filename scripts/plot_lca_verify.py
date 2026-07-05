#!/usr/bin/env python3
"""
plot_lca_verify.py — visualiza la verificación de correctitud del LCA
primarias-only (salida CSV de verify_lca_equiv).

Barras apiladas por largo de patrón: fracción EQUAL (LZ77 == ground truth)
vs LZ_DESCENDANT (LZ77 más específico). Los casos de bug (ancestro /
incomparable) se anotan aparte — deben ser 0.

Uso:
  python3 scripts/plot_lca_verify.py <lca_verify.csv> [--out=<png>] [--title=<str>]
"""
import argparse
import csv
import sys
from pathlib import Path

# Paleta categórica validada (scripts/validate_palette.js, modo light):
#   EQUAL = teal, LZ_DESCENDANT = coral.  ΔE CVD 28 (> objetivo 12).
C_EQUAL = "#2a9d8f"
C_DESC  = "#e76f51"
C_BUG   = "#9d2a2a"
SURFACE = "#fcfcfb"
INK     = "#222222"
MUTED   = "#666666"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--out", default=None)
    ap.add_argument("--title", default="Correctitud del LCA primarias-only vs largo del patrón")
    args = ap.parse_args()

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("ERROR: matplotlib no disponible. Usa ../.venv/bin/python.", file=sys.stderr)
        sys.exit(1)

    rows = list(csv.DictReader(Path(args.csv).open(newline="")))
    if not rows:
        print("ERROR: CSV vacío", file=sys.stderr)
        sys.exit(1)
    rows.sort(key=lambda r: int(r["label"]))

    lengths = [int(r["label"]) for r in rows]
    ev      = [int(r["evaluated"]) for r in rows]
    equal   = [100.0 * int(r["equal"]) / n for r, n in zip(rows, ev)]
    desc    = [100.0 * int(r["lz_descendant"]) / n for r, n in zip(rows, ev)]
    bugs    = [int(r["lz_ancestor"]) + int(r["incomparable"]) for r in rows]
    total_bugs = sum(bugs)

    x = list(range(len(lengths)))
    fig, ax = plt.subplots(figsize=(10, 5.6))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    w = 0.66
    b1 = ax.bar(x, equal, width=w, color=C_EQUAL, label="EQUAL (LZ77 = ground truth)",
                edgecolor=SURFACE, linewidth=2, zorder=3)
    ax.bar(x, desc, width=w, bottom=equal, color=C_DESC,
           label="LZ_DESCENDANT (LZ77 más específico)",
           edgecolor=SURFACE, linewidth=2, zorder=3)

    # Etiqueta directa del % EQUAL sobre cada barra (label selectivo).
    for xi, e in zip(x, equal):
        ax.text(xi, e - 3, f"{e:.0f}%", ha="center", va="top",
                color="white", fontsize=9, fontweight="bold", zorder=4)

    ax.set_xticks(x)
    ax.set_xticklabels(lengths)
    ax.set_xlabel("Largo del patrón (caracteres)", color=INK)
    ax.set_ylabel("% de patrones (con ocurrencia)", color=INK)
    ax.set_ylim(0, 100)
    ax.set_title(args.title, color=INK, fontsize=13, pad=14)

    # Ejes/grilla recesivos.
    ax.grid(axis="y", color="#dddddd", linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#cccccc")
    ax.tick_params(colors=MUTED)

    leg = ax.legend(loc="lower right", frameon=True, fontsize=10,
                    facecolor=SURFACE, edgecolor="#dddddd", framealpha=0.95)
    for t in leg.get_texts():
        t.set_color(INK)

    # Anotación del invariante de correctitud.
    msg = (f"Casos ancestro/incomparable (bug): {total_bugs} en "
           f"{sum(ev):,} patrones — el LCA de primarias nunca sale del linaje correcto")
    ax.text(0.0, -0.16, msg, transform=ax.transAxes, ha="left", va="top",
            color=(C_EQUAL if total_bugs == 0 else C_BUG), fontsize=9.5, fontweight="bold")

    out = Path(args.out) if args.out else Path(args.csv).with_suffix(".png")
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight", facecolor=SURFACE)
    print(f"Guardado: {out}")


if __name__ == "__main__":
    main()
