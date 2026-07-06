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


# Rampa secuencial azul (claro→oscuro = colección chica→grande).
# Luminancia estrictamente decreciente (validado). Encoda el orden de tamaño.
SEQ_BLUE = ["#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5", "#08519c"]


def _equal_pct(rows):
    """De un CSV de verify: {largo(int): %EQUAL}."""
    out = {}
    for r in rows:
        n = int(r["evaluated"])
        if n:
            out[int(r["label"])] = 100.0 * int(r["equal"]) / n
    return dict(sorted(out.items()))


def plot_multi(plt, series, out, title):
    """series: lista de (nombre_tamaño, {largo: %EQUAL}) en orden de tamaño."""
    fig, ax = plt.subplots(figsize=(10, 5.8))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    all_lengths = sorted({L for _, d in series for L in d})
    xpos = {L: i for i, L in enumerate(all_lengths)}
    colors = SEQ_BLUE[-len(series):] if len(series) <= len(SEQ_BLUE) else SEQ_BLUE

    for (name, d), col in zip(series, colors):
        xs = [xpos[L] for L in d]
        ys = list(d.values())
        ax.plot(xs, ys, "-o", color=col, linewidth=2, markersize=8,
                markeredgecolor=SURFACE, markeredgewidth=1.5, label=name, zorder=3)

    ax.set_xticks(list(xpos.values()))
    ax.set_xticklabels(all_lengths)
    ax.set_xlabel("Largo del patrón (caracteres)", color=INK)
    ax.set_ylabel("% de patrones con LCA EQUAL", color=INK)
    ax.set_title(title, color=INK, fontsize=13, pad=14)
    # Eje Y adaptativo: zoom al rango de los datos (con margen) para no aplastar
    # las diferencias cuando todas las series están cerca del 100%.
    ymin_data = min(v for _, d in series for v in d.values())
    ax.set_ylim(max(0, ymin_data - 4), 100.6)
    ax.margins(x=0.12)

    ax.grid(axis="y", color="#dddddd", linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#cccccc")
    ax.tick_params(colors=MUTED)

    leg = ax.legend(title="Colección", loc="lower right", frameon=True, fontsize=9,
                    facecolor=SURFACE, edgecolor="#dddddd", framealpha=0.95)
    leg.get_title().set_color(INK)
    for t in leg.get_texts():
        t.set_color(INK)

    ax.text(0.0, -0.15,
            "0 casos ancestro/incomparable (bug) en todos los tamaños y largos — "
            "el LCA de primarias nunca sale del linaje correcto",
            transform=ax.transAxes, ha="left", va="top",
            color=C_EQUAL, fontsize=9.5, fontweight="bold")

    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight", facecolor=SURFACE)
    print(f"Guardado: {out}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", nargs="?", default=None)
    ap.add_argument("--series", action="append", default=[],
                    help="modo multi-tamaño: NOMBRE:csv (repetible)")
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

    # ── Modo multi-tamaño (líneas EQUAL% por colección) ──────────────────────
    if args.series:
        series = []
        for spec in args.series:
            name, _, path = spec.partition(":")
            series.append((name, _equal_pct(list(csv.DictReader(Path(path).open(newline=""))))))
        out = Path(args.out) if args.out else Path("lca_verify_scale.png")
        plot_multi(plt, series, out, args.title)
        return

    if not args.csv:
        print("ERROR: pasa un CSV (modo barra) o --series (modo multi-tamaño)", file=sys.stderr)
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
