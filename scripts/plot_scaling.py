#!/usr/bin/env python3
"""
Genera gráficos comparativos de escalado (LZ77 vs sr-index) desde scaling.csv.

Uso:
  python3 scripts/plot_scaling.py <scaling.csv> [--out-dir=<dir>]
  python3 scripts/plot_scaling.py <scaling.csv> [--out=<size.png>]   # alias retro-compat

El CSV debe ser el unificado generado por run_scaling_bench.sh.
Columnas usadas: dataset, n_bytes, index, total_bytes, total_bpc, build_time_s

Produce varias imágenes, cada una con LZ77 vs sr-index:
  scaling_size.png         tamaño del índice (bytes, escala log) + línea del texto original
  scaling_bpc.png          bits por carácter (escala lineal)
  scaling_build.png        tiempo de construcción (segundos, escala log)
  scaling_query.png        tiempo de locate extremal (us/query, escala lineal)
  scaling_size_last5.png   barras agrupadas: tamaño (bytes, ESCALA LINEAL) últimos 5 datasets
  scaling_query_last5.png  barras agrupadas: locate (us/query, lineal) últimos 5 datasets

Eje X: número de iteración (1, 2, 3, ...), anotado con el tamaño del dataset.
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

# Estilo común a todas las figuras (consistente con la versión previa).
COLORS = {"lz77": "tomato", "sr": "seagreen"}
NAMES = {"lz77": "Índice LZ77", "sr": "sr-index"}


def read_csv(path: Path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def human_bytes(x, _pos=None):
    if x >= 1e9:
        return f"{x/1e9:.0f} GB"
    if x >= 1e6:
        return f"{x/1e6:.0f} MB"
    if x >= 1e3:
        return f"{x/1e3:.0f} KB"
    return f"{x:.0f} B"


def _float_or_none(value):
    if value in (None, ""):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def plot_metric(plt, ticker, *, ns, labels, by_index, key, title, ylabel,
                yscale, y_formatter, text_line, out_path):
    """Dibuja una métrica (key in {'size','bpc','build'}) y la guarda.

    Retorna True si se generó la figura, False si no había datos.
    """
    iters = list(range(1, len(ns) + 1))
    any_series = False

    fig, ax = plt.subplots(figsize=(10, 5))

    # Línea de referencia: tamaño del texto original (solo en el plot de tamaño).
    if text_line:
        ax.plot(iters, ns, marker="o", linewidth=2, linestyle="--",
                color="steelblue", label="Texto original")
        any_series = True

    for idx_name, color in COLORS.items():
        series = by_index.get(idx_name)
        if not series:
            continue
        ys = [series.get(n, {}).get(key) for n in ns]
        iters_filt = [i for i, y in zip(iters, ys) if y is not None]
        ys_filt = [y for y in ys if y is not None]
        if not ys_filt:
            continue
        ax.plot(iters_filt, ys_filt, marker="s", linewidth=2,
                color=color, label=NAMES.get(idx_name, idx_name))
        any_series = True

    if not any_series:
        plt.close(fig)
        print(f"AVISO: sin datos para '{key}', se omite {out_path.name}",
              file=sys.stderr)
        return False

    ax.set_xticks(iters)
    ax.set_xticklabels([f"{i}\n({lbl})" for i, lbl in zip(iters, labels)],
                       fontsize=8)
    ax.set_xlabel("Iteración (tamaño del dataset)", fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, fontsize=13)

    if y_formatter is not None:
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(y_formatter))
    if yscale:
        ax.set_yscale(yscale)

    ax.legend(fontsize=10)
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    ax.grid(axis="x", linestyle=":", alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return True


def plot_grouped_bars(plt, ticker, *, ns, labels, by_index, key, title, ylabel,
                      y_formatter, value_fmt, out_path):
    """Barras agrupadas LZ77 vs sr-index en escala LINEAL.

    Pensado para comparar directamente magnitudes (no logaritmos) de unos pocos
    datasets. `ns`/`labels` deben venir ya recortados (p. ej. últimos 5).
    Retorna True si se generó la figura.
    """
    x = list(range(len(ns)))
    width = 0.38

    present = [(name, COLORS[name]) for name in COLORS
               if any(by_index.get(name, {}).get(n, {}).get(key) is not None
                      for n in ns)]
    if not present:
        print(f"AVISO: sin datos para '{key}' (barras), se omite {out_path.name}",
              file=sys.stderr)
        return False

    fig, ax = plt.subplots(figsize=(10, 5.5))
    n_series = len(present)
    # Centrar el grupo de barras alrededor de cada tick.
    offsets = [(-(n_series - 1) / 2 + i) * width for i in range(n_series)]

    for (idx_name, color), off in zip(present, offsets):
        series = by_index.get(idx_name, {})
        ys = [series.get(n, {}).get(key) for n in ns]
        xs = [xi + off for xi in x]
        bars = ax.bar(xs, [y if y is not None else 0 for y in ys], width,
                      color=color, label=NAMES.get(idx_name, idx_name))
        for rect, y in zip(bars, ys):
            if y is None:
                continue
            ax.annotate(value_fmt(y), xy=(rect.get_x() + rect.get_width() / 2,
                                          rect.get_height()),
                        xytext=(0, 3), textcoords="offset points",
                        ha="center", va="bottom", fontsize=7)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_xlabel("Dataset", fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, fontsize=13)
    if y_formatter is not None:
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(y_formatter))
    ax.legend(fontsize=10)
    ax.grid(axis="y", linestyle="--", alpha=0.4)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", help="CSV unificado de run_scaling_bench.sh")
    parser.add_argument("--out-dir", default=None,
                        help="directorio de salida (default: directorio del CSV)")
    parser.add_argument("--out", default=None,
                        help="alias retro-compat: nombre del PNG de tamaño")
    args = parser.parse_args()

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.ticker as ticker
    except ImportError:
        print("ERROR: matplotlib no disponible. Usa el venv del proyecto: "
              "../.venv/bin/python3 scripts/plot_scaling.py ...", file=sys.stderr)
        sys.exit(1)

    csv_path = Path(args.csv)
    rows = read_csv(csv_path)
    if not rows:
        print("ERROR: CSV vacío", file=sys.stderr)
        sys.exit(1)

    out_dir = Path(args.out_dir) if args.out_dir else csv_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    # by_index[index][n_bytes] = {"size": ..., "bpc": ..., "build": ...}
    by_index = defaultdict(dict)
    datasets_order = []  # (n_bytes, label) único por n_bytes
    seen = set()

    for r in rows:
        try:
            n = int(r["n_bytes"])
            idx = r["index"]
        except (KeyError, ValueError):
            continue
        label = r.get("dataset", "").replace("scaling_", "")
        by_index[idx][n] = {
            "size": _float_or_none(r.get("total_bytes")),
            "bpc": _float_or_none(r.get("total_bpc")),
            "build": _float_or_none(r.get("build_time_s")),
            "query": _float_or_none(r.get("locate_ext_us")),
        }
        if n not in seen:
            seen.add(n)
            datasets_order.append((n, label))

    if not datasets_order:
        print("ERROR: no se hallaron filas con columnas 'n_bytes' e 'index' "
              "(¿CSV en formato antiguo? re-corre run_scaling_bench.sh)",
              file=sys.stderr)
        sys.exit(1)

    datasets_order.sort(key=lambda x: x[0])
    ns = [x[0] for x in datasets_order]
    labels = [x[1] for x in datasets_order]

    # Plot de tamaño: respeta --out si se pasó (retro-compat), si no out_dir.
    size_out = Path(args.out) if args.out else out_dir / "scaling_size.png"

    specs = [
        dict(key="size", title="Texto original vs tamaño del índice (LZ77 y sr-index)",
             ylabel="Tamaño (bytes)", yscale="log", y_formatter=human_bytes,
             text_line=True, out_path=size_out),
        dict(key="bpc", title="Bits por carácter del índice (LZ77 y sr-index)",
             ylabel="bpc (bits/carácter)", yscale=None, y_formatter=None,
             text_line=False, out_path=out_dir / "scaling_bpc.png"),
        dict(key="build", title="Tiempo de construcción del índice (LZ77 y sr-index)",
             ylabel="Tiempo (s)", yscale="log", y_formatter=None,
             text_line=False, out_path=out_dir / "scaling_build.png"),
        dict(key="query", title="Tiempo de locate extremal (LZ77 y sr-index)",
             ylabel="Tiempo (μs/query)", yscale=None, y_formatter=None,
             text_line=False, out_path=out_dir / "scaling_query.png"),
    ]

    saved = []
    for spec in specs:
        if plot_metric(plt, ticker, ns=ns, labels=labels, by_index=by_index, **spec):
            saved.append(spec["out_path"])

    # Comparaciones directas (escala lineal) de los últimos 5 datasets.
    ns_last5 = ns[-5:]
    labels_last5 = labels[-5:]
    bar_specs = [
        dict(key="size", title="Tamaño del índice — últimos 5 datasets (escala lineal)",
             ylabel="Tamaño (bytes)", y_formatter=human_bytes,
             value_fmt=human_bytes, out_path=out_dir / "scaling_size_last5.png"),
        dict(key="query",
             title="Tiempo de locate extremal — últimos 5 datasets (escala lineal)",
             ylabel="Tiempo (μs/query)", y_formatter=None,
             value_fmt=lambda v: f"{v:.0f}",
             out_path=out_dir / "scaling_query_last5.png"),
    ]
    for spec in bar_specs:
        if plot_grouped_bars(plt, ticker, ns=ns_last5, labels=labels_last5,
                             by_index=by_index, **spec):
            saved.append(spec["out_path"])

    if not saved:
        print("ERROR: ninguna figura generada", file=sys.stderr)
        sys.exit(1)

    print("Guardado:")
    for p in saved:
        print(f"  {p}")


if __name__ == "__main__":
    main()
