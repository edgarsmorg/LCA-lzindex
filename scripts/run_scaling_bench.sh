#!/usr/bin/env bash
# run_scaling_bench.sh — escalado LZ77 + sr-index: espacio y tiempo de query.
#
# Pipeline por dataset:
#   build_lz_index    → índice LZ77 serializado
#   inspect_lz_index  → espacio LZ77
#   build_sr_index    → índice SR serializado
#   inspect_sr_index  → espacio SR
#   bench_locate_compare → tiempo de query (ambos)
#
# Salida:
#   $out/raw/build.csv   — tiempos de construcción (LZ + SR)
#   $out/raw/size.csv    — espacio por componente  (LZ + SR)
#   $out/raw/locate.csv  — tiempo de query: locate_extremal + paso LCA (LZ + SR)
#   $out/scaling.csv     — CSV unificado (join de los 3)
#   $out/logs/           — stdout/stderr por dataset
#
# Uso:
#   scripts/run_scaling_bench.sh [--manifest=bench/datasets_scaling.json]
#                                [--build-dir=build] [--out=results/scaling/RUN]
#                                [--data-out=data/scaling] [--s=16] [--skip-gen] [--rebuild]
set -euo pipefail

manifest="bench/datasets_scaling.json"
data_out="data/scaling"
build_dir="build"
out=""
s_val=16
skip_gen=0
rebuild=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest=*)  manifest="${1#*=}";  shift ;;
    --build-dir=*) build_dir="${1#*=}"; shift ;;
    --out=*)       out="${1#*=}";       shift ;;
    --data-out=*)  data_out="${1#*=}";  shift ;;
    --s=*)         s_val="${1#*=}";     shift ;;
    --skip-gen)    skip_gen=1;          shift ;;
    --rebuild)     rebuild=1;           shift ;;
    *) echo "argumento desconocido: $1" >&2; exit 1 ;;
  esac
done

if [[ -z "$out" ]]; then
  out="results/scaling/$(date -u +%Y%m%dT%H%M%SZ)"
fi

# ── Verificar binarios requeridos ────────────────────────────────────────────
required_lz=(
  "$build_dir/build_lz_index"
  "$build_dir/inspect_lz_index"
)
required_sr=(
  "$build_dir/build_sr_index"
  "$build_dir/inspect_sr_index"
  "$build_dir/bench_locate_compare"
)

for bin in "${required_lz[@]}"; do
  if [[ ! -x "$bin" ]]; then
    echo "ERROR: falta binario LZ: $bin" >&2
    echo "  cmake --build $build_dir --target build_lz_index inspect_lz_index" >&2
    exit 2
  fi
done

sr_available=1
for bin in "${required_sr[@]}"; do
  if [[ ! -x "$bin" ]]; then
    echo "AVISO: binario SR no disponible ($bin) — se medirá solo LZ77" >&2
    sr_available=0
    break
  fi
done

mkdir -p "$out/raw" "$out/logs" "$out/indexes" "$data_out"

build_csv="$out/raw/build.csv"
size_csv="$out/raw/size.csv"
locate_csv="$out/raw/locate.csv"

# ── 1. Generar datasets ──────────────────────────────────────────────────────
if [[ $skip_gen -eq 0 ]]; then
  echo "Generando datasets desde $manifest..."
  python3 scripts/gen_bench_datasets.py \
    --manifest "$manifest" \
    --out "$data_out" \
    2>&1 | tee "$out/logs/generate.log"
else
  echo "Saltando generación (--skip-gen)"
fi

# ── 2. Extraer nombres de datasets ───────────────────────────────────────────
mapfile -t datasets < <(python3 - "$manifest" <<'PY'
import json, sys
from pathlib import Path
m = json.loads(Path(sys.argv[1]).read_text())
for spec in m.get("datasets", []):
    print(spec["name"])
PY
)

echo ""
echo "Datasets: ${#datasets[@]}"
echo "SR disponible: $sr_available  (s=$s_val)"
echo ""

# ── 3. Pipeline por dataset ──────────────────────────────────────────────────
for dataset in "${datasets[@]}"; do
  ref="$data_out/$dataset/reference.txt"
  pat="$data_out/$dataset/patterns.txt"
  lz_prefix="$out/indexes/$dataset/lz77"
  sr_dir="$out/indexes/$dataset/sr_s${s_val}"

  if [[ ! -f "$ref" ]]; then
    echo "AVISO: no existe $ref, saltando $dataset" >&2
    continue
  fi

  mkdir -p "$out/indexes/$dataset"
  echo "── $dataset ──────────────────────────────"

  # build LZ
  if [[ $rebuild -eq 1 || ! -f "${lz_prefix}.meta" ]]; then
    echo -n "  build_lz_index  ... "
    "$build_dir/build_lz_index" "$ref" "$lz_prefix" \
      --name="$dataset" --csv="$build_csv" \
      > "$out/logs/${dataset}_build_lz.log" 2>&1
    echo "OK"
  else
    echo "  build_lz_index  ... (ya existe, saltando)"
  fi

  # inspect LZ (espacio)
  echo -n "  inspect_lz_index... "
  "$build_dir/inspect_lz_index" "$ref" "$lz_prefix" \
    --name="$dataset" --csv="$size_csv" \
    > "$out/logs/${dataset}_inspect_lz.log" 2>&1
  echo "OK"

  if [[ $sr_available -eq 1 ]]; then
    # build SR
    if [[ $rebuild -eq 1 || ! -d "$sr_dir" ]] || \
       [[ -z "$(ls -A "$sr_dir" 2>/dev/null)" ]]; then
      mkdir -p "$sr_dir"
      echo -n "  build_sr_index  ... "
      "$build_dir/build_sr_index" "$ref" "$sr_dir" \
        --s="$s_val" --name="$dataset" --csv="$build_csv" \
        > "$out/logs/${dataset}_build_sr.log" 2>&1
      echo "OK"
    else
      echo "  build_sr_index  ... (ya existe, saltando)"
    fi

    # inspect SR (espacio)
    echo -n "  inspect_sr_index... "
    "$build_dir/inspect_sr_index" "$ref" "$sr_dir" \
      --s="$s_val" --name="$dataset" --csv="$size_csv" \
      > "$out/logs/${dataset}_inspect_sr.log" 2>&1
    echo "OK"

    # bench locate (ambos índices, 3 reps)
    echo -n "  bench_locate    ... "
    "$build_dir/bench_locate_compare" "$ref" "$pat" "$lz_prefix" "$sr_dir" \
      --s="$s_val" --name="$dataset" --reps=3 --csv="$locate_csv" \
      > "$out/logs/${dataset}_locate.log" 2>&1
    echo "OK"
  fi

  echo ""
done

# ── 4. Merge → scaling.csv unificado ────────────────────────────────────────
echo "Generando $out/scaling.csv..."

python3 - "$build_csv" "$size_csv" "$locate_csv" "$out/scaling.csv" "$sr_available" <<'PY'
import csv, sys
from collections import defaultdict
from pathlib import Path
from statistics import mean

build_path, size_path, locate_path, out_path, sr_flag = sys.argv[1:]
sr_available = int(sr_flag)

def read(path):
    p = Path(path)
    if not p.exists():
        return []
    return list(csv.DictReader(p.open(newline="")))

# Índice por (dataset, index, s)
build_map  = {}
size_map   = {}
locate_map = defaultdict(lambda: defaultdict(list))  # key -> metric -> [valores]
npat_map   = {}

for r in read(build_path):
    key = (r["dataset"], r["index"], r.get("s", "0"))
    build_map[key] = r

for r in read(size_path):
    key = (r["dataset"], r["index"], r.get("s", "0"))
    size_map[key] = r

def _f(v):
    try: return float(v)
    except (TypeError, ValueError): return None

for r in read(locate_path):
    key = (r["dataset"], r["index"], r.get("s", "0"))
    for m in ("us_per_query", "lca_us", "total_us"):
        v = _f(r.get(m))
        if v is not None:
            locate_map[key][m].append(v)
    npat_map[key] = r.get("n_patterns", "")

header = ("dataset,n_bytes,index,s,build_time_s,total_bytes,total_bpc,"
          "n_patterns,locate_ext_us,lca_us,locate_lca_us")
rows = []

def _avg(key, metric):
    vals = locate_map.get(key, {}).get(metric, [])
    return f"{mean(vals):.6f}" if vals else ""

for key, sz in size_map.items():
    dataset, index, s = key
    n_bytes    = sz.get("n_bytes", "")
    total_b    = sz.get("serialized_bytes", "")
    bpc        = sz.get("bpc", "")
    build_s    = build_map.get(key, {}).get("build_seconds", "")
    n_pat      = npat_map.get(key, "")
    loc_us     = _avg(key, "us_per_query")   # tiempo de locate_extremal
    lca_us     = _avg(key, "lca_us")         # tiempo del paso LCA (pos->genoma->lca)
    total_us   = _avg(key, "total_us")       # locate + lca

    rows.append((dataset, n_bytes, index, s, build_s,
                 total_b, bpc, n_pat, loc_us, lca_us, total_us))

# Ordenar: primero por n_bytes, luego por index
try:
    rows.sort(key=lambda r: (int(r[1]) if r[1] else 0, r[2]))
except Exception:
    pass

with open(out_path, "w", newline="") as f:
    f.write(header + "\n")
    for r in rows:
        f.write(",".join(str(x) for x in r) + "\n")

print(f"  filas escritas: {len(rows)}")
PY

echo ""
echo "=== Resultados ==="
echo "  Logs    : $out/logs/"
echo "  Raw CSV : $out/raw/"
echo "  Unified : $out/scaling.csv"
if [[ -f "$out/scaling.csv" ]]; then
  echo ""
  cat "$out/scaling.csv"
fi
