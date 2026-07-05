#!/usr/bin/env bash
# run_classify_bench.sh — benchmark de la pipeline de clasificación CON MEMs.
#
# Por dataset (mutated_copies con tree/genomes/reads generados):
#   build_lz_index           → índice LZ77 serializado
#   build_sr_index           → índice SR serializado
#   bench_classify_compare   → read → MEMs → locate → LCA (LZ vs SR)
#
# Salida:
#   $out/classify.csv  — una fila por (dataset, index, rep)
#   $out/logs/         — stdout/stderr por dataset
#
# Uso:
#   scripts/run_classify_bench.sh [--manifest=bench/datasets_classify.json]
#                                 [--build-dir=build] [--out=results/classify/RUN]
#                                 [--data-out=data/classify] [--s=16] [--reps=3]
#                                 [--min-mem=31] [--skip-gen] [--rebuild]
set -euo pipefail

manifest="bench/datasets_classify.json"
data_out="data/classify"
build_dir="build"
out=""
s_val=16
reps=3
min_mem=31
skip_gen=0
rebuild=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest=*)  manifest="${1#*=}";  shift ;;
    --build-dir=*) build_dir="${1#*=}"; shift ;;
    --out=*)       out="${1#*=}";       shift ;;
    --data-out=*)  data_out="${1#*=}";  shift ;;
    --s=*)         s_val="${1#*=}";     shift ;;
    --reps=*)      reps="${1#*=}";      shift ;;
    --min-mem=*)   min_mem="${1#*=}";   shift ;;
    --skip-gen)    skip_gen=1;          shift ;;
    --rebuild)     rebuild=1;           shift ;;
    *) echo "argumento desconocido: $1" >&2; exit 1 ;;
  esac
done

if [[ -z "$out" ]]; then
  out="results/classify/$(date -u +%Y%m%dT%H%M%SZ)"
fi

required=(
  "$build_dir/build_lz_index"
  "$build_dir/build_sr_index"
  "$build_dir/bench_classify_compare"
)
for bin in "${required[@]}"; do
  if [[ ! -x "$bin" ]]; then
    echo "ERROR: falta binario: $bin" >&2
    echo "  cmake --build $build_dir --target build_lz_index build_sr_index bench_classify_compare" >&2
    exit 2
  fi
done

mkdir -p "$out/logs" "$out/indexes" "$data_out"
csv="$out/classify.csv"

# ── 1. Generar datasets (ref + tree + genomes + reads) ───────────────────────
if [[ $skip_gen -eq 0 ]]; then
  echo "Generando datasets desde $manifest..."
  python3 scripts/gen_bench_datasets.py --manifest "$manifest" --out "$data_out" \
    2>&1 | tee "$out/logs/generate.log"
else
  echo "Saltando generación (--skip-gen)"
fi

mapfile -t datasets < <(python3 - "$manifest" <<'PY'
import json, sys
from pathlib import Path
m = json.loads(Path(sys.argv[1]).read_text())
for spec in m.get("datasets", []):
    print(spec["name"])
PY
)

echo ""
echo "Datasets: ${#datasets[@]}  (s=$s_val, reps=$reps, min_mem=$min_mem)"
echo ""

# ── 2. Pipeline por dataset ──────────────────────────────────────────────────
for dataset in "${datasets[@]}"; do
  dir="$data_out/$dataset"
  ref="$dir/reference.txt"
  tree="$dir/tree.tsv"
  genomes="$dir/genomes.tsv"
  reads="$dir/reads.tsv"
  lz_prefix="$out/indexes/$dataset/lz77"
  sr_dir="$out/indexes/$dataset/sr_s${s_val}"

  for req in "$ref" "$tree" "$genomes" "$reads"; do
    if [[ ! -f "$req" ]]; then
      echo "AVISO: falta $req, saltando $dataset" >&2
      continue 2
    fi
  done

  mkdir -p "$out/indexes/$dataset"
  echo "── $dataset ──────────────────────────────"

  if [[ $rebuild -eq 1 || ! -f "${lz_prefix}.meta" ]]; then
    echo -n "  build_lz_index ... "
    "$build_dir/build_lz_index" "$ref" "$lz_prefix" --name="$dataset" \
      > "$out/logs/${dataset}_build_lz.log" 2>&1
    echo "OK"
  else
    echo "  build_lz_index ... (ya existe)"
  fi

  if [[ $rebuild -eq 1 || ! -d "$sr_dir" || -z "$(ls -A "$sr_dir" 2>/dev/null)" ]]; then
    mkdir -p "$sr_dir"
    echo -n "  build_sr_index ... "
    "$build_dir/build_sr_index" "$ref" "$sr_dir" --s="$s_val" --name="$dataset" \
      > "$out/logs/${dataset}_build_sr.log" 2>&1
    echo "OK"
  else
    echo "  build_sr_index ... (ya existe)"
  fi

  echo -n "  bench_classify ... "
  "$build_dir/bench_classify_compare" "$ref" "$tree" "$genomes" "$reads" \
    "$lz_prefix" "$sr_dir" --s="$s_val" --name="$dataset" --reps="$reps" \
    --min-mem="$min_mem" --csv="$csv" \
    > "$out/logs/${dataset}_classify.log" 2>&1
  echo "OK"
  echo ""
done

echo "=== Resultados ==="
echo "  CSV  : $csv"
echo "  Logs : $out/logs/"
if [[ -f "$csv" ]]; then echo ""; cat "$csv"; fi
