#!/usr/bin/env bash
set -euo pipefail

manifest="bench/datasets.json"
data_out="data/bench"
out=""
s_list="4,16,64"
rebuild=0
build_dir="build"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest) manifest="$2"; shift 2 ;;
    --out) out="$2"; shift 2 ;;
    --data-out) data_out="$2"; shift 2 ;;
    --s-list) s_list="$2"; shift 2 ;;
    --build-dir) build_dir="$2"; shift 2 ;;
    --rebuild) rebuild=1; shift ;;
    *) echo "unknown argument: $1" >&2; exit 1 ;;
  esac
done

if [[ -z "$out" ]]; then
  out="results/bench/$(date -u +%Y%m%dT%H%M%SZ)"
fi

required=(
  "$build_dir/build_sr_index"
  "$build_dir/inspect_sr_index"
  "$build_dir/build_lz_index"
  "$build_dir/inspect_lz_index"
  "$build_dir/bench_locate_compare"
)
for bin in "${required[@]}"; do
  if [[ ! -x "$bin" ]]; then
    echo "missing binary: $bin" >&2
    echo "build with: cmake --build $build_dir --target build_sr_index inspect_sr_index build_lz_index inspect_lz_index bench_locate_compare" >&2
    exit 2
  fi
done

mkdir -p "$out/logs" "$out/indexes" "$data_out"

python3 scripts/gen_bench_datasets.py --manifest "$manifest" --out "$data_out" > "$out/logs/generate.log"

python3 - "$manifest" "$out/run_config.json" "$data_out" "$s_list" <<'PY'
import json
import sys
from pathlib import Path

manifest = json.loads(Path(sys.argv[1]).read_text())
config = {
    "manifest": sys.argv[1],
    "data_out": sys.argv[3],
    "s_list": [int(x) for x in sys.argv[4].split(",") if x],
    "datasets": manifest.get("datasets", []),
}
Path(sys.argv[2]).write_text(json.dumps(config, indent=2) + "\n")
PY

mapfile -t datasets < <(python3 - "$manifest" <<'PY'
import json
import sys
from pathlib import Path
manifest = json.loads(Path(sys.argv[1]).read_text())
for spec in manifest.get("datasets", []):
    print(spec["name"])
PY
)

IFS=',' read -r -a s_values <<< "$s_list"

for dataset in "${datasets[@]}"; do
  text="$data_out/$dataset/reference.txt"
  patterns="$data_out/$dataset/patterns.txt"
  lz_prefix="$out/indexes/$dataset/lz77"
  mkdir -p "$out/indexes/$dataset"

  if [[ $rebuild -eq 1 || ! -f "${lz_prefix}.meta" ]]; then
    "$build_dir/build_lz_index" "$text" "$lz_prefix" --name="$dataset" --csv="$out/build.csv" > "$out/logs/${dataset}_build_lz.log" 2>&1
  fi
  "$build_dir/inspect_lz_index" "$text" "$lz_prefix" --name="$dataset" --csv="$out/size.csv" > "$out/logs/${dataset}_inspect_lz.log" 2>&1

  for s in "${s_values[@]}"; do
    sr_dir="$out/indexes/$dataset/sr_s$s"
    mkdir -p "$sr_dir"
    if [[ $rebuild -eq 1 || ! -f "$sr_dir/alphabet_${dataset}.sdsl" ]]; then
      "$build_dir/build_sr_index" "$text" "$sr_dir" --s="$s" --name="$dataset" --csv="$out/build.csv" > "$out/logs/${dataset}_build_sr_s${s}.log" 2>&1
    fi
    "$build_dir/inspect_sr_index" "$text" "$sr_dir" --s="$s" --name="$dataset" --csv="$out/size.csv" > "$out/logs/${dataset}_inspect_sr_s${s}.log" 2>&1
    "$build_dir/bench_locate_compare" "$text" "$patterns" "$lz_prefix" "$sr_dir" --s="$s" --name="$dataset" --reps=3 --csv="$out/locate.csv" > "$out/logs/${dataset}_locate_s${s}.log" 2>&1
  done
done

echo "results: $out"
echo "csv: $out/build.csv $out/size.csv $out/locate.csv"
