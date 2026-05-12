#!/usr/bin/env bash
# bench_vs_sr.sh — compara LZ77-Index (WmMinRmq) vs sr-index (s=16) en E. coli
#
# Uso:
#   bash scripts/bench_vs_sr.sh [--skip-lz77]
#
# --skip-lz77: omite la construcción del LZ77-index (si ya se corrió antes y
#              quieres solo ver el tamaño del sr-index + su velocidad).
#
# Para comparar WtMinRmq vs WmMinRmq en espacio/tiempo, correr:
#   measure_index <texto> --patterns=<patrones> --variant=both

set -e
REPO="$(cd "$(dirname "$0")/.." && pwd)"
BUILD="$REPO/build"
DATA="$REPO/data"
SR_DIR="$REPO/results/baseline_ecoli"
PATTERNS="$DATA/patterns_ecoli_50.txt"
TEXT="$DATA/Escherichia_Coli"
N=112689515   # tamaño de E. coli (bytes)

# ── Verificaciones ────────────────────────────────────────────────────────────
for f in "$TEXT" "$PATTERNS" "$BUILD/measure_index" "$BUILD/sr_locate"; do
    [[ -f "$f" ]] || { echo "ERROR: falta $f"; exit 1; }
done
[[ -d "$SR_DIR" ]] || { echo "ERROR: falta $SR_DIR (correr scripts/run_sr_index.sh primero)"; exit 1; }

N_PAT=$(wc -l < "$PATTERNS")
echo "================================================================"
echo " Benchmark: LZ77-Index vs SR-Index (s=16)"
echo " texto  : $TEXT"
echo " patrones: $PATTERNS  ($N_PAT × ~50bp)"
echo "================================================================"

# ── LZ77-Index ────────────────────────────────────────────────────────────────
if [[ "$1" != "--skip-lz77" ]]; then
    echo ""
    echo ">>> LZ77-Index: construcción + bench locate (WmMinRmq)"
    echo "    (puede tardar ~8-15 min en 108 MB)"
    /usr/bin/time -v "$BUILD/measure_index" "$TEXT" --patterns="$PATTERNS" --variant=wm \
        2> "$SR_DIR/lz77_time.txt"
    echo ""
    echo "--- Pico de memoria LZ77 construction ---"
    grep "Maximum resident" "$SR_DIR/lz77_time.txt" || true
fi

# ── SR-Index (s=16): tamaño en disco ─────────────────────────────────────────
echo ""
echo ">>> SR-Index s=16: tamaño en disco"
SR_SIZE=$(du -sb \
    "$SR_DIR"/16_bwt_run_first_text_pos_by_last_Escherichia_Coli.txt.sdsl \
    "$SR_DIR"/16_bwt_run_first_text_pos_sorted_to_last_idx_Escherichia_Coli.txt.sdsl \
    "$SR_DIR"/16_bwt_run_first_text_pos_sorted_valid_area_Escherichia_Coli.txt.sdsl \
    "$SR_DIR"/16_bwt_run_first_text_pos_sorted_valid_mark_Escherichia_Coli.txt.sdsl \
    "$SR_DIR"/16_bwt_run_last_idx_Escherichia_Coli.txt.sdsl \
    "$SR_DIR"/16_bwt_run_last_text_pos_Escherichia_Coli.txt.sdsl \
    "$SR_DIR"/bwt_rle_Escherichia_Coli.txt.sdsl \
    2>/dev/null | awk '{s+=$1} END {print s}')
SR_MB=$(echo "scale=1; $SR_SIZE / 1048576" | bc)
SR_BPC=$(echo "scale=3; $SR_SIZE * 8 / $N" | bc)
echo "  SR-Index (s=16): ${SR_MB} MB  (${SR_BPC} bpc)"

# ── SR-Index: velocidad locate ────────────────────────────────────────────────
echo ""
echo ">>> SR-Index s=16: bench locate ($N_PAT patrones)"
T0=$(date +%s%N)
"$BUILD/sr_locate" \
    --data_dir="$SR_DIR" \
    --data_name="Escherichia_Coli.txt" \
    --patterns="$PATTERNS" > /dev/null
T1=$(date +%s%N)
TOTAL_MS=$(( (T1 - T0) / 1000000 ))
TOTAL_US=$(( (T1 - T0) / 1000 ))
US_PER_QUERY=$(echo "scale=1; $TOTAL_US / $N_PAT" | bc)
echo "  ${N_PAT} patrones: ${TOTAL_MS} ms total"
echo "  Velocidad: ${US_PER_QUERY} µs/query"

echo ""
echo "================================================================"
echo " RESUMEN"
echo "================================================================"
echo "  SR-Index (s=16):   ${SR_MB} MB  = ${SR_BPC} bpc"
echo "  SR locate:         ${US_PER_QUERY} µs/query (full enumerate)"
echo "  (ver salida de measure_index para métricas del LZ77-Index)"
