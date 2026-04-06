#!/bin/bash

# Script para ejecutar benchmarks del sr-index
# Uso: ./scripts/run_sr_index.sh <data_file> [min_s] [max_s]

set -e

SR_BUILD="./external/sr-index/build/benchmark/sr-index"
DATA="${1:-data/test_plain.txt}"
MIN_S="${2:-4}"
MAX_S="${3:-64}"
RESULTS_DIR="results"

if [ ! -f "$DATA" ]; then
    echo "Error: archivo de datos no encontrado: $DATA"
    exit 1
fi

mkdir -p "$RESULTS_DIR"
BASENAME=$(basename "$DATA" .txt)

echo "=== sr-index benchmark ==="
echo "Datos: $DATA ($(wc -c < "$DATA") bytes)"
echo "Parámetro s: $MIN_S .. $MAX_S"
echo ""

# Benchmark de construcción
echo "--- Construcción ---"
"$SR_BUILD/bm_construct_ri" \
    -data="$DATA" \
    -min_s="$MIN_S" \
    -max_s="$MAX_S" \
    --benchmark_format=json \
    --benchmark_out="$RESULTS_DIR/${BASENAME}_construct.json" \
    2>&1 | grep -v "^$"

echo ""
echo "✓ Resultados guardados en $RESULTS_DIR/${BASENAME}_construct.json"
