#!/usr/bin/env bash
# run_compare_classifiers.sh — ejecuta compare_classifiers sobre data/repetitive
# para sr∈{16,32} y deja los CSV en results/compare/.
#
# Requisitos previos:
#   1. python3 scripts/tree_json_to_tsv.py         → tree.tsv + genomes.tsv
#   2. python3 scripts/gen_classification_reads.py  → reads.tsv
#   3. cmake --build build --target compare_classifiers
#
# Uso:
#   bash scripts/run_compare_classifiers.sh [min_mem=31]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJ_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

DATA_DIR="$PROJ_DIR/data/repetitive"
SR_BUILD_DIR="$DATA_DIR/sr_index_build"
BIN="$PROJ_DIR/build/compare_classifiers"
MIN_MEM="${1:-31}"

REF="$DATA_DIR/reference.txt"
TREE="$DATA_DIR/tree.tsv"
GENOMES="$DATA_DIR/genomes.tsv"
READS="$DATA_DIR/reads.tsv"

if [[ ! -f "$BIN" ]]; then
    echo "ERROR: $BIN no existe. Ejecutar:"
    echo "  cmake --build build --target compare_classifiers"
    exit 1
fi

for f in "$REF" "$TREE" "$GENOMES" "$READS"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: no existe $f"
        echo "Ejecutar los scripts de preparación de datos primero."
        exit 1
    fi
done

mkdir -p "$PROJ_DIR/results/compare"

for SR in 16 32; do
    echo "============================================================"
    echo "  sr = $SR"
    echo "============================================================"
    OUT="$PROJ_DIR/results/compare/compare_repetitive_sr${SR}.csv"
    "$BIN" "$REF" "$TREE" "$GENOMES" "$READS" "$SR_BUILD_DIR" \
           "$SR" "--out=$OUT" "--min-mem=$MIN_MEM"
    echo ""
done

echo "=== Resumen de resultados ==="
for SR in 16 32; do
    OUT="$PROJ_DIR/results/compare/compare_repetitive_sr${SR}.csv"
    if [[ -f "$OUT" ]]; then
        echo "sr=$SR  — $(wc -l < "$OUT") líneas en $OUT"
    fi
done
