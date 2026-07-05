#!/bin/bash

# Descarga (opcional) un corpus del Pizza&Chili Repetitive Corpus.
# https://pizzachili.dcc.uchile.cl/repcorpus.html
#
# NOTA: la evaluación de la memoria NO usa corpus reales de Pizza&Chili, sino
# los datasets sintéticos de copias mutadas generados por gen_bench_datasets.py:
#   python3 scripts/gen_bench_datasets.py --manifest bench/datasets_scaling.json --out data/scaling
# Este script queda solo como utilidad para experimentos ad-hoc.
#
# Uso: scripts/download_corpus.sh [dest_root] [corpus_name]
#   dest_root    raíz donde crear data/ (default: .)
#   corpus_name  nombre en repcorpus/real (default: cere)

set -e

DATA_DIR="${1:-.}/data"
CORPUS="${2:-cere}"
mkdir -p "$DATA_DIR"

echo "Descargando '$CORPUS' del Pizza&Chili Repetitive Corpus..."
wget -O "$DATA_DIR/$CORPUS.gz" \
  "https://pizzachili.dcc.uchile.cl/repcorpus/real/$CORPUS.gz" 2>&1 | tail -5

echo "Descomprimiendo '$CORPUS'..."
gunzip -f "$DATA_DIR/$CORPUS.gz"

echo "✓ Corpus descargado en $DATA_DIR"
ls -lh "$DATA_DIR/$CORPUS"
