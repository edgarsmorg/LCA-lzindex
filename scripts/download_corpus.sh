#!/bin/bash

# Script para descargar Pizza&Chili repetitive corpus
# https://pizzachili.dcc.uchile.cl/repcorpus.html

set -e

DATA_DIR="${1:-.}/data"
mkdir -p "$DATA_DIR"

echo "Descargando corpus del Pizza&Chili Repetitive..."

# Escherichia_coli (~108 MB)
echo "Descargando Escherichia_coli..."
wget -O "$DATA_DIR/Escherichia_coli.gz" \
  "https://pizzachili.dcc.uchile.cl/repcorpus/real/Escherichia_coli.gz" 2>&1 | tail -5

echo "Descomprimiendo Escherichia_coli..."
gunzip -f "$DATA_DIR/Escherichia_coli.gz"

echo "✓ Corpus descargado en $DATA_DIR"
ls -lh "$DATA_DIR"/Escherichia_coli
