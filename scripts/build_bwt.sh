#!/bin/bash

# Script para construir BWT usando ropebwt3

set -e

ROPEBWT3="${1:-./external/ropebwt3/ropebwt3}"
INPUT="${2:-data/test_seq.fa}"
OUTPUT="${3:-data/index.fmd}"
THREADS="${4:-$(nproc)}"

if [ ! -f "$ROPEBWT3" ]; then
    echo "Error: ropebwt3 no encontrado en $ROPEBWT3"
    exit 1
fi

if [ ! -f "$INPUT" ]; then
    echo "Error: archivo de entrada no encontrado: $INPUT"
    exit 1
fi

echo "Construyendo BWT para: $INPUT"
echo "Salida: $OUTPUT"
echo "Threads: $THREADS"

# Construir índice dinámico (.fmr)
TEMP_INDEX="${OUTPUT%.fmd}.fmr"
echo "Paso 1: Construyendo índice dinámico..."
"$ROPEBWT3" build -t"$THREADS" -bo "$TEMP_INDEX" "$INPUT"

# Convertir a formato estático (.fmd)
echo "Paso 2: Convirtiendo a formato estático..."
"$ROPEBWT3" build -i "$TEMP_INDEX" -do "$OUTPUT"

# Limpiar temporal
rm -f "$TEMP_INDEX"

echo "✓ BWT construida: $(ls -lh $OUTPUT | awk '{print $5, $9}')"
