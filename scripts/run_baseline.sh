#!/bin/bash
# ============================================================
# run_baseline.sh — Pipeline baseline completo usando sr-index
#
# Dados un texto de referencia (genomas concatenados en orden
# DFS) y lecturas query en FASTA, este script:
#   1. Construye el índice BWT con ropebwt3
#   2. Construye el sr-index con bm_construct_ri
#   3. Extrae MEMs de las lecturas con ropebwt3
#   4. Busca las posiciones de cada MEM con sr_locate
#
# Uso:
#   ./scripts/run_baseline.sh <reference.fa> <queries.fa> [outdir] [min_mem_len] [sr_s]
#
# Ejemplo:
#   ./scripts/run_baseline.sh data/synthetic/reference.fa data/test_seq.fa results/baseline_run
#
# Requisitos (compilados en external/ y build/):
#   - external/ropebwt3/ropebwt3
#   - external/sr-index/build/benchmark/sr-index/bm_construct_ri
#   - build/sr_locate
# ============================================================

set -euo pipefail

# ── Argumentos ─────────────────────────────────────────────
REFERENCE="${1:-}"
QUERIES="${2:-}"
OUTDIR="${3:-results/baseline_run}"
MIN_MEM_LEN="${4:-19}"   # mínimo largo de MEM (ropebwt3 default: 19)
SR_S="${5:-16}"          # parámetro de muestreo del sr-index (potencia de 2)

if [[ -z "$REFERENCE" || -z "$QUERIES" ]]; then
    echo "Uso: $0 <reference.fa> <queries.fa> [outdir] [min_mem_len] [sr_s]"
    echo ""
    echo "  reference.fa   Genomas concatenados en orden DFS (FASTA)"
    echo "  queries.fa     Lecturas a clasificar (FASTA)"
    echo "  outdir         Directorio de salida (default: results/baseline_run)"
    echo "  min_mem_len    Largo mínimo de MEM (default: 19)"
    echo "  sr_s           Parámetro de muestreo sr-index (default: 16)"
    exit 1
fi

# ── Rutas de binarios ───────────────────────────────────────
ROPEBWT3="./external/ropebwt3/ropebwt3"
BM_CONSTRUCT="./external/sr-index/build/benchmark/sr-index/bm_construct_ri"
SR_LOCATE="./build/sr_locate"

for BIN in "$ROPEBWT3" "$BM_CONSTRUCT" "$SR_LOCATE"; do
    if [[ ! -x "$BIN" ]]; then
        echo "ERROR: binario no encontrado o no ejecutable: $BIN"
        echo "Compila el proyecto primero: cmake -S . -B build && cmake --build build"
        exit 1
    fi
done

if [[ ! -f "$REFERENCE" ]]; then
    echo "ERROR: archivo de referencia no encontrado: $REFERENCE"
    exit 1
fi

if [[ ! -f "$QUERIES" ]]; then
    echo "ERROR: archivo de queries no encontrado: $QUERIES"
    exit 1
fi

# ── Setup ───────────────────────────────────────────────────
mkdir -p "$OUTDIR"
BASENAME=$(basename "$REFERENCE" .fa)
REF_TXT="$OUTDIR/${BASENAME}.txt"
BWT_FMR="$OUTDIR/${BASENAME}.fmr"
BWT_FMD="$OUTDIR/${BASENAME}.fmd"
MEMS_FILE="$OUTDIR/mems.bed"
PATTERNS_FILE="$OUTDIR/patterns.txt"
LOCATIONS_FILE="$OUTDIR/locations.txt"
THREADS=$(nproc)

echo "============================================================"
echo "  Baseline sr-index — pipeline completo"
echo "============================================================"
echo "  Referencia:    $REFERENCE"
echo "  Queries:       $QUERIES"
echo "  Salida:        $OUTDIR"
echo "  Min MEM len:   $MIN_MEM_LEN"
echo "  sr param s:    $SR_S"
echo "  Threads:       $THREADS"
echo "============================================================"
echo ""

# ── Paso 1: Convertir FASTA a texto plano ──────────────────
echo "[1/5] Convirtiendo referencia a texto plano..."
# Elimina cabeceras FASTA y newlines internos; separa secuencias con '#'
python3 - <<PYEOF
import re, sys

fasta = open("$REFERENCE").read()
seqs = re.split(r'>.*\n', fasta)[1:]  # descarta cabecera inicial vacía
# Unir cada secuencia (quitar newlines internos)
clean = [''.join(s.split()) for s in seqs]
# Concatenar con separador '#' y terminar en '#'
result = '#'.join(clean) + '#'
with open("$REF_TXT", 'w') as f:
    f.write(result)
print(f"  Texto plano: {len(result)} chars, {len(clean)} secuencias")
PYEOF

# ── Paso 2: Construir BWT con ropebwt3 ─────────────────────
if [[ -f "$BWT_FMD" ]]; then
    echo "[2/5] BWT ya existe ($BWT_FMD) — omitiendo construcción."
else
    echo "[2/5] Construyendo BWT con ropebwt3..."
    "$ROPEBWT3" build -t"$THREADS" -bo "$BWT_FMR" "$REFERENCE"
    "$ROPEBWT3" build -i "$BWT_FMR" -do "$BWT_FMD"
    rm -f "$BWT_FMR"
    echo "  BWT construida: $(du -h "$BWT_FMD" | cut -f1)"
fi

# ── Paso 3: Construir sr-index ─────────────────────────────
# bm_construct_ri guarda los archivos .sdsl en el directorio actual (cwd)
# → lo ejecutamos desde OUTDIR para no contaminar el directorio raíz
SDSL_CHECK="$OUTDIR/${SR_S}_bwt_run_last_text_pos_${BASENAME}.txt.sdsl"
if [[ -f "$SDSL_CHECK" ]]; then
    echo "[3/5] sr-index ya construido (caché .sdsl existe) — omitiendo."
else
    echo "[3/5] Construyendo sr-index (s=$SR_S)..."
    (cd "$OUTDIR" && \
     "$OLDPWD/$BM_CONSTRUCT" \
         -data="$OLDPWD/$REF_TXT" \
         -min_s="$SR_S" \
         -max_s="$SR_S" \
         --benchmark_format=json \
         --benchmark_out="$OLDPWD/$OUTDIR/sr_construct.json" 2>&1 | grep -v "^$"
    )
    echo "  sr-index construido. Benchmark: $OUTDIR/sr_construct.json"
fi

# ── Paso 4: Extraer MEMs con ropebwt3 ──────────────────────
echo "[4/5] Extrayendo MEMs (min_len=$MIN_MEM_LEN)..."
"$ROPEBWT3" mem -t"$THREADS" -l"$MIN_MEM_LEN" \
    "$BWT_FMD" "$QUERIES" 2>/dev/null \
    > "$MEMS_FILE"

MEM_COUNT=$(wc -l < "$MEMS_FILE")
echo "  MEMs encontrados: $MEM_COUNT"

if [[ "$MEM_COUNT" -eq 0 ]]; then
    echo "  AVISO: no se encontraron MEMs. Prueba reducir min_mem_len."
    exit 0
fi

# Extraer las secuencias de los MEMs del FASTA de queries
# Formato BED: query_name  start  end  count
# → extraemos la subcadena query[start:end]
echo "  Extrayendo secuencias de MEMs desde $QUERIES..."
python3 - <<PYEOF
import re

# Cargar todas las secuencias del FASTA de queries
seqs = {}
current = None
for line in open("$QUERIES"):
    line = line.rstrip()
    if line.startswith('>'):
        current = line[1:].split()[0]
        seqs[current] = []
    elif current:
        seqs[current].append(line)
seqs = {k: ''.join(v) for k, v in seqs.items()}

# Extraer secuencias únicas de MEMs
patterns = set()
for line in open("$MEMS_FILE"):
    parts = line.rstrip().split('\t')
    if len(parts) < 3:
        continue
    name, start, end = parts[0], int(parts[1]), int(parts[2])
    if name in seqs:
        patterns.add(seqs[name][start:end])

with open("$PATTERNS_FILE", 'w') as f:
    for p in sorted(patterns):
        f.write(p + '\n')

print(f"  Patrones únicos extraídos: {len(patterns)}")
PYEOF

# ── Paso 5: Locate con sr_locate ───────────────────────────
echo "[5/5] Buscando posiciones con sr_locate..."
"$SR_LOCATE" \
    --data_dir="$OUTDIR" \
    --data_name="${BASENAME}.txt" \
    --patterns="$PATTERNS_FILE" \
    > "$LOCATIONS_FILE"

PATTERN_COUNT=$(grep -c "^>" "$LOCATIONS_FILE" || true)
TOTAL_HITS=$(grep -vc "^>" "$LOCATIONS_FILE" || true)
echo "  Patrones buscados:   $PATTERN_COUNT"
echo "  Ocurrencias totales: $TOTAL_HITS"

echo ""
echo "============================================================"
echo "  Resultados en: $OUTDIR/"
echo "    mems.bed       — MEMs crudos de ropebwt3"
echo "    patterns.txt   — secuencias únicas de MEMs"
echo "    locations.txt  — posiciones en la referencia"
echo "    sr_construct.json — benchmark de construcción"
echo "============================================================"
echo ""
echo "Próximo paso: mapear posiciones a genomas y computar LCA"
echo "  Ver: scripts/test_synthetic_lca.py para ejemplo completo"
