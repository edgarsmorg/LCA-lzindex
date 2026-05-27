#!/usr/bin/env bash
# bench_sr_sweep.sh — sweep de s para SR-Index (3 variantes) vs LZ77-Index
#
# Uso:
#   bash scripts/bench_sr_sweep.sh <texto> <sr_index_dir> <data_name> <patterns> [min_s] [max_s]
#
# Ejemplos:
#   bash scripts/bench_sr_sweep.sh \
#       data/repetitive/reference.txt \
#       data/repetitive/sr_index_build \
#       reference.txt \
#       data/repetitive/queries.txt
#
#   bash scripts/bench_sr_sweep.sh \
#       data/Escherichia_Coli \
#       data/sr_index_build \
#       Escherichia_Coli \
#       data/patterns_1000.txt

set -e
REPO="$(cd "$(dirname "$0")/.." && pwd)"
BUILD="$REPO/build"
BM_LOCATE="$REPO/external/sr-index/build/benchmark/sr-index/bm_locate_ri"

TEXT="$1";    SR_DIR="$2";   DATA_NAME="$3"
PATTERNS="$4"; MIN_S="${5:-4}"; MAX_S="${6:-64}"

if [[ -z "$TEXT" || -z "$SR_DIR" || -z "$DATA_NAME" || -z "$PATTERNS" ]]; then
    echo "Uso: $0 <texto> <sr_index_dir> <data_name> <patterns> [min_s=4] [max_s=64]"
    exit 1
fi

N=$(wc -c < "$TEXT")
N_PAT=$(wc -l < "$PATTERNS")

echo "========================================================================"
echo " Benchmark: SR-Index sweep s=$MIN_S..$MAX_S  +  LZ77-Index"
echo " texto    : $TEXT  ($N bytes)"
echo " patrones : $PATTERNS  ($N_PAT patrones)"
echo "========================================================================"
echo ""

# ── LZ77-Index ────────────────────────────────────────────────────────────────
echo "Construyendo LZ77-Index (WmMinRmq)..."
LZ77_OUT=$("$BUILD/measure_index" "$TEXT" --patterns="$PATTERNS" --variant=wm 2>/dev/null)

LZ77_Z=$(echo     "$LZ77_OUT" | grep "^z (frases"      | grep -oP '\d+' | head -1)
LZ77_TOTAL=$(echo "$LZ77_OUT" | grep "TOTAL"            | grep -oP '\d+ (MB|KB)' | head -1)
LZ77_US=$(echo    "$LZ77_OUT" | grep "locate_extremal " | grep -oP '[0-9]+\.[0-9]+' | head -1)

LZ77_BYTES=$(echo "$LZ77_TOTAL" | awk '{if($2=="MB") print $1*1048576; else print $1*1024}')
LZ77_BPC=$(awk "BEGIN{printf \"%.4f\", $LZ77_BYTES*8/$N}")
LZ77_KB=$(awk  "BEGIN{printf \"%.1f\",  $LZ77_BYTES/1024}")
LZ77_MS=$(awk  "BEGIN{printf \"%.3f\",  $LZ77_US/1000}")

echo "LZ77-Index: z=$LZ77_Z frases, ${LZ77_KB} KB = ${LZ77_BPC} bpc, ${LZ77_US} µs/query"
echo ""

# ── SR-Index sweep ────────────────────────────────────────────────────────────
echo "Corriendo bm_locate_ri (s=$MIN_S..$MAX_S, reps=3)..."
BM_OUT=$("$BM_LOCATE" \
    --data_dir="$SR_DIR" \
    --data_name="$DATA_NAME" \
    --patterns="$PATTERNS" \
    --pattern_code=PLAIN \
    --min_s="$MIN_S" --max_s="$MAX_S" \
    --reps=3 2>&1)

echo ""
echo "========================================================================"
printf "%-22s  %8s  %8s  %12s  %12s\n" \
    "Índice" "KB" "bpc" "µs/query" "ns/occ"
printf "%s\n" "------------------------------------------------------------------------"

# ── Fila LZ77 ─────────────────────────────────────────────────────────────────
printf "%-22s  %8s  %8s  %12s  %12s\n" \
    "LZ77 (locate_ext)" "${LZ77_KB} KB" "$LZ77_BPC" "$LZ77_US" "—"

# ── Filas SR-Index (parsear output bm_locate_ri con awk) ─────────────────────
# Cada línea de datos tiene counters: Bits_x_Symbol=B Size(bytes)=Sk
#   Time_x_Pattern=Tms  Time_x_Occurrence=Xns
# "Size(bytes)" es el índice propio (no Collection_Size).
echo "$BM_OUT" | awk '
function to_kb(s,   v) {
    v = s
    if (v ~ /k$/) { sub(/k$/,"",v); return v * 1000 / 1024 }
    if (v ~ /M$/) { sub(/M$/,"",v); return v * 1000000 / 1024 }
    if (v ~ /G$/) { sub(/G$/,"",v); return v * 1000000000 / 1024 }
    return v / 1024
}
function to_us(s,   v) {
    v = s
    if (v ~ /ms$/) { sub(/ms$/,"",v); return v * 1000 }
    if (v ~ /µs$/ || v ~ /us$/) { sub(/[µu]s$/,"",v); return v }
    if (v ~ /ns$/) { sub(/ns$/,"",v); return v / 1000 }
    return v
}
function to_ns(s,   v) {
    v = s
    if (v ~ /ns$/) { sub(/ns$/,"",v); return v }
    if (v ~ /µs$/ || v ~ /us$/) { sub(/[µu]s$/,"",v); return v * 1000 }
    if (v ~ /ms$/) { sub(/ms$/,"",v); return v * 1000000 }
    return v
}
/^(SR-Index|R-Index)/ {
    variant = $1
    bpc = ""; size_kb = ""; txp_us = ""; txo_ns = ""
    for (i=1; i<=NF; i++) {
        if ($i ~ /^Bits_x_Symbol=/)         { split($i,a,"="); bpc    = sprintf("%.4f", a[2]) }
        # capturar solo " Size(bytes)=" (con espacio antes, no Collection_Size)
        if ($i ~ /^ *Size\(bytes\)=/ || ($i ~ /^Size\(bytes\)=/ && $(i-1) !~ /Collection/)) {
            split($i,a,"="); size_kb = sprintf("%.1f", to_kb(a[2]))
        }
        if ($i ~ /^Time_x_Pattern=/)        { split($i,a,"="); txp_us = sprintf("%.1f", to_us(a[2])) }
        if ($i ~ /^Time_x_Occurrence=/)     { split($i,a,"="); txo_ns = sprintf("%.1f", to_ns(a[2])) }
    }
    printf "%-22s  %8s  %8s  %12s  %12s\n", \
        variant, size_kb " KB", bpc, txp_us, txo_ns
}'

echo ""
echo "Nota: LZ77 locate_extremal ≠ SR locate_all."
echo "  SR mide tiempo de localizar TODAS las ocurrencias (depende de occ)."
echo "  LZ77 mide tiempo de obtener solo min/max posición (O(m^2 log^2 z))."
