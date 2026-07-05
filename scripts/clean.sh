#!/bin/bash
# Limpia artefactos regenerables del proyecto (caché, datasets generados y
# resultados de benchmarks). NUNCA borra archivos versionados en git
# (metadata.json de datasets, notas .md de results/, fuentes).
#
# Uso: ./scripts/clean.sh [--dry-run] [--data] [--results] [--build] [--all]
#   (sin flags)  Solo caché suelta de sr-index/sdsl y temporales.
#   --data       Datasets sintéticos generados (data/scaling, classify, sweep).
#   --results    Salidas de benchmarks regenerables (results/{bench,classify,scaling}).
#   --build      Directorio build/ del proyecto.
#   --all        Equivale a --data --results --build.
#   --dry-run    Solo muestra qué se borraría, sin borrar nada.

set -e

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DRY_RUN=0
CLEAN_DATA=0
CLEAN_RESULTS=0
CLEAN_BUILD=0

for arg in "$@"; do
  case "$arg" in
    --dry-run) DRY_RUN=1 ;;
    --data)    CLEAN_DATA=1 ;;
    --results) CLEAN_RESULTS=1 ;;
    --build)   CLEAN_BUILD=1 ;;
    --all)     CLEAN_DATA=1; CLEAN_RESULTS=1; CLEAN_BUILD=1 ;;
    *) echo "Opción desconocida: $arg" >&2; exit 1 ;;
  esac
done

delete() {
  if [ $DRY_RUN -eq 1 ]; then
    echo "  [dry-run] rm -rf $1"
  else
    rm -rf "$1"
  fi
}

# Borra por patrón bajo un directorio, contando coincidencias.
delete_glob() {  # $1 = dir base, resto = args de find
  local base="$1"; shift
  local count=0
  [ -d "$base" ] || { echo "  0 (no existe $base)"; return; }
  while IFS= read -r -d '' f; do
    delete "$f"
    ((count++)) || true
  done < <(find "$base" "$@" -print0 2>/dev/null)
  echo "  $count"
}

total_before=$(du -sb "$ROOT" --exclude=".git" --exclude="external" 2>/dev/null | awk '{print $1}')

echo "=== clean.sh — limpieza de artefactos regenerables ==="
[ $DRY_RUN -eq 1 ] && echo "  (modo dry-run — no se borra nada)"
echo ""

# ── 1. Caché suelta de sr-index/sdsl en la raíz + temporales ───────────────
# sr-index deja estos en el cwd al construir/benchmarquear ad-hoc.
echo "▸ Caché sr-index/sdsl y temporales (raíz):"
delete_glob "$ROOT" -maxdepth 1 \
  \( -name '*.sdsl' -o -name 'construction-*.json' -o -name 'construction-*.html' \
     -o -name '*.tmp' -o -name '*.temp' -o -name 'timing.txt' \) -type f

# ── 2. Datasets sintéticos generados (--data) ──────────────────────────────
if [ $CLEAN_DATA -eq 1 ]; then
  echo "▸ Datasets generados data/{scaling,classify}/:"
  count=0
  for d in "$ROOT/data/scaling" "$ROOT/data/classify"; do
    [ -e "$d" ] && { delete "$d"; ((count++)) || true; }
  done
  echo "  $count directorios"
  # sweep 100mb: solo el texto pesado; se conserva metadata.json versionado.
  echo "▸ reference/patterns del sweep 100mb (conserva metadata.json):"
  delete_glob "$ROOT/data/bench_100mb_mutation_sweep" \
    \( -name 'reference.txt' -o -name 'patterns.txt' \) -type f
fi

# ── 3. Resultados de benchmarks (--results) ────────────────────────────────
# Se conservan las notas .md de results/ y results/compare/ (versionadas).
if [ $CLEAN_RESULTS -eq 1 ]; then
  echo "▸ Salidas de benchmarks results/{bench,classify,scaling}/:"
  count=0
  for d in "$ROOT/results/bench" "$ROOT/results/classify" "$ROOT/results/scaling"; do
    [ -e "$d" ] && { delete "$d"; ((count++)) || true; }
  done
  echo "  $count directorios"
fi

# ── 4. build/ del proyecto (--build) ───────────────────────────────────────
if [ $CLEAN_BUILD -eq 1 ]; then
  echo "▸ Directorio build/:"
  if [ -d "$ROOT/build" ]; then
    delete "$ROOT/build"
    echo "  1 directorio"
  else
    echo "  (no existe)"
  fi
fi

echo ""
if [ $DRY_RUN -eq 0 ]; then
  total_after=$(du -sb "$ROOT" --exclude=".git" --exclude="external" 2>/dev/null | awk '{print $1}')
  freed=$(( total_before - total_after ))
  echo "✓ Limpieza completa — espacio liberado: $(numfmt --to=iec $freed 2>/dev/null || echo "${freed} bytes")"
else
  echo "Ejecuta sin --dry-run para borrar. Sugerencia: --all para todo lo regenerable."
fi
