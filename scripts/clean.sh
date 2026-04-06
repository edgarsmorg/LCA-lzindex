#!/bin/bash
# Limpia archivos de caché generados por sr-index, sdsl-lite y ropebwt3.
# Uso: ./scripts/clean.sh [--dry-run] [--build]
#   --dry-run  Solo muestra qué se borraría, sin borrar nada
#   --build    También borra el directorio build/ del proyecto

set -e

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DRY_RUN=0
CLEAN_BUILD=0

for arg in "$@"; do
  case "$arg" in
    --dry-run) DRY_RUN=1 ;;
    --build)   CLEAN_BUILD=1 ;;
  esac
done

delete() {
  if [ $DRY_RUN -eq 1 ]; then
    echo "  [dry-run] rm $1"
  else
    rm -rf "$1"
  fi
}

total_before=$(du -sb "$ROOT" --exclude=".git" --exclude="external" 2>/dev/null | awk '{print $1}')

echo "=== clean.sh — limpieza de caché del proyecto ==="
[ $DRY_RUN -eq 1 ] && echo "  (modo dry-run — no se borra nada)"
echo ""

# ── 1. Caché sdsl-lite (.sdsl) generada en cwd por sr-index benchmarks ─────
echo "▸ Archivos sdsl-lite (.sdsl):"
count=0
while IFS= read -r -d '' f; do
  delete "$f"
  ((count++)) || true
done < <(find "$ROOT" -maxdepth 1 -name "*.sdsl" -print0 2>/dev/null)
echo "  $count archivos"

# ── 2. Logs de construcción de sr-index ────────────────────────────────────
echo "▸ Logs de construcción sr-index (construction-*.{json,html}):"
count=0
while IFS= read -r -d '' f; do
  delete "$f"
  ((count++)) || true
done < <(find "$ROOT" -maxdepth 1 \( -name "construction-*.json" -o -name "construction-*.html" \) -print0 2>/dev/null)
echo "  $count archivos"

# ── 3. Índices BWT dinámicos temporales (.fmr) ─────────────────────────────
echo "▸ Índices BWT temporales (.fmr):"
count=0
while IFS= read -r -d '' f; do
  delete "$f"
  ((count++)) || true
done < <(find "$ROOT/data" -name "*.fmr" -print0 2>/dev/null)
echo "  $count archivos"

# ── 4. Directorio sr_build en data/ ────────────────────────────────────────
echo "▸ Directorios sr_build/:"
count=0
while IFS= read -r -d '' d; do
  delete "$d"
  ((count++)) || true
done < <(find "$ROOT/data" -type d -name "sr_build" -print0 2>/dev/null)
echo "  $count directorios"

# ── 5. build/ del proyecto (opcional) ──────────────────────────────────────
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
  echo "Ejecuta sin --dry-run para borrar."
fi
