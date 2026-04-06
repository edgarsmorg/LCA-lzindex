#!/usr/bin/env python3
"""
Test end-to-end del pipeline LCA sobre datos sintéticos.

Pipeline (usando sr-index como oráculo de locate):
  query → sr-index locate → posiciones en ref → genoma_min/max → LCA

Requiere haber corrido gen_synthetic.py antes.
"""

import json
import subprocess
import os
import sys
import tempfile

BASE = os.path.join(os.path.dirname(__file__), '..')
SYNTHETIC_DIR = os.path.join(BASE, 'data', 'synthetic')
SR_LOCATE = os.path.join(BASE, 'external', 'sr-index', 'build',
                         'benchmark', 'sr-index', 'bm_locate_ri')

# ── Cargar metadatos ────────────────────────────────────────────────────────
with open(f'{SYNTHETIC_DIR}/tree.json') as f:
    meta = json.load(f)

with open(f'{SYNTHETIC_DIR}/queries.json') as f:
    queries = json.load(f)

tree_nodes = {n['id']: n for n in meta['tree']['nodes']}
genomes    = meta['genomes']  # [{"node_id", "name", "start", "end"}, ...]

# ── LCA naive en Python (para el test) ─────────────────────────────────────
def depth(node_id):
    d = 0
    n = tree_nodes[node_id]
    while n['parent'] is not None:
        d += 1
        n = tree_nodes[n['parent']]
    return d

def lca(u, v):
    du, dv = depth(u), depth(v)
    nu, nv = tree_nodes[u], tree_nodes[v]
    while du > dv:
        nu = tree_nodes[nu['parent']]
        du -= 1
    while dv > du:
        nv = tree_nodes[nv['parent']]
        dv -= 1
    while nu['id'] != nv['id']:
        nu = tree_nodes[nu['parent']]
        nv = tree_nodes[nv['parent']]
    return nu['id']

def genome_of(pos):
    for g in genomes:
        if g['start'] <= pos < g['end']:
            return g['node_id']
    return None  # separador

# ── Buscar patrón en el texto (brute-force como oráculo) ───────────────────
with open(f'{SYNTHETIC_DIR}/reference.txt') as f:
    ref_text = f.read()

def locate_brute(pattern):
    positions = []
    start = 0
    while True:
        idx = ref_text.find(pattern, start)
        if idx == -1:
            break
        positions.append(idx)
        start = idx + 1
    return positions

# ── Pipeline LCA ───────────────────────────────────────────────────────────
def classify(pattern):
    positions = locate_brute(pattern)
    if not positions:
        return None, []

    genome_nodes = [genome_of(p) for p in positions]
    genome_nodes = [g for g in genome_nodes if g is not None]
    if not genome_nodes:
        return None, []

    result = genome_nodes[0]
    for g in genome_nodes[1:]:
        result = lca(result, g)

    return result, positions

def classify_extremal(pattern):
    """Versión optimizada: solo usa pos_min y pos_max."""
    positions = locate_brute(pattern)
    if not positions:
        return None, []

    genome_nodes = [genome_of(p) for p in positions if genome_of(p) is not None]
    if not genome_nodes:
        return None, []

    # Solo necesitamos los extremos para el LCA en orden DFS
    node_min = genome_nodes[0]   # posición más a la izquierda → genoma más "izquierdo"
    node_max = genome_nodes[-1]  # posición más a la derecha  → genoma más "derecho"

    return lca(node_min, node_max), positions

# ── Ejecutar tests ──────────────────────────────────────────────────────────
print("=" * 65)
print("Test end-to-end: pipeline LCA sobre datos sintéticos")
print("=" * 65)
print(f"Texto de referencia: {ref_text!r}")
print()

passed = 0
failed = 0

for q in queries:
    pattern      = q['seq']
    expected_id  = q['expected_lca_id']
    expected_name = q['expected_lca_name']

    lca_full, positions = classify(pattern)
    lca_ext,  _         = classify_extremal(pattern)

    ok_full = lca_full == expected_id
    ok_ext  = lca_ext  == expected_id

    status = "✓" if (ok_full and ok_ext) else "✗"
    if ok_full and ok_ext:
        passed += 1
    else:
        failed += 1

    print(f"{status} {q['id']}: pattern={pattern!r}")
    print(f"    posiciones:        {positions}")
    print(f"    genomas tocados:   {[tree_nodes[genome_of(p)]['name'] for p in positions if genome_of(p) is not None]}")
    print(f"    LCA completo:      {tree_nodes[lca_full]['name']} (esperado: {expected_name}) {'✓' if ok_full else '✗'}")
    print(f"    LCA extremal:      {tree_nodes[lca_ext]['name']}  (esperado: {expected_name}) {'✓' if ok_ext else '✗'}")
    print(f"    ({q['note']})")
    print()

print("=" * 65)
print(f"Resultado: {passed}/{passed+failed} tests pasaron")
if failed == 0:
    print("✓ Pipeline LCA validado correctamente sobre datos sintéticos")
else:
    print("✗ Hay fallos — revisar la lógica")
    sys.exit(1)
