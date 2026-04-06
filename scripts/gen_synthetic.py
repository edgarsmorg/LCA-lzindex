#!/usr/bin/env python3
"""
Genera datos sintéticos para validar el pipeline LCA completo.

Árbol filogenético:
         root (0)
        /        \
      A (1)      B (2)
      / \        / \
   A1(3) A2(4) B1(5) B2(6)

Genomas diseñados para tener MEMs con LCAs conocidos:
  - "ACGTACGT" aparece en A1 y A2  → LCA = A (nodo 1)
  - "TGCATGCA" aparece en B1 y B2  → LCA = B (nodo 2)
  - "GGGGGGGG" aparece en A2 y B1  → LCA = root (nodo 0)
  - "AAAAACGT" aparece solo en A1   → LCA = A1 (nodo 3)
"""

import json
import os

OUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'data', 'synthetic')
os.makedirs(OUT_DIR, exist_ok=True)

# ── Árbol filogenético ──────────────────────────────────────────────────────
tree = {
    "nodes": [
        {"id": 0, "name": "root", "parent": None},
        {"id": 1, "name": "A",    "parent": 0},
        {"id": 2, "name": "B",    "parent": 0},
        {"id": 3, "name": "A1",   "parent": 1},
        {"id": 4, "name": "A2",   "parent": 1},
        {"id": 5, "name": "B1",   "parent": 2},
        {"id": 6, "name": "B2",   "parent": 2},
    ],
}

# ── Genomas en orden DFS (hojas izq→der) ───────────────────────────────────
# Diseñados para tener regiones compartidas específicas:
#   ACGTACGT → compartido por A1, A2 (grupo A)
#   TGCATGCA → compartido por B1, B2 (grupo B)
#   GGGGGGGG → compartido por A2, B1 (cruce → LCA = root)
#   AAAAACGT → único de A1
genomes_dfs = [
    (3, "A1", "AAAAACGTACGTACGTACGT"),          # A1: único "AAAA" + compartido "ACGTACGT"
    (4, "A2", "ACGTACGTACGTGGGGGGGG"),          # A2: "ACGTACGT" + "GGGGGGGG" (cruce)
    (5, "B1", "GGGGGGGGTTTTTTGCATGCA"),         # B1: "GGGGGGGG" (cruce) + "TGCATGCA"
    (6, "B2", "TGCATGCATGCACCCCCCCC"),          # B2: "TGCATGCA" + único "CCCCCCCC"
]

# ── Texto concatenado en orden DFS ─────────────────────────────────────────
# Separamos genomas con '#' (no-ACGT) para evitar MEMs cruzando límites
SEP = '#'
text_parts = []
genome_metadata = []
pos = 0

for node_id, name, seq in genomes_dfs:
    start = pos
    text_parts.append(seq)
    pos += len(seq)
    genome_metadata.append({
        "node_id": node_id,
        "name": name,
        "start": start,
        "end": pos,
        "length": len(seq),
    })
    text_parts.append(SEP)
    pos += 1  # separador

text = ''.join(text_parts)

# ── Queries de prueba con LCA esperado ────────────────────────────────────
queries = [
    {"id": "q1", "seq": "ACGTACGT",  "expected_lca_id": 1, "expected_lca_name": "A",    "note": "compartido A1+A2"},
    {"id": "q2", "seq": "TGCATGCA",  "expected_lca_id": 2, "expected_lca_name": "B",    "note": "compartido B1+B2"},
    {"id": "q3", "seq": "GGGGGGGG",  "expected_lca_id": 0, "expected_lca_name": "root", "note": "compartido A2+B1"},
    {"id": "q4", "seq": "AAAAACGT",  "expected_lca_id": 3, "expected_lca_name": "A1",   "note": "único A1"},
    {"id": "q5", "seq": "CCCCCCCC",  "expected_lca_id": 6, "expected_lca_name": "B2",   "note": "único B2"},
]

# ── Escritura de archivos ──────────────────────────────────────────────────
# 1. Texto concatenado (para sr-index y ropebwt3)
with open(f'{OUT_DIR}/reference.txt', 'w') as f:
    f.write(text)

# 2. FASTA del texto (para ropebwt3)
with open(f'{OUT_DIR}/reference.fa', 'w') as f:
    for node_id, name, seq in genomes_dfs:
        f.write(f'>{name}\n{seq}\n')

# 3. Árbol filogenético + metadatos
with open(f'{OUT_DIR}/tree.json', 'w') as f:
    json.dump({
        "tree": tree,
        "genomes": genome_metadata,
    }, f, indent=2)

# 4. Queries de prueba
with open(f'{OUT_DIR}/queries.json', 'w') as f:
    json.dump(queries, f, indent=2)

# 5. Resumen para verificación manual
print("=== Texto de referencia ===")
print(f"  Texto: {text!r}")
print(f"  Largo: {len(text)} chars")
print()
print("=== Genomas ===")
for g in genome_metadata:
    print(f"  [{g['start']:3d},{g['end']:3d}) {g['name']} ({g['node_id']}): {text[g['start']:g['end']]!r}")
print()
print("=== Queries y LCAs esperados ===")
for q in queries:
    # Buscar manualmente en qué genomas aparece
    ocurrencias = []
    for g in genome_metadata:
        if q['seq'] in text[g['start']:g['end']]:
            ocurrencias.append(g['name'])
    print(f"  {q['id']}: {q['seq']!r} → aparece en {ocurrencias} → LCA={q['expected_lca_name']} ({q['note']})")
print()
print(f"Archivos generados en {OUT_DIR}/")
