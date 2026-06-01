#!/usr/bin/env python3
"""
Genera corpus sintético para demostrar ventaja de LZ77-Index vs SR-Index.

Diseño:
  - Unidad base de UNIT_LEN bp (ACGT aleatorio)
  - Genoma ancestral = UNIT_REPS repeticiones de la unidad base
  - N_GENOMES copias con MUT_RATE SNPs independientes por copia
  - Árbol binario balanceado; concatenación en orden DFS

Por qué este corpus favorece al LZ77-Index:
  - z ≈ r ≈ N_GENOMES × MUT_RATE × GENOME_LEN   (pocas mutaciones → pocos z, r)
  - occ de un 20-mer de la unidad base ≈ UNIT_REPS × N_GENOMES     (alto)
  - locate_extremal LZ77: O(m^2 log^2 z)  — no depende de occ
  - locate_all SR-Index : O(m + occ × s)  — crece linealmente con occ
  → speedup esperado ≈ occ*s / (m^2 log^2 z) ≈ 8x para los parámetros por defecto
"""

import random, os, json, math, sys, textwrap

SEED        = 42
UNIT_LEN    = 1_000        # largo de la unidad base (bp)
UNIT_REPS   = 100          # repeticiones internas del genoma ancestral
N_GENOMES   = 500          # número de hojas en el árbol
MUT_RATE    = 0.0001       # SNPs por bp por copia (0.01%)
ALPHABET    = "ACGT"

OUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'data', 'repetitive')
os.makedirs(OUT_DIR, exist_ok=True)

random.seed(SEED)

GENOME_LEN  = UNIT_LEN * UNIT_REPS   # 100 KB
N_MUTS_PER  = int(GENOME_LEN * MUT_RATE)   # SNPs por copia

# ── 1. Unidad base ───────────────────────────────────────────────────────────
unit_base = ''.join(random.choice(ALPHABET) for _ in range(UNIT_LEN))
ancestor  = unit_base * UNIT_REPS    # 100KB perfectamente periódico

print(f"Unidad base      : {UNIT_LEN} bp")
print(f"Genoma ancestral : {GENOME_LEN:,} bp  ({UNIT_REPS} repeticiones)")
print(f"Copias           : {N_GENOMES}")
print(f"SNPs por copia   : {N_MUTS_PER}  ({MUT_RATE*100:.3f}%)")
print(f"z estimado       : {N_GENOMES * N_MUTS_PER * 2:,}")
print(f"occ(20-mer base) : {UNIT_REPS * N_GENOMES:,}  (100 reps × {N_GENOMES} copias)")
print()

# ── 2. Árbol balanceado (binario, N_GENOMES hojas) ───────────────────────────
def build_tree(n_leaves):
    """
    Retorna lista de nodos: cada nodo es dict con 'id', 'children', 'genome_idx'.
    Árbol binario completo; hojas extra asignadas en el último nivel.
    """
    nodes = []
    counter = [0]

    def new_node(genome_idx=None):
        nid = counter[0]; counter[0] += 1
        node = {"id": nid, "children": [], "genome_idx": genome_idx}
        nodes.append(node)
        return nid

    def build(depth, max_depth, leaves_left, leaf_counter):
        if depth == max_depth or leaves_left[0] == 0:
            idx = leaf_counter[0]; leaf_counter[0] += 1
            leaves_left[0] -= 1
            return new_node(genome_idx=idx)
        nid = new_node()
        for _ in range(2):
            if leaves_left[0] > 0:
                child = build(depth + 1, max_depth, leaves_left, leaf_counter)
                nodes[nid]["children"].append(child)
        return nid

    depth = math.ceil(math.log2(n_leaves))
    leaves_left = [n_leaves]
    leaf_counter = [0]
    build(0, depth, leaves_left, leaf_counter)
    return nodes

nodes = build_tree(N_GENOMES)

# DFS order de hojas (orden de concatenación)
def dfs_leaves(nodes, root=0):
    stack = [root]
    order = []
    while stack:
        nid = stack.pop()
        n = nodes[nid]
        if not n["children"]:
            order.append(n["genome_idx"])
        else:
            for c in reversed(n["children"]):
                stack.append(c)
    return order

leaf_order = dfs_leaves(nodes)  # orden DFS de los genomas

# ── 3. Generar copias con mutaciones ─────────────────────────────────────────
print("Generando genomas...", end=" ", flush=True)
genomes = []
for i in range(N_GENOMES):
    genome = list(ancestor)
    positions = random.sample(range(GENOME_LEN), N_MUTS_PER)
    for pos in positions:
        orig = genome[pos]
        alt  = random.choice([c for c in ALPHABET if c != orig])
        genome[pos] = alt
    genomes.append(''.join(genome))
print("OK")

# ── 4. Texto concatenado en DFS order ────────────────────────────────────────
SEP = '$'   # separador único fuera de ACGT para aislar genomas en el SA
text = SEP.join(genomes[leaf_order[i]] for i in range(N_GENOMES))
n_bytes = len(text)
print(f"Texto total      : {n_bytes:,} bytes  ({n_bytes/1e6:.1f} MB)")

with open(os.path.join(OUT_DIR, 'reference.txt'), 'w') as f:
    f.write(text)

# ── 5. Árbol JSON ─────────────────────────────────────────────────────────────
tree_data = {
    "nodes": nodes,
    "root": 0,
    "n_leaves": N_GENOMES,
    "leaf_dfs_order": leaf_order,
    "params": {
        "unit_len": UNIT_LEN,
        "unit_reps": UNIT_REPS,
        "n_genomes": N_GENOMES,
        "mut_rate": MUT_RATE,
        "n_muts_per_copy": N_MUTS_PER,
        "genome_len": GENOME_LEN,
        "n_total_bytes": n_bytes,
        "z_estimate": N_GENOMES * N_MUTS_PER * 2,
        "occ_estimate_20mer": UNIT_REPS * N_GENOMES,
    }
}
with open(os.path.join(OUT_DIR, 'tree.json'), 'w') as f:
    json.dump(tree_data, f, indent=2)

# ── 6. Queries: 20-mers de la unidad base ────────────────────────────────────
# Tomamos 1000 posiciones del centro del unit_base donde es improbable mutar
# (evitar bordes donde una SNP podría afectar ambos lados de la unidad)
N_QUERIES = 1000
QLEN      = 20
pop_size   = UNIT_LEN - QLEN   # 980 posiciones disponibles en la unidad base
positions  = random.choices(range(pop_size), k=N_QUERIES)   # con reemplazo
queries    = [unit_base[p:p+QLEN] for p in positions]

with open(os.path.join(OUT_DIR, 'queries.txt'), 'w') as f:
    f.write('\n'.join(queries) + '\n')

# ── 7. Resumen esperado ───────────────────────────────────────────────────────
z_est  = N_GENOMES * N_MUTS_PER * 2
r_est  = z_est   # en este corpus z ≈ r
occ    = UNIT_REPS * N_GENOMES
m      = QLEN
logz   = math.log2(max(z_est, 2))
lz_ops = m*m * logz**2
sr_ops = m + occ * 16

print()
print("=== Predicción teórica para 20-mers de la unidad base ===")
print(f"  z estimado            : {z_est:,}")
print(f"  r estimado            : {r_est:,}  (z ≈ r ✓)")
print(f"  occ(20-mer)           : {occ:,}")
print(f"  LZ77 locate_extremal  : {lz_ops:,.0f} ops  O(m^2 log^2 z)")
print(f"  SR   locate_all→min/max: {sr_ops:,.0f} ops  O(m + occ*16)")
print(f"  Speedup esperado      : {sr_ops/lz_ops:.1f}x  a favor de LZ77")
print()
print(f"Archivos generados en {OUT_DIR}/:")
print(f"  reference.txt  ({n_bytes/1e6:.1f} MB)")
print(f"  tree.json")
print(f"  queries.txt    ({N_QUERIES} patrones de {QLEN} bp)")
