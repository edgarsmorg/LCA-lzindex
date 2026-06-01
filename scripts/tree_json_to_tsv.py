#!/usr/bin/env python3
"""
Convierte data/repetitive/tree.json a dos TSVs consumibles desde C++:

  tree.tsv     — node_id, parent_id, name, genome_idx
                 parent_id = -1 para la raíz; genome_idx = -1 para nodos internos.

  genomes.tsv  — leaf_dfs_rank, start, end, node_id
                 start/end = rango [start, end) en el texto de referencia.

Uso:
  python3 scripts/tree_json_to_tsv.py [data_dir]
  data_dir por defecto: data/repetitive
"""

import json, sys, os

def main():
    data_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        os.path.dirname(__file__), '..', 'data', 'repetitive')
    data_dir = os.path.normpath(data_dir)

    tree_path    = os.path.join(data_dir, 'tree.json')
    tree_tsv     = os.path.join(data_dir, 'tree.tsv')
    genomes_tsv  = os.path.join(data_dir, 'genomes.tsv')
    ref_path     = os.path.join(data_dir, 'reference.txt')

    with open(tree_path) as f:
        t = json.load(f)

    nodes       = t['nodes']          # list of {id, children, genome_idx}
    root        = t['root']
    leaf_order  = t['leaf_dfs_order'] # genome_idx en orden DFS
    params      = t['params']
    genome_len  = params['genome_len']
    n_leaves    = params['n_genomes']

    # Validar tamaño del texto de referencia si existe
    if os.path.exists(ref_path):
        ref_size = os.path.getsize(ref_path)
        expected = n_leaves * genome_len + (n_leaves - 1)
        assert ref_size == expected, (
            f"Tamaño de reference.txt ({ref_size}) ≠ esperado ({expected})")
        print(f"reference.txt: {ref_size:,} bytes ✓")

    # Construir mapa parent: node_id → parent_id
    parent = {root: -1}
    for n in nodes:
        for child in n['children']:
            parent[child] = n['id']

    # Mapa genome_idx → node_id (para hojas)
    genome_to_node = {}
    for n in nodes:
        if n['genome_idx'] is not None:
            genome_to_node[n['genome_idx']] = n['id']

    # tree.tsv
    with open(tree_tsv, 'w') as f:
        f.write('node_id\tparent_id\tname\tgenome_idx\n')
        for n in nodes:
            nid   = n['id']
            pid   = parent.get(nid, -1)
            gidx  = n['genome_idx'] if n['genome_idx'] is not None else -1
            name  = f'genome_{gidx}' if gidx >= 0 else f'internal_{nid}'
            f.write(f'{nid}\t{pid}\t{name}\t{gidx}\n')

    print(f"tree.tsv: {len(nodes)} nodos")

    # genomes.tsv
    stride = genome_len + 1  # genome_len bytes + 1 separador '$'
    with open(genomes_tsv, 'w') as f:
        f.write('leaf_dfs_rank\tstart\tend\tnode_id\n')
        for rank, genome_idx in enumerate(leaf_order):
            node_id = genome_to_node[genome_idx]
            start   = rank * stride
            end     = start + genome_len      # exclusivo
            f.write(f'{rank}\t{start}\t{end}\t{node_id}\n')

    # Verificar que el último end coincide con el tamaño total del texto
    last_end = (n_leaves - 1) * stride + genome_len
    expected_size = n_leaves * genome_len + (n_leaves - 1)
    assert last_end == expected_size, f"last_end={last_end} ≠ expected={expected_size}"
    print(f"genomes.tsv: {n_leaves} rangos, último end={last_end:,} ✓")

if __name__ == '__main__':
    main()
