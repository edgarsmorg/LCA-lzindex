#!/usr/bin/env python3
"""
Genera data/repetitive/reads.tsv con reads de ~250 bp y ground-truth taxonómico.

Cada read es un fragmento extraído de un genoma específico con 1% de SNPs.
Los reads son lo suficientemente largos para que MEMExtractor con min_len=31 funcione.

Columnas de salida:
  read_id, true_genome_idx, true_node_id, sequence

Uso:
  python3 scripts/gen_classification_reads.py [data_dir] [n_reads] [read_len]
"""

import json, random, sys, os

SEED     = 42
N_READS  = 2000
READ_LEN = 250
SNP_RATE = 0.01
ALPHABET = 'ACGT'


def apply_snps(seq, rate, rng):
    seq = list(seq)
    for i in range(len(seq)):
        if rng.random() < rate:
            orig = seq[i]
            alts = [c for c in ALPHABET if c != orig]
            seq[i] = rng.choice(alts)
    return ''.join(seq)


def main():
    data_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        os.path.dirname(__file__), '..', 'data', 'repetitive')
    n_reads  = int(sys.argv[2]) if len(sys.argv) > 2 else N_READS
    read_len = int(sys.argv[3]) if len(sys.argv) > 3 else READ_LEN
    data_dir = os.path.normpath(data_dir)

    ref_path     = os.path.join(data_dir, 'reference.txt')
    tree_path    = os.path.join(data_dir, 'tree.json')
    genomes_path = os.path.join(data_dir, 'genomes.tsv')
    out_path     = os.path.join(data_dir, 'reads.tsv')

    # Cargar genomes.tsv (ya generado por tree_json_to_tsv.py)
    if not os.path.exists(genomes_path):
        print(f"ERROR: {genomes_path} no existe. Ejecutar tree_json_to_tsv.py primero.")
        sys.exit(1)

    genomes = []  # lista de (start, end, node_id, genome_idx_dfs_rank)
    with open(genomes_path) as f:
        next(f)  # header
        for line in f:
            rank, start, end, nid = line.strip().split('\t')
            genomes.append((int(start), int(end), int(nid)))

    # Cargar tree.json para obtener genome_idx por node_id
    with open(tree_path) as f:
        t = json.load(f)
    params     = t['params']
    leaf_order = t['leaf_dfs_order']  # dfs_rank → genome_idx

    # Leer referencia
    print(f"Cargando reference.txt ({os.path.getsize(ref_path):,} bytes)...", end=" ", flush=True)
    with open(ref_path, 'r') as f:
        reference = f.read()
    print("OK")

    rng = random.Random(SEED)

    n_genomes   = len(genomes)
    min_margin  = 0  # leer desde inicio del genoma
    max_margin  = 0  # hasta el final del genoma

    with open(out_path, 'w') as out:
        out.write('read_id\ttrue_genome_idx\ttrue_node_id\tsequence\n')
        for i in range(n_reads):
            # Elegir genoma aleatorio
            rank   = rng.randint(0, n_genomes - 1)
            start, end, node_id = genomes[rank]
            genome_len = end - start  # exclusivo: [start, end)

            if genome_len < read_len:
                # Genoma demasiado corto para este read_len; usar el fragmento entero
                frag_start = start
                frag_len   = genome_len
            else:
                max_pos    = genome_len - read_len
                frag_start = start + rng.randint(0, max_pos)
                frag_len   = read_len

            seq = reference[frag_start : frag_start + frag_len]

            # Verificar que no cruzamos un separador '$'
            assert '$' not in seq, f"Read {i} cruza separador en rank={rank}, pos={frag_start}"

            # Aplicar SNPs
            seq = apply_snps(seq, SNP_RATE, rng)

            genome_idx = leaf_order[rank]
            out.write(f'{i}\t{genome_idx}\t{node_id}\t{seq}\n')

    print(f"reads.tsv: {n_reads} reads de {read_len} bp ({SNP_RATE*100:.0f}% SNPs)")
    print(f"  cobertura: {n_reads} / {n_genomes} genomas (aprox {n_reads/n_genomes:.1f} reads/genoma)")

if __name__ == '__main__':
    main()
