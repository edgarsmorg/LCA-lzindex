#!/usr/bin/env python3
"""Generate deterministic plain-text benchmark datasets from a JSON manifest."""

import argparse
import json
import random
from pathlib import Path

DNA = "ACGT"


def random_dna(rng, n):
    return "".join(rng.choice(DNA) for _ in range(n))


def mutate_text(rng, text, rate):
    out = []
    for ch in text:
        if rng.random() < rate:
            choices = [c for c in DNA if c != ch]
            out.append(rng.choice(choices))
        else:
            out.append(ch)
    return "".join(out)


def build_text(spec):
    rng = random.Random(spec.get("seed", 1))
    dtype = spec["type"]
    if dtype == "random":
        return random_dna(rng, int(spec["size_bytes"]))
    if dtype == "periodic":
        unit = spec.get("unit", "ACGT")
        n = int(spec["size_bytes"])
        return (unit * ((n + len(unit) - 1) // len(unit)))[:n]
    if dtype == "mutated_copies":
        genome_len = int(spec["genome_len"])
        copies = int(spec["copies"])
        rate = float(spec.get("mutation_rate", 0.01))
        ancestor = random_dna(rng, genome_len)
        return "".join(mutate_text(rng, ancestor, rate) for _ in range(copies))
    raise ValueError(f"unknown dataset type: {dtype}")


def build_patterns(spec, text):
    rng = random.Random(int(spec.get("seed", 1)) + 1000003)
    count = int(spec.get("pattern_count", 1000))
    plen = int(spec.get("pattern_len", 16))
    patterns = []
    max_pos = max(0, len(text) - plen)
    for i in range(count):
        if i % 10 == 9:
            patterns.append(random_dna(rng, plen))
        else:
            pos = rng.randint(0, max_pos)
            patterns.append(text[pos:pos + plen])
    return patterns


# ── Artefactos filogenéticos para datasets mutated_copies ─────────────────────
# Cada copia mutada es un genoma = hoja de un árbol binario balanceado. El orden
# izquierda→derecha de las hojas coincide con el orden de concatenación (= orden
# DFS), requisito del diseño primary-only del LZ77-index.

def build_balanced_tree(copies):
    """Árbol binario balanceado sobre `copies` hojas.

    Devuelve (nodes, rank_to_node, root_id) donde:
      nodes        — lista de dicts {id, parent, genome_idx} en orden de id (0..N-1);
                     genome_idx = rank de la hoja, o -1 para nodos internos.
      rank_to_node — rank de genoma (0..copies-1) → node_id de su hoja.
    Los node_id son secuenciales 0..N-1 (requisito de PhyloTree::build), y todo
    padre tiene id menor que sus hijos.
    """
    nodes = []
    rank_to_node = [None] * copies
    counter = [0]

    def new_id():
        i = counter[0]
        counter[0] += 1
        return i

    def build(lo, hi, parent):
        nid = new_id()
        if hi - lo == 1:
            nodes.append({"id": nid, "parent": parent, "genome_idx": lo})
            rank_to_node[lo] = nid
        else:
            nodes.append({"id": nid, "parent": parent, "genome_idx": -1})
            mid = (lo + hi) // 2
            build(lo, mid, nid)
            build(mid, hi, nid)
        return nid

    root = build(0, copies, -1)
    return nodes, rank_to_node, root


def write_tree_tsv(path, nodes):
    lines = ["node_id\tparent_id\tname\tgenome_idx"]
    for n in nodes:
        gidx = n["genome_idx"]
        name = f"genome_{gidx}" if gidx >= 0 else f"internal_{n['id']}"
        lines.append(f"{n['id']}\t{n['parent']}\t{name}\t{gidx}")
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def write_genomes_tsv(path, genome_len, copies, rank_to_node):
    """Rangos contiguos [rank*genome_len, (rank+1)*genome_len), sin separador."""
    lines = ["leaf_dfs_rank\tstart\tend\tnode_id"]
    for rank in range(copies):
        start = rank * genome_len
        end = start + genome_len
        lines.append(f"{rank}\t{start}\t{end}\t{rank_to_node[rank]}")
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def write_reads_tsv(path, spec, text, genome_len, copies, rank_to_node):
    rng = random.Random(int(spec.get("seed", 1)) + 7000003)
    n_reads = int(spec.get("n_reads", 2000))
    read_len = int(spec.get("read_len", 250))
    snp_rate = float(spec.get("read_mut_rate", 0.01))

    lines = ["read_id\ttrue_genome_idx\ttrue_node_id\tsequence"]
    for i in range(n_reads):
        rank = rng.randint(0, copies - 1)
        g_start = rank * genome_len
        if genome_len <= read_len:
            frag_start, frag_len = g_start, genome_len
        else:
            frag_start = g_start + rng.randint(0, genome_len - read_len)
            frag_len = read_len
        seq = mutate_text(rng, text[frag_start:frag_start + frag_len], snp_rate)
        lines.append(f"{i}\t{rank}\t{rank_to_node[rank]}\t{seq}")
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def write_phylo_artifacts(spec, text, out_dir):
    """Para mutated_copies: emite tree.tsv, genomes.tsv y reads.tsv."""
    genome_len = int(spec["genome_len"])
    copies = int(spec["copies"])
    nodes, rank_to_node, _ = build_balanced_tree(copies)
    write_tree_tsv(out_dir / "tree.tsv", nodes)
    write_genomes_tsv(out_dir / "genomes.tsv", genome_len, copies, rank_to_node)
    write_reads_tsv(out_dir / "reads.tsv", spec, text, genome_len, copies, rank_to_node)


def generate_dataset(spec, out_root):
    name = spec["name"]
    out_dir = out_root / name
    out_dir.mkdir(parents=True, exist_ok=True)

    text = build_text(spec)
    patterns = build_patterns(spec, text)

    (out_dir / "reference.txt").write_text(text, encoding="ascii")
    (out_dir / "patterns.txt").write_text("\n".join(patterns) + "\n", encoding="ascii")

    # Artefactos filogenéticos solo para mutated_copies (tienen estructura de genomas).
    if spec["type"] == "mutated_copies":
        write_phylo_artifacts(spec, text, out_dir)

    metadata = dict(spec)
    metadata["n_bytes"] = len(text)
    metadata["reference"] = "reference.txt"
    metadata["patterns"] = "patterns.txt"
    (out_dir / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="ascii")
    return out_dir


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", default="bench/datasets.json")
    parser.add_argument("--out", default="data/bench")
    args = parser.parse_args()

    manifest = json.loads(Path(args.manifest).read_text(encoding="utf-8"))
    out_root = Path(args.out)
    for spec in manifest.get("datasets", []):
        out_dir = generate_dataset(spec, out_root)
        print(out_dir)


if __name__ == "__main__":
    main()
