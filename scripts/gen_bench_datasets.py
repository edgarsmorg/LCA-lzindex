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


def generate_dataset(spec, out_root):
    name = spec["name"]
    out_dir = out_root / name
    out_dir.mkdir(parents=True, exist_ok=True)

    text = build_text(spec)
    patterns = build_patterns(spec, text)

    (out_dir / "reference.txt").write_text(text, encoding="ascii")
    (out_dir / "patterns.txt").write_text("\n".join(patterns) + "\n", encoding="ascii")
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
