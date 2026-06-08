#!/usr/bin/env python3
"""Print compact aggregate tables from microbench CSV outputs."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def read_csv(path):
    if not path.exists():
        return []
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def mean(values):
    values = [float(v) for v in values if v not in ("", None)]
    return sum(values) / len(values) if values else 0.0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("results_dir")
    args = parser.parse_args()
    root = Path(args.results_dir)

    print("SIZE")
    for row in read_csv(root / "size.csv"):
        print(f"{row['dataset']}\t{row['index']}\ts={row['s']}\tbpc={row['bpc']}\tbytes={row['serialized_bytes']}")

    locate = defaultdict(list)
    for row in read_csv(root / "locate.csv"):
        locate[(row["dataset"], row["index"], row["s"])].append(row["us_per_query"])
    print("\nLOCATE mean us/query")
    for key, vals in sorted(locate.items()):
        print(f"{key[0]}\t{key[1]}\ts={key[2]}\t{mean(vals):.3f}")


if __name__ == "__main__":
    main()
