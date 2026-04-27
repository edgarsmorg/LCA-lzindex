#!/usr/bin/env python3
"""
Genera datos sintéticos GRANDES para medir overhead de WtMinRmq con z controlado.

Estrategia: generar genomas altamente repetitivos (DNA-like).
Patrón de 1KB se repite con mutaciones ocasionales.
Esperar z natural ≈ 10-100k dependiendo del tamaño base.
"""

import random
import os
import json

OUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'data', 'synthetic')
os.makedirs(OUT_DIR, exist_ok=True)

def generate_dna(size_bytes, pattern_len=1000, mutation_rate=0.01):
    """Generar secuencia DNA altamente repetitiva para bajar z."""
    alphabet = "ACGT"
    pattern = ''.join(random.choice(alphabet) for _ in range(pattern_len))

    seq = ""
    pos = 0
    while pos < size_bytes:
        # Agregar patrón con mutaciones ocasionales
        mutated = list(pattern)
        for i in range(len(mutated)):
            if random.random() < mutation_rate:
                mutated[i] = random.choice(alphabet)
        seq += ''.join(mutated)
        pos += pattern_len

    return seq[:size_bytes]

# Generar 5 datasets con z aproximado esperado:
# ~1 MB → z ≈ 10k
# ~5 MB → z ≈ 50k
# ~10 MB → z ≈ 100k
# ~20 MB → z ≈ 200k
# ~30 MB → z ≈ 300k

configs = [
    {"name": "z_10k", "size": 1_000_000, "desc": "1 MB, expect z≈10k"},
    {"name": "z_50k", "size": 5_000_000, "desc": "5 MB, expect z≈50k"},
    {"name": "z_100k", "size": 10_000_000, "desc": "10 MB, expect z≈100k"},
]

for cfg in configs:
    seq = generate_dna(cfg['size'])

    # Guardar como FASTA
    fasta_path = f'{OUT_DIR}/{cfg["name"]}.fa'
    with open(fasta_path, 'w') as f:
        f.write(f'>{cfg["name"]}\n')
        # Escribir en líneas de 80 chars (estándar FASTA)
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + '\n')

    print(f"✓ {cfg['name']:15} ({cfg['desc']:25}): {fasta_path}")
    print(f"  Size: {len(seq):,} bytes")

print(f"\nNext: build indexes with:")
for cfg in configs:
    print(f"  ./build/bin/measure_index data/synthetic/{cfg['name']}.fa")
