# Uso del Baseline: sr-index para clasificación taxonómica

Este documento describe paso a paso cómo usar el **sr-index** como baseline para
verificar la correctitud del pipeline de clasificación taxonómica. El flujo completo
parte de genomas de referencia y lecturas query, y termina con posiciones en el texto
de referencia que luego se mapean a un árbol filogenético para computar el LCA.

---

## Visión general del pipeline

```
genomas.fa (FASTA, orden DFS)
      │
      ▼
 [ropebwt3 build]
      │
      ├──────────────────────────────────────┐
      ▼                                      ▼
index.fmd                            reference.txt (plano)
(para MEM extraction)                (para sr-index)
      │                                      │
      │                                 [bm_construct_ri]
      │                                      │
      ▼                                      ▼
 [ropebwt3 mem]                       caché .sdsl
      │                                      │
      ▼                                      │
 mems.bed (query_name, start, end, count)    │
      │                                      │
      ▼                                      │
 patterns.txt (secuencias únicas de MEMs) ───┤
      │                                      │
      └───────────► [sr_locate] ◄────────────┘
                         │
                         ▼
                  locations.txt (posición en ref por patrón)
                         │
                         ▼
               mapeo pos → genoma_id → LCA
```

---

## Requisitos previos

### Compilar el proyecto

```bash
# Desde la raíz del repositorio
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel

# Verificar binarios disponibles
ls build/sr_locate                                  # locate queries
ls external/ropebwt3/ropebwt3                      # BWT + MEM extraction
ls external/sr-index/build/benchmark/sr-index/bm_construct_ri  # construcción sr-index
ls external/sr-index/build/benchmark/sr-index/bm_locate_ri     # benchmark locate
```

### Estructura de datos de entrada

El texto de referencia debe ser la **concatenación de los genomas en orden DFS** del
árbol filogenético, separados por el caracter `#` (que no aparece en los genomas):

```
genoma_A1 + '#' + genoma_A2 + '#' + genoma_B1 + '#' + genoma_B2 + '#'
```

El orden DFS es crítico: garantiza que las posiciones extremales (min/max) en el texto
correspondan a los genomas extremos en el árbol, lo que permite computar el LCA
correctamente sin enumerar todas las ocurrencias.

---

## Paso 1: Construir la BWT (ropebwt3)

ropebwt3 acepta FASTA directamente e indexa **forward + reverse complement** por
defecto. Para un índice solo-forward (recomendado para clasificación taxonómica sin
reverse complement), agregar `-R`:

```bash
# Construir índice dinámico (.fmr)
./external/ropebwt3/ropebwt3 build -t$(nproc) -bo data/mi_ref.fmr data/genomas.fa

# Convertir a formato estático (.fmd) — más rápido para queries
./external/ropebwt3/ropebwt3 build -i data/mi_ref.fmr -do data/mi_ref.fmd

# Limpiar temporal
rm data/mi_ref.fmr
```

**Formatos:**
- `.fmr` — B+-tree dinámico, soporta inserciones online
- `.fmd` — memoria-mappeable, más rápido para queries, recomendado en producción

---

## Paso 2: Construir el sr-index (bm_construct_ri)

El sr-index necesita el texto en **formato plano** (sin cabeceras FASTA). El parámetro
`s` controla el tradeoff espacio/tiempo: mayor `s` = menos espacio, más lento el locate.

```bash
# Convertir FASTA a texto plano con separadores '#'
python3 -c "
import re, sys
fasta = open('data/genomas.fa').read()
seqs = [s.replace('\n','').replace(' ','') for s in re.split(r'>.*\n', fasta)[1:]]
open('data/reference.txt', 'w').write('#'.join(seqs) + '#')
"

# Construir sr-index con s=16 (buen balance espacio/tiempo)
# IMPORTANTE: ejecutar desde el directorio donde quieres los archivos .sdsl
mkdir -p data/sr_build
(cd data/sr_build && \
 ../../external/sr-index/build/benchmark/sr-index/bm_construct_ri \
     -data=../../data/reference.txt \
     -min_s=16 \
     -max_s=16 \
     --benchmark_format=json \
     --benchmark_out=../../results/construct.json)
```

### Selección del parámetro s

| s  | Espacio (relativo) | Locate time | Recomendación |
|----|-------------------|-------------|---------------|
| 4  | Mayor             | Más rápido  | Benchmarks de velocidad |
| 16 | Moderado          | Moderado    | **Balance recomendado** |
| 64 | Menor             | Más lento   | Espacio crítico |
| 128| Mínimo            | Lento       | Solo si el espacio es muy limitado |

La variante `SrIndexValidArea` es la recomendada (mejor tradeoff espacio/tiempo).

### Archivos generados (.sdsl)

`bm_construct_ri` genera ~30 archivos `.sdsl` en el directorio donde se ejecuta.
Estos son la caché del índice — si ya existen, no es necesario reconstruir.

```bash
# Verificar que el índice está construido (archivo clave)
ls data/sr_build/16_bwt_run_last_text_pos_reference.txt.sdsl
```

---

## Paso 3: Extraer MEMs (ropebwt3 mem)

Los MEMs (Maximal Exact Matches) son las semillas que conectan una lectura query con
el texto de referencia.

```bash
# Extraer MEMs con longitud mínima 19 (default)
./external/ropebwt3/ropebwt3 mem \
    -t$(nproc) \
    -l19 \
    data/mi_ref.fmd \
    data/queries.fa \
    > data/mems.bed 2>/dev/null
```

**Formato de salida (BED-like):**
```
query_name   start   end   count
seq1         0       45    3
seq1         12      67    1
seq2         5       38    7
```
- `start`, `end`: posiciones en la secuencia **query** (0-based, half-open `[start, end)`)
- `count`: número de ocurrencias en la referencia

**Parámetros importantes:**
```
-l INT    Longitud mínima del MEM (default: 19, recomendado ≥ 31 para datos reales)
-c INT    Intervalo BWT mínimo (default: 1; usar -c2 para ignorar singletons)
-R        Indexar solo forward strand (sin reverse complement)
-p INT    Imprimir hasta INT posiciones en el texto de referencia (requiere .fmd estático)
```

---

## Paso 4: Locate con sr_locate

`sr_locate` toma patrones (uno por línea) y devuelve sus posiciones en el texto de
referencia usando el sr-index pre-construido.

```bash
# Extraer secuencias de MEMs desde el FASTA de queries
python3 - <<'PYEOF'
import re

seqs = {}
current = None
for line in open('data/queries.fa'):
    line = line.rstrip()
    if line.startswith('>'):
        current = line[1:].split()[0]
        seqs[current] = []
    elif current:
        seqs[current].append(line)
seqs = {k: ''.join(v) for k, v in seqs.items()}

patterns = set()
for line in open('data/mems.bed'):
    parts = line.rstrip().split('\t')
    if len(parts) < 3:
        continue
    name, start, end = parts[0], int(parts[1]), int(parts[2])
    if name in seqs:
        patterns.add(seqs[name][start:end])

with open('data/patterns.txt', 'w') as f:
    for p in sorted(patterns):
        f.write(p + '\n')

print(f"Patrones únicos: {len(patterns)}")
PYEOF

# Ejecutar sr_locate
./build/sr_locate \
    --data_dir=data/sr_build \
    --data_name=reference.txt \
    --patterns=data/patterns.txt \
    > data/locations.txt
```

**Formato de salida:**
```
>ACGTACGTACGT
4
21
45
>TTGAAGGAGTCT
100
```
Cada línea `>PATRON` es seguida por sus posiciones en el texto de referencia (una por línea).

---

## Paso 5: Mapear posiciones a genomas y LCA

Con las posiciones en `locations.txt`, el paso final es:
1. Determinar a qué genoma pertenece cada posición (usando los límites del árbol filogenético)
2. Tomar la posición mínima y máxima para cada patrón
3. Computar el LCA de los genomas extremos

```python
# Ejemplo simplificado (ver scripts/test_synthetic_lca.py para versión completa)
def genome_of(pos, genomes):
    for g in genomes:
        if g['start'] <= pos < g['end']:
            return g['node_id']
    return None  # separador '#'

def lca_from_locations(positions, genomes, tree):
    nodes = [genome_of(p, genomes) for p in positions]
    nodes = [n for n in nodes if n is not None]
    if not nodes:
        return None
    # Versión extremal (solo min y max)
    return lca(nodes[0], nodes[-1], tree)   # nodes[0]=pos_min, nodes[-1]=pos_max
```

---

## Uso con el script automatizado

Para ejecutar el pipeline completo de una sola vez:

```bash
./scripts/run_baseline.sh <reference.fa> <queries.fa> [outdir] [min_mem_len] [sr_s]

# Ejemplo con datos sintéticos
./scripts/run_baseline.sh \
    data/synthetic/reference.fa \
    data/test_seq.fa \
    results/baseline_synthetic \
    8 \
    16
```

---

## Benchmark comparativo

Para medir el rendimiento del sr-index (tiempo de locate y espacio):

```bash
# Construir con múltiples valores de s y comparar
(cd data/sr_build && \
 ../../external/sr-index/build/benchmark/sr-index/bm_locate_ri \
     -data_dir=. \
     -data_name=reference.txt \
     -patterns=../../data/patterns.txt \
     -min_s=4 \
     -max_s=64 \
     --benchmark_format=json \
     --benchmark_out=../../results/locate_benchmark.json)

# Medir espacio (bits por carácter)
python3 -c "
import os, json
n = os.path.getsize('data/reference.txt')
idx_size = sum(os.path.getsize(f) for f in __import__('glob').glob('data/sr_build/*.sdsl'))
print(f'n = {n} bytes')
print(f'índice = {idx_size} bytes ({idx_size/n:.2f}x)')
print(f'bpc = {idx_size*8/n:.2f}')
"
```

---

## Verificación de correctitud

Para verificar que sr_locate devuelve posiciones correctas, comparar contra búsqueda
brute-force en el texto plano:

```bash
python3 -c "
text = open('data/reference.txt').read()
for pattern in open('data/patterns.txt'):
    pattern = pattern.rstrip()
    if not pattern: continue
    positions = []
    start = 0
    while True:
        idx = text.find(pattern, start)
        if idx == -1: break
        positions.append(idx)
        start = idx + 1
    print(f'>{pattern}')
    for p in positions: print(p)
" > data/locations_brute.txt

diff data/locations.txt data/locations_brute.txt && echo "✓ sr_locate correcto" || echo "✗ hay diferencias"
```

---

## Troubleshooting

### sr_locate no encuentra el índice
```
Error cargando el índice
```
Verificar que `--data_dir` apunta al directorio con los archivos `.sdsl` y que
`--data_name` coincide exactamente con el nombre usado al construir (incluyendo `.txt`).

### bm_construct_ri genera .sdsl en cwd
`bm_construct_ri` siempre guarda los archivos `.sdsl` en el directorio de trabajo actual.
Ejecutar siempre desde el directorio destino:
```bash
(cd mi_directorio && /ruta/completa/bm_construct_ri -data=/ruta/completa/ref.txt ...)
```

### MEMs muy cortos / demasiados falsos positivos
Aumentar `-l` en `ropebwt3 mem`. Para datos reales de ADN, `-l31` es el mínimo práctico.
Para datos sintéticos pequeños, `-l8` o `-l10` puede ser necesario.

### Texto de referencia vs FASTA
`bm_construct_ri` requiere texto plano, **no FASTA**. Si se pasa un FASTA directamente,
el índice incluirá las cabeceras `>seq_name` como parte del texto, produciendo
resultados incorrectos.
