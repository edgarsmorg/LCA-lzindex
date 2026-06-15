# LCA-lzindex

Índice comprimido basado en LZ77 para clasificación taxonómica de lecturas de ADN.

Memoria de grado — Edgar S. Morales González, Ingeniería Civil en Computación, Universidad de Chile, DCC.
Profesor guía: Gonzalo Navarro.

---

## Qué hace

Dado un conjunto de genomas de referencia concatenados en orden DFS de un árbol filogenético, construye un índice que permite:

1. Encontrar las **posiciones extremales** (mín/máx) de las ocurrencias de un patrón de ADN en el texto concatenado.
2. Mapear esas posiciones extremales al **LCA** (último ancestro común) en el árbol filogenético.

La idea central es que para clasificación taxonómica **no necesitamos enumerar todas las ocurrencias**, solo las extremales, y estas quedan capturadas por las **ocurrencias primarias** del parsing LZ77. Esto permite descartar las estructuras del LZ-index original para secundarias (~4z·log n bits ahorrados).

El baseline de comparación es el [sr-index](https://github.com/duscob/sr-index) de Cobas, Gagie y Navarro.

---

## Requisitos

- GCC 12+, C++17
- CMake 3.20+
- [sdsl-lite (fork duscob)](https://github.com/duscob/sdsl-lite) — instalado en `$HOME/.local`
- [ropebwt3](https://github.com/lh3/ropebwt3) — en `external/ropebwt3/`
- [sr-index](https://github.com/duscob/sr-index) — en `external/sr-index/`

### Setup inicial

```bash
# Clonar el repo con todos los submódulos
git clone --recurse-submodules https://github.com/edgarsmorg/LCA-lzindex
# o, si ya se clonó sin --recurse-submodules:
git submodule update --init --recursive

# Instalar sdsl-lite (header-only + divsufsort) en $HOME/.local
cd external/sdsl-lite && ./install.sh $HOME/.local && cd ../..

# Compilar ropebwt3
make -j$(nproc) -C external/ropebwt3

# Compilar sr-index
cmake -S external/sr-index -B external/sr-index/build -DCMAKE_BUILD_TYPE=Release
cmake --build external/sr-index/build --parallel $(nproc)
```

---

## Compilación

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel $(nproc)
```

Flags usados automáticamente: `-O3 -DNDEBUG -march=native -msse4.2 -mbmi -mbmi2`.

---

## Tests

```bash
cd build && ctest --output-on-failure -j$(nproc)
```

Los tests verifican correctitud contra un FM-index estándar de sdsl-lite como ground truth (≥10k patrones aleatorios por dataset).

---

## Herramientas principales

Todos los binarios quedan en `build/`.

### Construir el LZ77-index

```bash
# Texto plano (sin headers FASTA, sin saltos de línea)
grep -v "^>" data/genomes.fa | tr -d '\n' > data/reference.txt

build/build_lz_index data/reference.txt indexes/migenoma
# Produce: indexes/migenoma.meta  .grid  .rcsa_fwd  .rcsa_rev
```

### Construir el sr-index (baseline)

```bash
build/build_sr_index data/reference.txt indexes/sr_s16 --s=16
```

### Medir espacio

```bash
# LZ77-index: desglose por componente (bpc total + grid + CSA)
build/inspect_lz_index data/reference.txt indexes/migenoma

# sr-index: footprint real del índice cargado (usa size_in_bytes(), NO suma .sdsl del dir)
build/inspect_sr_index data/reference.txt indexes/sr_s16 --s=16
```

### Benchmark de locate (LZ77 vs sr-index)

```bash
build/bench_locate_compare \
    data/reference.txt \
    data/patterns.txt \
    indexes/migenoma \
    indexes/sr_s16 \
    --s=16 --reps=3 --csv=results/locate.csv
```

### Pipeline LCA con LZ77-index

```bash
build/baseline_lca_lz \
    data/reference.txt \
    indexes/migenoma \
    data/patterns.txt \
    --mode=locate
```

### Pipeline LCA con sr-index (baseline end-to-end)

```bash
build/baseline_lca_sr \
    indexes/sr_s16 \
    migenoma \
    data/patterns.txt \
    --s=16 --mode=locate
```

---

## Datasets sintéticos

Los benchmarks de la tesis usan colecciones sintéticas de 100 MiB: 100 copias de un genoma ancestral aleatorio de 1 MiB con tasa de mutación variable (1 %–0.01 %).

```bash
python3 scripts/gen_synthetic_large.py  # ver --help para parámetros
```

Para datasets repetitivos conocidos (Pizza & Chili):

```bash
bash scripts/download_corpus.sh
```

---

## Estructura del repositorio

```
src/          Código del índice (lz77/, wavelet/, mem/, taxonomy/, baseline/)
tools/        Binarios auxiliares (build, inspect, bench, compare)
tests/        Tests unitarios con GoogleTest
scripts/      Pipeline end-to-end, generación de datos, benchmarks
data/         Datasets (gitignore)
results/      Resultados de benchmarks y resúmenes de sesión
docs/         Glosario, papers, memoria LaTeX
external/     Dependencias (sdsl-lite, ropebwt3, sr-index)
```

Nomenclatura y estructuras de datos: [docs/glosario.md](docs/glosario.md).

---

## Referencias clave

- Kreft & Navarro, *Self-indexing Based on LZ77*, TCS 2013 — diseño del índice
- Cobas, Gagie & Navarro, *sr-index*, 2024 — baseline
- Draesslerová et al., *MEMs for taxonomic classification*, SEA 2024 — caso de uso
