# CLAUDE.md — Memoria: Mapeo de secuencias de ADN a árboles filogenéticos usando compresión Lempel-Ziv

## Identidad del proyecto

Este es el proyecto de memoria de Edgar S. Morales Gonzalez para optar al título de Ingeniero Civil en Computación, Universidad de Chile, DCC. Profesor guía: Gonzalo Navarro. El código se escribe en **C++17** (mínimo), compilado con GCC 12+ sobre Linux (Ubuntu 22.04/24.04).

## Qué hace este proyecto

Construye un **índice comprimido basado en LZ77** para clasificación taxonómica de lecturas de ADN. La pipeline completa es:

1. Concatenar genomas de referencia en orden DFS del árbol filogenético
2. Construir la BWT con **ropebwt3**
3. Computar el parsing **LZ77** del texto (z frases)
4. Construir una **grilla bidimensional** donde cada frase es un punto (x, y)
5. Almacenar la grilla en un **Wavelet Tree** (`wt_int` de sdsl-lite)
6. Integrar **RMQ sucinto** (`rmq_succinct_sct`) en cada nivel del Wavelet Tree
7. Para cada lectura query: extraer MEMs → proyectar a rectángulos en la grilla → RMQ para min/max posición → mapear a genomas → LCA en el árbol filogenético

El **baseline de comparación** es el **sr-index** de Cobas, Gagie y Navarro.

## Decisión de diseño clave: SOLO ocurrencias primarias

Del LZ-index original de Kreft-Navarro, **descartamos deliberadamente** las estructuras para ocurrencias secundarias y reconstrucción de texto. Esto es posible porque para clasificación taxonómica **no necesitamos enumerar TODAS las ocurrencias** de un MEM, sino solo las posiciones extremales (más a la izquierda y más a la derecha) para computar el LCA en el árbol filogenético.

### Qué SÍ construimos (ocurrencias primarias)
- **Parsing LZ77**: las z frases con sus límites
- **Grilla 2D**: los z puntos (x, y) que codifican cruces de límites de frase
- **Wavelet Tree** sobre la grilla para consultas rectangulares
- **RMQ sucinto** integrado en el WT para localización extremal sin enumerar
- **Bitvectors rank/select** para mapear rangos BWT ↔ ejes de la grilla

### Qué NO construimos (ahorro de espacio)
- **Array de profundidad D[1..z]**: usado en el LZ-index original para propagar ocurrencias secundarias vía `prevLess`. NO lo necesitamos → ahorro de z·log(n) bits.
- **Punteros source inversos**: la maquinaria para rastrear recursivamente qué frases copian a qué otras. NO lo necesitamos → ahorro de z·log(n) bits.
- **Permutación P y su inversa**: usada en el LZ-index original para la extracción de texto (reconstruit T[i..j] sin almacenar T). NO lo necesitamos porque la clasificación taxonómica no requiere reconstruir subcadenas del texto.
- **Bitmap S de sources**: marca posiciones source en el texto para la propagación. NO lo necesitamos.
- **LCP de frases**: usado internamente para text extraction. NO lo necesitamos.

### Impacto en espacio
El LZ-index original ocupa ~2|LZ| + o(|LZ|) bits. Nuestra versión "primarias-only" ocupa significativamente menos: solo la grilla (z puntos en WT = z·log(z) bits + overhead), los bitvectors de mapeo (≈z·log(n/z) bits con Elias-Fano), y las RMQ (2z + o(z) bits por nivel). Las estructuras descartadas representaban aproximadamente la mitad del espacio del índice original.

### Justificación teórica
Para clasificación taxonómica, la pregunta no es "¿dónde aparecen TODAS las copias del MEM?" sino "¿cuáles son los genomas extremos donde aparece?". Las ocurrencias primarias ya cubren todos los cruces de límites de frase. Como la concatenación sigue el orden DFS del árbol filogenético, los genomas extremos (min/max posición) determinan el LCA correctamente, sin necesidad de encontrar ocurrencias secundarias intermedias que caerían entre los mismos extremos.

**IMPORTANTE**: Este argumento requiere que la concatenación de genomas siga estrictamente el orden DFS. Si el orden es incorrecto, las posiciones extremales de las primarias podrían no capturar el LCA correcto. Validar esto en los tests.

## Estructura del repositorio

```
LCA-lzindex/
├── CLAUDE.md              ← este archivo
├── CMakeLists.txt         ← build system principal
├── src/
│   ├── main.cpp           ← entry point (placeholder)
│   ├── index.hpp / .cpp   ← LZ77Index: pipeline pública; usa sri::RCSA (PIMPL)
│   ├── lz77/
│   │   ├── phrase.hpp     ← struct Phrase { start_pos, length, next_char }
│   │   ├── parser.hpp/cpp ← parsing LZ77 greedy desde SA (divsufsort + Kasai)
│   │   └── grid.hpp/cpp   ← Grid2D: wt_int + sd_vector fwd/rev + RMQ extremal
│   ├── wavelet/
│   │   ├── wt_rmq_min.hpp/cpp ← legacy: WT pointer-based + RMQ topológico
│   │   └── wm_rmq_min.hpp/cpp ← DEFAULT: wavelet matrix plana + RMQ por nivel
│   ├── mem/
│   │   ├── mem.hpp        ← struct MEM { query_start, ref_start, length }
│   │   └── extractor.hpp/cpp ← extracción de MEMs (sliding window sobre BWT)
│   ├── taxonomy/
│   │   ├── lca.hpp        ← PhyloTree con LCA naive (depth + subida)
│   │   └── classifier.hpp ← Classifier::classify_read (template sobre Index)
│   └── baseline/
│       └── sr_index_locator.hpp/cpp ← wrapper PIMPL del sr-index para comparación
├── tests/                  ← GoogleTest (FetchContent)
│   ├── test_lz77.cpp
│   ├── test_grid.cpp
│   ├── test_wt_rmq_min.cpp
│   ├── test_wm_rmq_min.cpp
│   ├── test_index.cpp
│   ├── test_mem.cpp
│   ├── test_classify.cpp
│   └── test_baseline_sr.cpp
├── tools/                  ← binarios auxiliares (bench y baseline)
│   ├── measure_index.cpp   ← desglose de tamaño por componente
│   ├── bench_special.cpp   ← micro-bench de query_special
│   ├── baseline_lca_sr.cpp ← pipeline LCA end-to-end con sr-index
│   ├── sr_locate.cpp       ← locate sobre sr-index puro
│   ├── makeData.c          ← generador sintético (Mersenne Twister)
│   └── mtwister.{c,h}
├── bench/                  ← (vacío — pendiente)
├── scripts/
│   ├── build_bwt.sh        ← wrapper ropebwt3 build
│   ├── download_corpus.sh  ← descarga Pizza&Chili
│   ├── clean.sh            ← limpia caché sdsl/sr-index/ropebwt3
│   ├── gen_synthetic.py
│   ├── gen_synthetic_large.py
│   ├── run_baseline.sh / run_sr_index.sh / bench_vs_sr.sh
│   └── test_synthetic_lca.py
├── data/                   ← gitignore
│   ├── Escherichia_Coli.fa
│   └── synthetic/
├── docs/
│   ├── glosario.md         ← nomenclatura del proyecto
│   └── *.pdf               ← papers de referencia
├── results/                ← resúmenes de sesión y análisis
└── external/               ← dependencias (clonadas, no submodules)
    ├── sdsl-lite/          ← duscob fork, header-only
    ├── ropebwt3/
    ├── sr-index/
    ├── big-bwt/
    └── uiHRDC/
```

## Dependencias y compilación

### Bibliotecas requeridas

| Biblioteca | Repositorio | Rol | Notas |
|---|---|---|---|
| **sdsl-lite (fork duscob)** | `https://github.com/duscob/sdsl-lite` | Wavelet Tree, RMQ, bitvectors, FM-index | Fork de Dustin Cobas con modificaciones para sr-index. Usar este para mantener equivalencia de estructuras con el baseline. NO usar simongog/sdsl-lite (sin mantenimiento desde 2016) ni xxsds/sdsl-lite (incompatible con sr-index) |
| **ropebwt3** | `https://github.com/lh3/ropebwt3` | Construcción de BWT, extracción de MEMs | Compilar con `make`. Produce archivos .fmr (dinámico) y .fmd (estático) |
| **sr-index** | `https://github.com/duscob/sr-index` | Baseline de comparación | Usa CMake + FetchContent, descarga su propio fork de sdsl-lite |
| **libdivsufsort** | Incluida en sdsl-lite | Construcción de suffix array | Necesaria para `-ldivsufsort -ldivsufsort64` |
| **Google Test** | Via FetchContent | Testing | Solo para tests |
| **Google Benchmark** | Via FetchContent | Benchmarking | Solo para benchmarks |

### Compilación

```bash
# Setup inicial
# NOTA: el fork duscob/sdsl-lite ES header-only (divsufsort incluido inline,
# sin libsdsl.a). Basta con tener los headers visibles en $HOME/.local/include.
git clone --recursive https://github.com/duscob/sdsl-lite.git
cd sdsl-lite && ./install.sh $HOME/.local

# Build del proyecto (CMakeLists usa SDSL_INCLUDE = $HOME/.local/include
# y SR_INDEX_INCLUDE = external/sr-index/include directamente).
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel $(nproc)

# Ejecutar tests (GoogleTest se descarga vía FetchContent al primer build)
cd build && ctest --output-on-failure

# Binarios auxiliares actuales (no hay todavía bench/ propio)
./build/measure_index <texto>          # desglose de espacio por componente
./build/baseline_lca_sr <args>         # pipeline LCA con sr-index
./build/sr_locate <args>               # locate sobre sr-index puro
```

### Flags de compilación obligatorios

Siempre compilar con:
- `-O3 -DNDEBUG` (rendimiento)
- `-march=native -msse4.2 -mbmi -mbmi2` (instrucciones SIMD requeridas por sdsl-lite)
- `-std=c++17` mínimo
- `-Wall -Wextra -pedantic` (solo en desarrollo)
- **NO usar** `-fsanitize=address` con benchmarks (distorsiona tiempos)

## Estructuras de datos clave — API de sdsl-lite

### Wavelet Tree con consultas 2D

```cpp
#include <sdsl/wavelet_trees.hpp>
using namespace sdsl;

// Construcción: R[i] = coordenada Y del punto con coord X = i
int_vector<> R = {/* z valores */};
wt_int<> wt;
construct_im(wt, R);

// Consulta rectangular [x1,x2] × [y1,y2]
auto result = wt.range_search_2d(x1, x2, y1, y2);
// result.first  = número de puntos encontrados
// result.second = vector<pair<size_t, size_t>> con (posición_x, valor_y)

// CUIDADO: Issue #143 reporta bug en range_search_2d para casos borde.
// Siempre verificar contra brute-force en tests.
```

### RMQ sucinto

```cpp
#include <sdsl/rmq_support.hpp>

int_vector<> values = {/* posiciones de texto asociadas a puntos */};
rmq_succinct_sct<> rmq(&values);  // 2n+o(n) bits
auto min_idx = rmq(i, j);         // O(1), retorna ÍNDICE del mínimo en [i..j]
// Luego: values[min_idx] es el valor mínimo
// NOTA: después de construir rmq, values puede liberarse si no se necesita más
```

### Bitvectors sparse con rank/select

```cpp
#include <sdsl/sd_vector.hpp>

bit_vector bv(n, 0);
// Marcar límites de frase
for (auto& phrase : phrases) bv[phrase.start] = 1;

sd_vector<> sdb(bv);              // Elias-Fano, óptimo para sparse
sd_vector<>::rank_1_type rank_bv(&sdb);
sd_vector<>::select_1_type select_bv(&sdb);

size_t phrase_id = rank_bv(pos);         // ¿en qué frase cae pos?
size_t phrase_start = select_bv(id + 1); // ¿dónde empieza la frase id? (1-indexed)
```

### CSA primaria: `sri::RCSA<>` (no `csa_wt`)

A partir del commit `8f946ff` el índice usa `sri::RCSA<>` del sr-index como
CSA forward/reverso, envuelto en `CountOnlyRCSA` (subclase definida en
`src/index.cpp`) que **solo** carga/serializa `Alphabet + PsiCoreRLE` y stubs
out de Locate. Motivación: en E. coli (108 MB) bajó la CSA de 9.20 → 3.76 bpc.
Locate vive afuera de la CSA (via Grid2D + WmMinRmq), por eso podemos descartar
samples/marks de la RCSA.

Para verificación rápida contra un FM-index estándar de sdsl-lite:

```cpp
#include <sdsl/suffix_arrays.hpp>
sdsl::csa_wt<sdsl::wt_huff<>, 32, 32> fm;
sdsl::construct(fm, "text_file.txt", 1);
auto cnt = sdsl::count(fm, "ACGT");
auto occ = sdsl::locate(fm, "ACGT");  // ground truth en tests
```

### Wavelet Matrix + RMQ por nivel (`WmMinRmq`, default)

Reemplazo de `WtMinRmq` (pointer-based, O(z) instancias sdsl, ~63 bpc en E. coli)
por una wavelet matrix plana: `sdsl::wm_int<>` + un `rmq_succinct_sct` por nivel.
Espacio O(z log σ) bits con O(log σ) instancias sdsl en lugar de O(z).

```cpp
#include "wavelet/wm_rmq_min.hpp"
lz77tax::WmMinRmq wm;
wm.build(y_values, text_pos, sigma);     // y_values en orden X
auto r = wm.range_argmin_2d(x_lo, x_hi, y_lo, y_hi, text_pos);
// r.argmin_global = índice en [0,n) con text_pos[idx] mínimo
```

`Grid2D::build()` toma `RmqVariant` (`Wm` default, `Wt` legacy para benchmarks
comparativos). `WtMinRmq` se conserva solo para reproducir mediciones previas.

## ropebwt3 — comandos clave

```bash
# Construir BWT incremental (formato dinámico .fmr)
ropebwt3 build -t$(nproc) -bo index.fmr genome1.fa genome2.fa

# Convertir a formato estático (más rápido para queries)
ropebwt3 build -i index.fmr -do index.fmd

# Extraer SMEMs con longitud mínima 31
ropebwt3 mem -t4 -l31 index.fmd reads.fa > mems.bed

# Exportar BWT a texto plano (útil para debugging)
ropebwt3 get -a index.fmd > bwt.txt

# NOTAS:
# - Por defecto indexa forward + reverse complement. Usar -R para solo forward.
# - Formato .fmd es memory-mappable (rápido para cargar).
# - La salida de `mem` es BED-like: query_name, start, end, num_hits.
# - Con -p se obtienen posiciones en el texto de referencia.
```

## sr-index — comandos clave

```bash
# Construir índice con variantes s=4..64
./build/tools/bm_construct_ri \
  --data=data/escherichia_coli.txt \
  --min_s=4 --max_s=64 \
  --benchmark_counters_tabular=true

# Ejecutar benchmark de locate
./build/tools/bm_locate_ri \
  --data=data/escherichia_coli.txt \
  --patterns=data/patterns.txt \
  --min_s=4 --max_s=64

# NOTAS:
# - Variante recomendada: SrIndexValidArea (mejor espacio/tiempo)
# - Parámetro s: mayor s = menos espacio, más lento locate (factor s)
# - Para comparación justa: ambos índices deben usar la MISMA BWT
```

## Parsing LZ77

### Desde Suffix Array (enfoque simple, recomendado para empezar)

```cpp
// Parsing LZ77 greedy left-to-right
// Requiere: SA, ISA, LCP arrays (construidos por sdsl-lite/libdivsufsort)
//
// DECISIÓN DE DISEÑO: almacenamos start_pos (posición en T donde inicia la frase)
// y length, pero NO almacenamos source (posición de la copia anterior).
// source solo se necesita para propagar ocurrencias secundarias, que descartamos.

struct Phrase {
    size_t start_pos;    // posición en T donde comienza esta frase
    size_t length;       // largo de la frase (0 si es literal)
    uint8_t next_char;   // carácter explícito al final
    // NO hay: size_t source — no propagamos secundarias
};

vector<Phrase> parse_lz77(const string& text) {
    vector<Phrase> phrases;
    size_t i = 0;
    while (i < text.size()) {
        auto [source, len] = longest_previous_factor(i, SA, ISA, LCP);
        if (len == 0) {
            phrases.push_back({i, 0, text[i]});
            i++;
        } else {
            phrases.push_back({i, len, text[i + len]});
            i += len + 1;
        }
        // NOTA: source se usa solo para calcular len, luego se descarta.
        // No lo almacenamos en la Phrase.
    }
    return phrases;
}
// VERIFICACIÓN: comparar z obtenido con estadísticas de Pizza&Chili
```

### Desde RLBWT (enfoque avanzado, espacio O(r))

Algoritmo de Policriti y Prezza (Algorithmica 2018). Tiempo O(n log r), espacio O(r) words. Más complejo de implementar pero no requiere materializar el texto completo. Referencia de implementación: `https://github.com/nicolaprezza/lz-rlbwt`.

## Inventario de estructuras: LZ-index original vs nuestra versión

| Estructura | LZ-index original (Kreft-Navarro) | Nuestra versión | Espacio ahorrado |
|---|---|---|---|
| Frases LZ77 (límites) | ✅ | ✅ | — |
| Grilla 2D de puntos | ✅ | ✅ | — |
| Wavelet Tree sobre grilla | ✅ | ✅ | — |
| Bitvectors rank/select (BWT↔frases) | ✅ | ✅ | — |
| RMQ sucinto en WT/WM | ❌ (no en original) | ✅ `WmMinRmq` por defecto (legacy `WtMinRmq` opcional) | +O(z log σ) bits |
| CSA forward/reverso | csa_wt completa | `sri::RCSA<>` con `CountOnlyRCSA` (sin samples de locate) | E. coli: 9.20 → 3.76 bpc |
| **Array de profundidad D[1..z]** | ✅ (para prevLess) | ❌ **DESCARTADO** | z·log(n) bits |
| **Punteros source inversos** | ✅ (para secundarias) | ❌ **DESCARTADO** | z·log(n) bits |
| **Permutación P + P⁻¹** | ✅ (para extract) | ❌ **DESCARTADO** | 2z·log(z) bits |
| **Bitmap S de sources** | ✅ (para secundarias) | ❌ **DESCARTADO** | n bits (!!!) |
| **LCP de frases** | ✅ (para extract) | ❌ **DESCARTADO** | z·log(n) bits |

**Resultado neto**: descartamos ~4z·log(n) + 2z·log(z) + n bits del original, y añadimos solo ~2z bits de RMQ. Para colecciones grandes (n >> z), el ahorro es masivo — especialmente el bitmap S que es O(n).

## Grilla bidimensional LZ77 — construcción paso a paso

1. **Obtener z frases** LZ77 del texto T[0..n-1]
2. **Marcar límites** de frase en bitvector B[0..n-1] (B[pos]=1 si pos es inicio de frase)
3. **Construir SA del texto** (puede reutilizarse el de sdsl-lite)
4. **Para cada frase k**: su coordenada Y = rank en SA de la posición donde inicia la frase
5. **Construir SA del texto reverso** T^R
6. **Para cada frase k**: su coordenada X = rank en SA reverso de la posición donde termina la frase
7. **Generar secuencia R[1..z]**: R[i] = coordenada Y del punto cuya coordenada X = i
8. **Construir `wt_int` sobre R**

### Búsqueda de un patrón P[1..m] — SOLO ocurrencias primarias

Para cada partición P[1..i] · P[i+1..m] (con i = 1..m-1):
- Buscar P[i+1..m] en la BWT directa → intervalo [sp, ep] → proyectar al eje Y
- Buscar P[1..i]^R en la BWT reversa → intervalo [sp', ep'] → proyectar al eje X
- Consulta rectangular `wt.range_search_2d(sp', ep', sp, ep)` → ocurrencias primarias
- **NO propagamos ocurrencias secundarias** — para clasificación taxonómica solo
  necesitamos las posiciones extremales, que se obtienen directamente con RMQ

### Caso especial: patrones end-aligned dentro de una frase

`Grid2D::query_special(sp_rev, ep_rev, plen)` (commit `316b501`) captura el caso
en que el patrón termina exactamente al final de una frase k y cabe dentro de ella
(`phrase_total_len_[k] >= plen`). Esto mitiga parcialmente el "caso patológico"
de MEMs internos a una frase (ver sección "Trampa del diseño primarias-only").

### Localización extremal con RMQ

En vez de listar todos los puntos del rectángulo:
```cpp
// Obtener la posición mínima (más a la izquierda) en el texto
auto min_idx = rmq_at_level(sp', ep', sp, ep);  // custom, integrado en WT
auto min_text_pos = text_positions[min_idx];

// Mapear a genoma
auto genome_id = rank_genome_bv(min_text_pos);

// Obtener la posición máxima (más a la derecha) análogamente
// Computar LCA(genome_min, genome_max) en el árbol filogenético
```

## Datasets de desarrollo y prueba

### Datasets actualmente en el repo

| Dataset | Tamaño | z (LZ77) | Estado | Uso |
|---------|--------|----------|--------|-----|
| Escherichia_Coli.fa | ~108 MB | ~1.75M | descargado | Dev + tests + benchmark principal |
| synthetic/z_10k | 1 MB | ~6.9k | generado | Verificación overhead WtMinRmq |
| synthetic/z_50k | 5 MB | ~27k | generado | Verificación overhead WtMinRmq |
| synthetic/z_100k | 10 MB | ~49k | generado | Verificación overhead WtMinRmq |
| synthetic (LCA toy) | 85 B | trivial | generado | Validación end-to-end del pipeline LCA |

Pendientes de bajar de Pizza & Chili (si se quiere ampliar la evaluación):
`cere` (~461 MB), `para` (~429 MB), `influenza` (~154 MB).

Descarga: `https://pizzachili.dcc.uchile.cl/repcorpus.html`

```bash
# Script de descarga (scripts/download_corpus.sh)
wget -P data/ https://pizzachili.dcc.uchile.cl/repcorpus/real/Escherichia_Coli.gz
wget -P data/ https://pizzachili.dcc.uchile.cl/repcorpus/real/cere.gz
wget -P data/ https://pizzachili.dcc.uchile.cl/repcorpus/real/para.gz
wget -P data/ https://pizzachili.dcc.uchile.cl/repcorpus/real/influenza.gz
cd data && gunzip *.gz
```

### Generación de patrones de prueba

```bash
# Extraer 100k patrones aleatorios de longitud m desde el texto
python3 -c "
import random, sys
text = open(sys.argv[1], 'rb').read()
m = int(sys.argv[2])
n = int(sys.argv[3])
for _ in range(n):
    pos = random.randint(0, len(text) - m)
    sys.stdout.buffer.write(text[pos:pos+m] + b'\n')
" data/escherichia_coli 50 100000 > data/patterns_50.txt
```

## Protocolo de benchmarking

### Métricas

1. **Espacio**: bits por carácter (bpc) = (tamaño_índice_bytes × 8) / n
2. **Tiempo de construcción**: wall-clock seconds
3. **Memoria peak de construcción**: MaxRSS vía `/usr/bin/time -v`
4. **Tiempo de count**: microsegundos por query
5. **Tiempo de locate**: microsegundos por ocurrencia

### Ejecución justa

```bash
# Antes de cada medición
sync; echo 3 | sudo tee /proc/sys/vm/drop_caches > /dev/null

# Medir con /usr/bin/time -v para MaxRSS
/usr/bin/time -v ./build/bench/bench_query 2> timing.txt

# Ambos índices deben:
# - Usar la MISMA BWT (construida con ropebwt3)
# - Buscar los MISMOS patrones
# - Compilarse con los MISMOS flags (-O3 -DNDEBUG -march=native)
```

## Testing — reglas obligatorias

1. **Todo componente tiene test unitario** antes de integrarse
2. **Verificar z contra Pizza&Chili** al implementar el parser LZ77
3. **Verificar count/locate contra FM-index de sdsl-lite** (ground truth)
4. **Mínimo 10,000 patrones aleatorios** por dataset en tests de correctitud
5. **Tests de regresión** para bugs encontrados (especialmente off-by-one en SA/BWT)
6. **No optimizar sin tests verdes** — correctitud antes que rendimiento

```bash
# Ejecutar todos los tests
cd build && ctest --output-on-failure -j$(nproc)

# Ejecutar un test específico con output detallado
./build/tests/test_lz77 --gtest_filter="*ParseEcoli*"
```

## Errores comunes y cómo evitarlos

### Off-by-one en SA/BWT
- sdsl-lite usa **0-indexed** SA por defecto
- ropebwt3 puede usar **1-indexed** en su formato interno
- `select` en sdsl-lite es **1-indexed** (select(1) = primera posición)
- `rank` en sdsl-lite es **0-indexed** (rank(i) = bits set en [0..i))
- **Siempre** escribir un test que verifique `T[SA[i]..] < T[SA[i+1]..]` para todo i

### Formato .fmd de ropebwt3
- Es un B+-tree serializado con run-length encoding tipo UTF-8
- **No parsear directamente** a menos que sea necesario
- Usar `ropebwt3 get -a` para exportar a texto plano como paso intermedio
- Si se necesita parseo directo, estudiar `rb3_fmt.c` en el repo de ropebwt3

### range_search_2d de sdsl-lite
- Bug reportado en Issue #143: puede retornar puntos fuera del rango en casos borde
- **Siempre** verificar resultados contra búsqueda brute-force
- Considerar implementar wrapper que filtre resultados espurios

### Serialización de sdsl
- Todas las estructuras sdsl-lite soportan `store_to_file()` y `load_from_file()`
- Usar para persistir el índice construido y no reconstruir cada vez
- `size_in_bytes()` para medir espacio

### Trampa del diseño primarias-only
- Al descartar ocurrencias secundarias, asumimos que las posiciones extremales de las
  primarias son suficientes para determinar el LCA correcto. Esto es válido SOLO si:
  1. La concatenación sigue orden DFS estricto del árbol filogenético
  2. Para el MEM buscado, existen ocurrencias primarias en los genomas extremos
- **Caso patológico**: si un MEM aparece SOLO dentro de una frase LZ77 (sin cruzar
  límites), no tiene ocurrencia primaria → sería invisible para nuestro índice.
  Esto ocurre cuando el MEM es más corto que la frase que lo contiene.
  **Mitigación parcial implementada**: `Grid2D::query_special` (commit `316b501`)
  captura patrones end-aligned que caben dentro de su frase contenedora.
  **Mitigación heurística**: MEMs largos (≥ longitud mínima de frase) siempre cruzan
  al menos un límite. Para MEMs cortos, la evidencia taxonómica es débil de todos modos.
- **Test crítico**: verificar que para MEMs de longitud ≥ L_min, el LCA obtenido
  solo con primarias coincide con el LCA obtenido enumerando TODAS las ocurrencias
  (usando FM-index como ground truth, o `SrIndexLocator` para corpus reales)

## Cronograma de 15 semanas

| Semana | Fase | Entregable | Estado |
|--------|------|-----------|--------|
| 1 | Setup | Entorno compilando, corpus descargado, BWT construida | ✅ |
| 2 | Teoría | Papers leídos, caps 1-2 escritos, sr-index operativo | ✅ (sr-index ok; escritura pendiente) |
| 3 | LZ77 | Parser LZ77 desde SA implementado y verificado | ✅ |
| 4 | LZ77 | Grilla 2D construida, bitvectors rank/select | ✅ |
| 5 | Grilla | Búsqueda brute-force sobre grilla, cap 3 escrito | ✅ (brute-force; escritura pendiente) |
| 6 | WT | `wt_int` integrado, count query funcional | ✅ |
| 7 | RMQ | RMQ en WT implementado, locate query funcional | ✅ (WtMinRmq + WmMinRmq) |
| 8 | Verif | Correctitud verificada, benchmarks preliminares | 🟡 (solo E. coli + synthetic) |
| 9 | MEMs | Extracción de MEMs conectada al índice | ✅ |
| 10 | Tax | Clasificación taxonómica (LCA) implementada | ✅ |
| 11 | Opt | Hot paths optimizados, benchmark completo | 🟡 (RCSA + wavelet matrix integrados) |
| 12 | Eval | sr-index benchmarked, comparación cuantitativa | 🟡 (`baseline_lca_sr` operativo) |
| 13 | Eval | Análisis de resultados, caps 4-5 escritos | ⏳ |
| 14 | Doc | Borrador completo, revisión con profesor | ⏳ |
| 15 | Doc | Memoria entregada, defensa preparada | ⏳ |

## Qué cortar si el tiempo se agota (en orden de prioridad)

1. Clasificación taxonómica (semanas 10-11) → demostrar solo MEMs
2. Extracción de MEMs (semana 9) → enfocarse solo en count/locate del LZ77-index
3. Comparación con sr-index → comparar solo con FM-index estándar de sdsl-lite
4. Benchmarks en todos los datasets → evaluar solo en Escherichia_coli y cere
5. **NUNCA CORTAR**: verificación de correctitud, al menos 1 benchmark comparativo, escritura

## Convenciones de código

- **Namespaces**: `lz77tax::` para todo el código del proyecto
- **Nombres**: `snake_case` para funciones y variables, `PascalCase` para tipos/clases
- **Headers**: `.hpp` para headers, `.cpp` para implementación
- **Includes**: `#pragma once` en todos los headers
- **Documentación**: Doxygen comments en interfaces públicas
- **Formato**: aún sin `.clang-format` en raíz — mantener estilo a mano (4 espacios,
  llaves K&R) hasta que se agregue uno

## Papers de referencia rápida

| Cita | Paper | Componente |
|------|-------|-----------|
| [KN13] | Kreft & Navarro, TCS 2013 | LZ77 self-index, grilla 2D, primarias/secundarias |
| [CGN24] | Cobas, Gagie & Navarro, 2024 | sr-index (baseline) |
| [D+24] | Draesslerová et al., SEA 2024 | MEMs para clasificación taxonómica |
| [Nav14] | Navarro, JDA 2014 | Wavelet Trees |
| [FN16] | Ferrada & Navarro, DCC 2016 | RMQ sucinto |
| [Nav20] | Navarro, ACM CS 2022 | Survey indexación repetitiva |
| [PP18] | Policriti & Prezza, Algorithmica 2018 | LZ77 desde RLBWT |
