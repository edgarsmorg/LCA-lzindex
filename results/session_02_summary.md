# Sesión 02 — Resumen y estado con certeza

**Fecha:** 2026-04-08
**Semanas cubiertas:** 2 (revisión y corrección) y 3 del cronograma

---

## Contexto: qué encontramos al inicio

La sesión 01 había documentado resultados como "pipeline LCA validado end-to-end",
pero al revisar el código se constató que:

- El locate de patrones era **brute-force Python** (`str.find()`), no un índice
- `sr_locate.cpp` había sido escrito por Codex 5.3 pero **nunca compilado**
- `scripts/test_synthetic_lca.py` tenía un **syntax error** (variable `SR_LOCATE`
  partida en dos líneas por Codex)
- El "hallazgo clave" de que extremal ≡ completo estaba validado en 5 casos
  triviales con búsqueda lineal — sigue siendo correcto conceptualmente, pero
  la evidencia era débil

---

## Lo que se hizo esta sesión

### Fase 1 — Reparar y completar el baseline sr-index

| Tarea | Estado |
|-------|--------|
| Fix syntax error en `test_synthetic_lca.py` | ✅ |
| Diagnosticar: sdsl duscob fork es **header-only** (no hay `libsdsl.a`) | ✅ |
| Compilar `sr_locate` con linking correcto (solo gflags + pthread) | ✅ |
| Validar `sr_locate` contra brute-force en datos sintéticos | ✅ |
| Probar baseline sobre E. coli real (108 MB) | ✅ |
| `scripts/run_baseline.sh` — pipeline end-to-end automatizado | ✅ |
| `results/baseline_usage.md` — guía completa para nuevos lectores | ✅ |

### Fase 2 — LZ77 Parser (Semana 3 del roadmap)

| Tarea | Estado |
|-------|--------|
| `src/lz77/parser.cpp` — implementación completa | ✅ |
| GTest via FetchContent + `tests/test_lz77.cpp` | ✅ |
| 10/10 tests unitarios pasan | ✅ |
| Verified z=6 en "abracadabra" (canónico) | ✅ |
| Verified z=14 en texto sintético (n=85) | ✅ |

---

## Certezas técnicas establecidas esta sesión

### Sobre el baseline sr-index

**Compilación correcta de sr_locate:**
El fork duscob/sdsl-lite es completamente header-only. `divsufsort` está inline en
`sdsl/divsufsort.hpp` (2665 líneas). No existe ningún `libsdsl.a`, `libdivsufsort.a`
ni `libdivsufsort64.a`. El linking solo requiere gflags y pthread.

**Resultados en E. coli (n=112,689,517):**

| Métrica | Valor |
|---------|-------|
| n (chars) | 112,689,517 |
| r (BWT runs, R-Index) | 15,044,488 |
| r' (muestras SR-Index s=16) | 1,928,549 |
| Construcción R-Index | ~73 seg |
| Construcción SR-Index s=16 | ~10 seg (incremental) |
| Locate 1000 queries longitud 50 | 220 ms total (~0.22 ms/query) |
| Ocurrencias promedio por query | 9.4 |

El ratio r/n ≈ 13.3% — E. coli tiene repetitividad moderada. Para `cere` y `para`
se esperan ratios de 1-3%, donde la compresión LZ-based es más efectiva.

### Sobre el LZ77 Parser

**Algoritmo implementado:** Parsing greedy left-to-right desde SA+LCP.
- SA construido con `divsufsort` (header en sdsl)
- ISA por inversión directa del SA
- LCP con el algoritmo de Kasai
- Búsqueda del factor previo más largo: escaneo izquierdo/derecho desde ISA[i],
  manteniendo mínimo corriente de LCP; para cuando min cae a 0

**Propiedades correctas verificadas en tests:**
1. Las frases parten el texto exactamente (sin gaps ni overlaps en la cobertura)
2. `phrases[0].start_pos == 0`, avance coherente entre frases consecutivas
3. `next_char == text[start_pos + length]` para cada frase
4. Toda copia tiene fuente válida en el texto previo (incluyendo copias solapadas,
   que son válidas en LZ77 estándar)

**Copias solapadas:** son parte del estándar LZ77. Ejemplo: en "aaaa", la segunda
frase copia 2 caracteres desde la posición 0 aunque la copia se solapa con el
destino. Esto es correcto y esperado.

**Resultados en textos conocidos:**

| Texto | n | z | Nota |
|-------|---|---|------|
| "abracadabra" | 11 | 6 | Canónico, matches literatura |
| "ACGTACGTACGT" | 12 | 5 | 4 literales + 1 copia |
| Texto sintético | 85 | 14 | 4 genomas, separadores '#' |

**Frases del texto sintético (n=85, z=14):**
```
[0]   len=0  'A'         → literal A
[1]   len=4  'C'         → copia 4 chars + 'C'  (i=1..5 cubierto)
[6]   len=0  'G'         → literal G
[7]   len=0  'T'         → literal T
[8]   len=12 '#'         → copia 12 chars + '#'
[21]  len=12 'G'         → copia 12 chars + 'G'
[34]  len=7  '#'         → copia 7 chars + '#'
[42]  len=8  'T'         → copia 8 chars + 'T'
[51]  len=5  'G'         → copia 5 chars + 'G'
[57]  len=1  'A'         → copia 1 char + 'A'
[59]  len=4  '#'         → copia 4 chars + '#'
[64]  len=8  'T'         → copia 8 chars + 'T'
[73]  len=3  'C'         → copia 3 chars + 'C'
[77]  len=7  '#'         → copia 7 chars + '#'
```

---

## Lo que NO está verificado todavía

| Ítem | Pendiente |
|------|-----------|
| Valor de z en E. coli (108 MB) | No medido aún (parser es O(n²) worst case, podría ser lento) |
| z de E. coli vs Pizza&Chili reference | Sin confrontar con valor publicado |
| Parser sobre textos > 2GB | Limitado a int32_t en SA (divsufsort por defecto) |
| Rendimiento del parser a escala | No benchmarkeado |
| Validación de la propiedad extremal con índice real | Solo verificada con brute-force Python sobre 5 queries |

---

## Estado actual del repositorio

```
12 commits en main

src/
├── lz77/
│   ├── phrase.hpp           — struct Phrase{start_pos, length, next_char} ✓
│   ├── parser.hpp           — interfaz LZ77Parser ✓
│   └── parser.cpp           — IMPLEMENTADO: SA+ISA+LCP+greedy ✓
│   └── grid.{hpp,cpp}       — stub (pendiente sesión 03)
├── mem/                     — stubs (pendiente)
├── taxonomy/
│   ├── lca.hpp              — PhyloTree + LCA naive ✓
│   └── classifier.hpp       — Classifier + classify_extremal() ✓
└── bwt/, wavelet/           — placeholders vacíos

tests/
└── test_lz77.cpp            — 10/10 pasan ✓

scripts/
├── run_baseline.sh          — pipeline end-to-end sr-index ✓
└── test_synthetic_lca.py    — 5/5 pasan (brute-force) ✓

tools/
├── sr_locate.cpp            — locate en batch sobre sr-index ✓
└── sr_index_template.cpp    — template de construcción

results/
├── baseline_usage.md        — guía de uso del baseline ✓
├── session_01_summary.md    — resumen semanas 1-2 (con advertencias)
├── session_02_summary.md    — este archivo
└── pipeline_sr_index_lca.md — arquitectura del baseline
```

---

## Deuda técnica identificada

- [ ] Medir z del parser en E. coli y confrontar con Pizza&Chili
- [ ] Benchmarkear tiempo del parser en textos > 1 MB
- [ ] Considerar optimización con PSV/NSV o RMQ para O(n) garantizado
- [ ] Validar propiedad extremal con sr_locate real (no solo brute-force Python)
- [ ] `.clang-format` con estilo Google aún pendiente
- [ ] `docs/memoria.tex` estructura inicial pendiente

---

## Roadmap — Próxima sesión

### Sesión 03 — Semana 4: Grilla 2D + Bitvectors

**Objetivo:** Construir los z puntos (x,y) a partir de SA directo y SA reverso,
con bitvectors `sd_vector` para el mapeo BWT↔frases.

**Tareas:**
1. Implementar `Grid2D::build()` en `src/lz77/grid.cpp`:
   - SA del texto directo → coordenadas Y (rank del inicio de cada frase)
   - SA del texto reverso T^R → coordenadas X (rank del final de cada frase)
   - Secuencia `R[i]` = coordenada Y del punto con X=i
2. Bitvectors `sd_vector` para límites de frases (rank/select)
3. Búsqueda brute-force `range_search_2d` como ground truth
4. Test: `Grid2D::point_count() == z`

**Archivos a tocar:**
- `src/lz77/grid.hpp` — completar interfaz
- `src/lz77/grid.cpp` — implementación
- `tests/test_grid.cpp` — nuevo

---

## Comandos útiles

```bash
# Compilar y correr tests
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build --parallel
cd build && ctest --output-on-failure

# Verificar pipeline sintético
python3 scripts/test_synthetic_lca.py

# Correr baseline completo
./scripts/run_baseline.sh data/synthetic/reference.fa data/test_seq.fa results/test 8 16

# Locate en batch sobre E. coli (ya construido)
./build/sr_locate --data_dir=results/baseline_ecoli --data_name=Escherichia_Coli.txt \
    --patterns=data/patterns_ecoli_50.txt > /tmp/locs.txt
```
