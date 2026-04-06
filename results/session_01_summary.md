# Sesión 01 — Resumen y Roadmap

**Fecha:** 2026-04-06  
**Semanas cubiertas:** 1 y 2 del cronograma de 15 semanas

---

## Lo que hicimos

### Semana 1 — Setup del entorno

| Tarea | Estado |
|-------|--------|
| Estructura de directorios (`src/`, `tests/`, `bench/`, `scripts/`, `data/`, `docs/`, `results/`) | ✅ |
| Headers de todos los módulos (`phrase.hpp`, `parser.hpp`, `grid.hpp`, `mem.hpp`, `extractor.hpp`, `index.hpp`) | ✅ |
| Stubs de implementación que compilan (`parser.cpp`, `grid.cpp`, `index.cpp`) | ✅ |
| CMakeLists.txt con flags correctos (`-O3 -march=native -msse4.2 -mbmi -mbmi2 -std=c++17`) | ✅ |
| ropebwt3 compilado y verificado (build BWT sobre datos de prueba) | ✅ |
| Scripts `build_bwt.sh`, `download_corpus.sh` | ✅ |
| `.gitignore` completo para C++ + sdsl + corpus | ✅ |
| Permisos de Claude Code configurados (`defaultMode: auto`) | ✅ |

### Semana 2 — Teoría y baseline

| Tarea | Estado |
|-------|--------|
| sdsl-lite (fork duscob) compilado e instalado en `~/.local` | ✅ |
| sr-index compilado con todos sus benchmarks (`bm_construct_ri`, `bm_locate_ri`, `bm_count_ri`) | ✅ |
| sr-index verificado sobre texto de prueba (n=85, r=32 runs) | ✅ |
| `PhyloTree` con LCA naive implementado (`src/taxonomy/lca.hpp`) | ✅ |
| `Classifier` con `classify_extremal()` implementado (`src/taxonomy/classifier.hpp`) | ✅ |
| Dataset sintético generado: 4 genomas, árbol de profundidad 2, 5 queries | ✅ |
| Pipeline LCA end-to-end validado: **5/5 tests pasan** | ✅ |
| Bosquejo de arquitectura con tamaños reales (`results/pipeline_sr_index_lca.md`) | ✅ |
| `scripts/clean.sh` para limpiar caché sdsl/sr-index | ✅ |

### Hallazgo clave de la sesión

> **El LCA extremal (solo `pos_min` y `pos_max`) produce resultados idénticos al LCA completo (todas las posiciones)** — validado empíricamente en los 5 casos sintéticos. Esto confirma la premisa central del diseño primarias-only del LZ77-index.

---

## Estado actual del repositorio

```
7 commits en main
src/
├── index.{hpp,cpp}          — interfaz pública (stub)
├── lz77/
│   ├── phrase.hpp           — struct Phrase (sin source, primarias-only)
│   ├── parser.{hpp,cpp}     — LZ77Parser (stub, pendiente implementar)
│   └── grid.{hpp,cpp}       — Grid2D (stub, pendiente implementar)
├── mem/
│   ├── mem.hpp              — struct MEM
│   └── extractor.hpp        — MEMExtractor (interfaz)
└── taxonomy/
    ├── lca.hpp              — PhyloTree + LCA naive ✓
    └── classifier.hpp       — Classifier + classify_extremal() ✓

external/ (compilados)
├── ropebwt3/                — ropebwt3 ejecutable listo
├── sdsl-lite/               — instalado en ~/.local
└── sr-index/                — bm_construct_ri, bm_locate_ri, bm_count_ri listos

data/synthetic/              — árbol + genomas + queries de referencia
results/
├── pipeline_sr_index_lca.md — arquitectura del baseline documentada
└── session_01_summary.md    — este archivo
```

---

## Roadmap — Próximas sesiones

### Sesión 02 — Semana 3: LZ77 Parser

**Objetivo:** Parser LZ77 implementado y verificado contra Pizza&Chili.

**Tareas:**
1. Implementar `LZ77Parser::parse()` usando SA + LCP desde sdsl-lite
   - Construir SA con `sdsl::construct`
   - Búsqueda del factor más largo previo con PSV/NSV sobre LCP
   - Output: vector de `Phrase{start_pos, length, next_char}`
2. Test unitario en `tests/test_lz77.cpp`:
   - Verificar que la concatenación de frases reconstruye el texto
   - Verificar z (número de frases) contra valores conocidos
3. Verificar z sobre el texto sintético (n=85) manualmente
4. Descargar corpus real (Escherichia_coli) y verificar z contra Pizza&Chili

**Archivos a tocar:**
- `src/lz77/parser.cpp` — implementación completa
- `tests/test_lz77.cpp` — nuevo archivo de tests
- `CMakeLists.txt` — agregar GTest vía FetchContent

---

### Sesión 03 — Semana 4: Grilla 2D + Bitvectors

**Objetivo:** `Grid2D` construida y consultas rectangulares brute-force funcionando.

**Tareas:**
1. Implementar `Grid2D::build()`:
   - SA del texto directo → coordenadas Y
   - SA del texto reverso → coordenadas X
   - Secuencia R[i] = Y del punto con X=i
2. Bitvectors `sd_vector` para mapeo BWT ↔ frases
3. Brute-force `range_search_2d` para verificación
4. Test: `Grid2D::point_count() == z`

---

### Sesión 04 — Semanas 5-6: Wavelet Tree

**Objetivo:** `wt_int` de sdsl-lite integrado, `count` query funcional.

**Tareas:**
1. Construir `wt_int` sobre la secuencia R
2. Consultas rectangulares `wt.range_search_2d(x1,x2,y1,y2)`
3. Workaround bug Issue #143 (filtrar resultados fuera de rango)
4. Verificar count contra FM-index de sdsl (ground truth)

---

### Sesión 05 — Semana 7: RMQ + locate extremal

**Objetivo:** Localización extremal sin enumerar todas las ocurrencias.

**Tareas:**
1. Integrar `rmq_succinct_sct` sobre posiciones de texto en la grilla
2. Implementar `locate_extremal(pattern)` → `(pos_min, pos_max)`
3. Conectar con `Classifier::classify_extremal()`
4. Benchmark comparativo: enumerate-all vs extremal sobre datos sintéticos

---

### Deuda técnica pendiente

- [ ] `.clang-format` con estilo Google
- [ ] Corpus real descargado (URLs de Pizza&Chili a verificar)
- [ ] `scripts/run_sr_index.sh` — ejecutar desde subdirectorio para no contaminar cwd con .sdsl
- [ ] LaTeX `docs/memoria.tex` — estructura inicial caps 1-2

---

## Comandos útiles para la próxima sesión

```bash
# Limpiar caché antes de empezar
./scripts/clean.sh

# Compilar el proyecto
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build --parallel

# Verificar pipeline sintético sigue verde
python3 scripts/test_synthetic_lca.py

# Correr sr-index sobre datos sintéticos
./scripts/run_sr_index.sh data/synthetic/reference.txt 4 16
```
