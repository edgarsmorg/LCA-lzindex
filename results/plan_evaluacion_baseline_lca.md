# Plan: Evaluación del baseline sr-index — ¿resuelve el problema y es equivalente al LZ77-Index?

## Contexto

El proyecto compara un LZ77-Index "primarias-only" contra el sr-index de Cobas-Gagie-Navarro como baseline. Hoy, la comparación está **incompleta**:

- `tools/baseline_lca_sr.cpp` solo mide tiempo de `sr_loc.locate_extremal(pattern)` sobre patrones crudos. **No construye árbol filogenético, no usa MEMExtractor, no calcula LCA**: el "baseline" no resuelve hoy el problema de clasificación taxonómica end-to-end, solo una porción.
- `tests/test_baseline_sr.cpp:123` verifica que `SrIndexLocator::locate_extremal ≡ sdsl::locate(fm,…)`, pero **nunca compara contra el LZ77-Index en pipeline LCA**.
- `tests/test_classify.cpp` corre el pipeline LCA completo, pero **solo con LZ77-Index**.

Sin un baseline que resuelva el mismo problema, no podemos afirmar que el sr-index "funciona" para clasificación taxonómica, ni cuantificar qué tan equivalente es la salida del LZ77-Index. La memoria necesita esa comparación para justificar las decisiones de diseño primarias-only.

Outcome buscado: una pipeline reproducible que tome un dataset con árbol + orden DFS + reads, ejecute clasificación con **ambos índices** y emita (a) tasa de coincidencia por categoría y (b) CSV de divergencias para análisis caso a caso.

## Decisiones de diseño

1. **Un solo tool C++ `compare_classifiers`** que corre ambos índices en una pasada. Comparte árbol, genome ranges, reads y MEMExtractor → garantiza que la comparación es justa por construcción.
2. **Árbol y rangos de genoma en TSV plano**, no JSON, para evitar agregar `nlohmann/json` como dependencia. Un script Python convierte `tree.json` → `tree.tsv` + `genomes.tsv`.
3. **Reads sintéticos nuevos** de ~250 bp con ground-truth `genome_idx`. Las `queries.txt` actuales (20-mers de la unidad base) NO sirven para clasificación: son demasiado cortas para `min_mem_len=31` y colapsan al LCA=root.
4. **Reutilizar `Classifier::classify_read<Index>` que ya es template** (`src/taxonomy/classifier.hpp:81`). Tanto `LZ77Index` como `SrIndexLocator` ya exponen `locate_extremal(pattern) → (min, max)` → la sustitución es directa.
5. **NO reconstruir sr-index sobre `data/repetitive/`** todavía: primero validar la cadena en el dataset sintético de 4 hojas que ya usa `test_classify.cpp`. Luego escalar.

## Archivos a crear

### `scripts/tree_json_to_tsv.py`
Lee `data/repetitive/tree.json` (formato de `scripts/gen_repetitive.py`) y emite:
- `data/repetitive/tree.tsv` — columnas `node_id`, `parent_id`, `name`, `genome_idx` (-1 para internos). Re-mapea ids si la raíz no es 0 (`PhyloTree::build` en `src/taxonomy/lca.hpp` requiere ids 0..n-1 consecutivos).
- `data/repetitive/genomes.tsv` — columnas `leaf_dfs_rank`, `start`, `end`, `node_id`. Calculado desde `leaf_dfs_order` y `params.genome_len` con separador `$` de 1 byte. Verificación: `end[-1] + 1 == len(reference.txt)`.

### `scripts/gen_classification_reads.py`
Genera `data/repetitive/reads.tsv`: `read_id`, `true_genome_idx`, `true_node_id`, `sequence`. N=2000 reads de 250 bp, posición aleatoria dentro de un genoma aleatorio, 1% SNPs (para que los MEMs ≥31 sigan apareando). Semilla fija (42).

### `tools/compare_classifiers.cpp`
Pipeline:
1. Carga `reference.txt` → `LZ77Index lz_idx; lz_idx.build(text)`.
2. Construye `sdsl::csa_wt<> fm` sobre el mismo texto + `MEMExtractor ext(fm)` (reutiliza patrón de `tests/test_classify.cpp:75-76`).
3. Carga sr-index pre-construido: `SrIndexLocator sr_loc; sr_loc.load(data_name, data_dir, sr)`.
4. Carga `tree.tsv` → `PhyloTree`; `genomes.tsv` → `vector<Classifier::GenomeRange>`; arma `Classifier`.
5. Para cada read en `reads.tsv`:
   - `lz = classifier.classify_read(seq, lz_idx, ext, 31)`
   - `sr = classifier.classify_read(seq, sr_loc, ext, 31)`
   - Categoría:
     - `EQUAL`: `lz == sr`
     - `LZ_DESCENDANT_OF_SR`: `tree.lca(lz, sr) == sr && lz != sr` (limitación primary-only **esperada**)
     - `SR_DESCENDANT_OF_LZ`: caso **no esperado** — potencial bug
     - `INCOMPARABLE`: ramas distintas — potencial bug
     - `EITHER_UNCLASSIFIED`: uno de los dos retornó -1
6. Compara también contra `true_node_id` (precisión absoluta vs ground-truth de generación).
7. Emite `results/compare/compare_<dataset>_sr<s>.csv` (una fila por read) y resumen agregado en stdout con porcentajes por categoría.

### `scripts/run_compare_classifiers.sh`
Wrapper que ejecuta el tool para s∈{16,32} sobre cada dataset disponible y deja los CSV en `results/compare/`.

## Archivos a modificar

### `CMakeLists.txt`
Agregar target `compare_classifiers` análogo al target existente `baseline_lca_sr`: linkea `src/index.cpp`, `src/lz77/*.cpp`, `src/wavelet/*.cpp`, `src/mem/extractor.cpp` y el objeto PIMPL `sr_index_locator_obj`. `src/taxonomy/{lca,classifier}.hpp` son header-only.

### `tests/test_baseline_sr.cpp` (o nuevo `tests/test_compare_indexes.cpp` si hay colisión de linkeo)
Fixture nuevo `ClassifyCompareTest` que replica el árbol de 4 hojas de `tests/test_classify.cpp:47-77`. Tests:
- `SrIndex_LCA_IsExact`: para cada patrón compartido, `classify_read(p, sr_loc, ext, 8) == bf_classify(p)`. El sr-index ve todas las ocurrencias → debe dar el LCA real exacto (esto es la **definición** de que el baseline está bien).
- `SrIndex_AgreesWithLZ77_OrIsAncestor`: para los mismos patrones, `tree.lca(lz_result, sr_result) == sr_result` (LZ77 está en subárbol del sr-index — propiedad débil ya probada en `tests/test_classify.cpp:209`, aquí extendida cruzando los dos índices).

## Orden de implementación y validación

1. **Conversor de árbol**: correr `tree_json_to_tsv.py`; validar a mano que `start[i+1] - end[i] == 1` para todo i y que `end[-1] + 1 == len(reference.txt)`.
2. **Generador de reads**: correr `gen_classification_reads.py`; `grep` algunos reads en `reference.txt` para confirmar posición/genoma.
3. **Tests unitarios** (fixture `ClassifyCompareTest`): compilar y correr `ctest`. **Deben pasar antes** de continuar — si el sr-index no da LCA exacto sobre 4 hojas, no tiene sentido escalar.
4. **Tool `compare_classifiers`**: compilar; correr primero en modo `--self-test` sobre el dataset de 4 hojas (replicado in-memory) para sanity-check; luego sobre `data/repetitive` (requiere construir antes el sr-index para `data/repetitive/reference.txt` — extender `scripts/run_sr_index.sh` para apuntarlo a ese dataset).
5. **Wrapper `run_compare_classifiers.sh`**: correr para s∈{16, 32}; los resultados de clasificación deben ser **idénticos** entre ambos s (s solo afecta tiempo/espacio del locate, no qué ocurrencias retorna).

## Verificación end-to-end (criterios de aceptación)

- En el fixture sintético: `SrIndex_LCA_IsExact` pasa 100% y `SrIndex_AgreesWithLZ77_OrIsAncestor` pasa 100%.
- En `data/repetitive`:
  - `EQUAL + LZ_DESCENDANT_OF_SR` ≥ 95% de reads clasificados.
  - `SR_DESCENDANT_OF_LZ == 0` y `INCOMPARABLE == 0`. Cualquier no-cero es bug y se diagnostica re-ejecutando esos `read_id` con `sdsl::locate` brute-force para identificar si el problema está en el sr-index, en el LZ77-Index, o en `MEMExtractor`.
  - Para reads donde `sr_node == true_node` (sr-index clavó el genoma exacto), el LZ77 debe estar en el mismo subárbol — verificación invariante.
- El CSV `results/compare/compare_repetitive_sr16.csv` permite filtrar por `category != EQUAL` para el análisis cualitativo caso-a-caso pedido (qué reads divergen, qué MEMs producen, en qué profundidad del árbol cae cada índice).
