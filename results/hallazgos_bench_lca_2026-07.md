# Hallazgos de evaluación — suite 10 KiB, tiempo de LCA y verificación de correctitud

> Fecha: 2026-07-06. Documento para el agente de redacción de la memoria.
> Resume los resultados nuevos de esta iteración. Todos los CSV/PNG citados viven
> bajo `results/scaling/10kb/` (regenerables; `results/` está en `.gitignore`).

## TL;DR (para el escrito)

1. **Nueva suite de escalado con unidad base pequeña** (`bench/datasets_scaling_10kb.json`):
   genoma base de **10 KiB** repetido en copias mutadas (0.1 %) hasta cubrir de
   **10 MB a 1.8 GB** (8 tamaños). A diferencia de `datasets_scaling*.json`
   (genoma base 1 MiB), aquí la unidad repetida es 100× más chica → muchas más
   copias y más límites de frase LZ77 por byte (grilla 2D más densa).
2. **El índice LZ77 gana claramente en locate sobre el sr-index en colecciones
   repetitivas**: su tiempo de `locate_extremal` es casi plano al escalar,
   mientras el sr-index explota porque enumera todas las ocurrencias. A 1.8 GB,
   **≈110× más rápido** (277 µs vs 30 393 µs por consulta). El LZ77 también es
   más chico (0.15 vs 0.18 bpc).
3. **El paso de LCA (posición → genoma → árbol) es despreciable**: sub-microsegundo
   a ~3.7 µs por consulta, dominado por el locate, e independiente del índice.
4. **Correctitud del diseño primarias-only, verificada empíricamente contra el
   FM-index de sdsl (ground truth con TODAS las ocurrencias)**: el LCA de las
   ocurrencias primarias **nunca sale del linaje correcto** (0 casos ancestro /
   incomparable en cientos de miles de patrones), y **converge a exacto (EQUAL)
   al crecer la colección**: en 1.2/1.5/1.8 GB es 100 % en todos los largos de
   patrón. → *Descartar las ocurrencias secundarias no cuesta exactitud
   taxonómica en el régimen de interés (colecciones grandes y repetitivas).*

---

## 1. Qué se mide exactamente

- **`locate_extremal(P)`**: dado un patrón/MEM `P`, devuelve las **dos posiciones
  extremas en el texto** donde ocurre (`pos_min`, `pos_max`). NO enumera todas las
  ocurrencias ni las cuenta.
  - LZ77: obtiene los extremos **directo con RMQ** sobre la grilla de primarias.
    Costo función del **largo** de `P`, no de cuántas veces ocurre.
  - sr-index (baseline): `locate(P)` enumera **todas** las ocurrencias y toma
    `min`/`max`. Costo `O(occ)` → explota en textos muy repetitivos.
- **Paso de LCA**: `pos_min`/`pos_max` → genoma (búsqueda binaria sobre los rangos
  DFS, O(log g)) → `PhyloTree::lca`. Es lo que convierte "encontré el MEM" en
  "sé entre qué genomas está". Se cronometra **aparte** del locate.
- La clasificación end-to-end de una read = extraer MEMs → `locate_extremal` por
  MEM → LCA. Aquí se evalúa por-patrón (sin la extracción de MEMs).

## 2. Escalado LZ77 vs sr-index (suite 10 KiB)

Fuente: `results/scaling/10kb/scaling.csv`. Patrones: 10 000 de largo 32, 3 reps.
Figuras: `results/scaling/10kb/plots/scaling_{size,bpc,build,query,lca,total}.png`
(+ `*_last5.png`).

| tamaño | bpc LZ77 / sr | build s LZ77 / sr | **locate µs** LZ77 / sr | **LCA µs** LZ77 / sr |
|--------|---------------|-------------------|-------------------------|----------------------|
| 10 MB  | 0.260 / 0.346 | 0.8 / 0.6         | 34.7 / 51.1             | 0.10 / 0.07          |
| 50 MB  | 0.207 / 0.275 | 4.9 / 5.0         | 55.7 / 269.2            | 0.19 / 0.16          |
| 100 MB | 0.219 / 0.231 | 10.6 / 10.9       | 68.8 / 564.7            | 0.24 / 0.27          |
| 500 MB | 0.168 / 0.190 | 59.7 / 60.3       | 138.8 / 3 344.3         | 0.39 / 0.44          |
| 1 GB   | 0.157 / 0.193 | 137.7 / 136.5     | 188.0 / 8 910.8         | 0.52 / 1.28          |
| 1.2 GB | 0.153 / 0.190 | 177.7 / 180.4     | 196.8 / 12 508.1        | 0.59 / 2.30          |
| 1.5 GB | 0.149 / 0.187 | 232.1 / 257.2     | 210.4 / 19 657.0        | 0.63 / 3.25          |
| 1.8 GB | 0.150 / 0.184 | 297.1 / 367.9     | **277.4 / 30 389.5**    | 0.78 / 3.66          |

Lectura:
- **Locate**: LZ77 crece suave (34→277 µs); el sr-index se dispara (51→30 389 µs)
  porque enumera ~188 K ocurrencias por patrón (una por copia). Es la evidencia
  central a favor del diseño primarias-only con RMQ.
- **Espacio**: LZ77 consistentemente menor (≈0.15 vs 0.18 bpc a escala grande).
- **LCA**: sub-µs a ~3.7 µs → despreciable frente al locate; crece ~O(log g) más
  efectos de caché. Similar entre índices (trabajo posición→árbol, no del índice).
- **Construcción**: comparable entre ambos.

## 3. Verificación de correctitud (primarias-only vs ground truth)

Herramienta nueva: `tools/verify_lca_equiv.cpp`. Contrasta por patrón:
- **LZ77**: `locate_extremal` (extremos de primarias) → `LCA(genoma(min), genoma(max))`.
- **Ground truth**: **FM-index de sdsl** (`csa_wt`) enumera TODAS las ocurrencias,
  mapea cada una a su genoma y calcula el LCA del conjunto completo.

Categorías: `EQUAL` (idéntico), `LZ_DESCENDANT` (LZ77 más específico, mismo linaje),
`LZ_ANCESTOR`/`INCOMPARABLE` (**bugs** — LZ77 fuera del linaje correcto).

Fundamento: como la referencia sigue orden DFS, el LCA de todas las ocurrencias
== LCA de sus extremos. El LZ77 solo ve primarias (un subconjunto), así que su LCA
es EQUAL o un **descendiente** del real — nunca un ancestro ni incomparable.

### 3a. Convergencia con el largo del patrón (10 MB, figura de barras)

`results/scaling/10kb/lca_verify/lca_verify_10MB.png` — 5 000 patrones/largo:

| largo | 16 | 24 | 32 | 48 | 64 | 96 | 128 | 192 | 256 |
|-------|----|----|----|----|----|----|-----|-----|-----|
| EQUAL % | 37 | 50 | 61 | 75 | 84 | 93 | 96 | 98 | 98 |
| bugs   | 0  | 0  | 0  | 0  | 0  | 0  | 0   | 0   | 0   |

A mayor largo de MEM, más cruza límites de frase → las primarias capturan los
extremos → más EQUAL. **10 MB es demasiado chico para ser representativo**
(dismisado por decisión del autor); se usa solo para ilustrar la dependencia con
el largo.

### 3b. Escalado 100 MB → 1.8 GB (figura principal de correctitud)

`results/scaling/10kb/lca_verify/lca_verify_scale_100mb.png` — % EQUAL vs largo,
una línea por colección. Patrones capados por costo del ground truth (100/500 MB:
1 500–3 000; 1–1.8 GB: 300–500):

| colección (copias)   | @16  | @32  | @64  | @128 | @256 | bugs |
|----------------------|------|------|------|------|------|------|
| 100 MB (10 240)      | 99.3 | 99.3 | 98.1 | 96.5 | 94.0 | 0    |
| 500 MB (51 200)      | 100  | 99.9 | 99.7 | 99.7 | 99.2 | 0    |
| 1 GB (104 858)       | 100  | 100  | 100  | 100  | 99.8 | 0    |
| 1.2 GB (125 829)     | 100  | 100  | 100  | 100  | 100  | 0    |
| 1.5 GB (157 286)     | 100  | 100  | 100  | 100  | 100  | 0    |
| 1.8 GB (188 744)     | 100  | 100  | 100  | 100  | 100  | 0    |

**Resultado clave**: a más copias, más EQUAL; desde ~500 MB es ≥99 % y en
1.2/1.5/1.8 GB es **exacto (100 %) en todos los largos**. Intuición: más copias →
más límites de frase → las primarias casi siempre alcanzan los genomas extremos.
Cuando difieren, el LZ77 da un **descendiente** (clasificación más específica sobre
el linaje correcto), nunca un error de rama.

Verificación de invariante en toy con ground-truth conocido:
`tests/test_compare_indexes.cpp` (13 casos, LZ siempre en subárbol del sr-index,
sr == FM-index).

## 4. Cambios de código relevantes (para no citar cosas viejas)

**Añadido**:
- `tools/verify_lca_equiv.cpp` — verificación LCA primarias vs FM-index (con
  salida CSV y barrido multi-largo `--pat=largo:archivo`).
- `scripts/plot_lca_verify.py` — barras apiladas (un tamaño) y líneas multi-tamaño
  (`--series`). Paleta validada con el skill dataviz.
- `bench/datasets_scaling_10kb.json` — suite de escalado 10 KiB → 1.8 GB.
- `tools/bench_locate_compare.cpp` — ahora reporta `lca_us` y `total_us` además de
  `locate_ext_us`.
- `scripts/{plot_scaling,run_scaling_bench,run_classify_bench}.py/.sh`,
  `tools/{baseline_lca_lz,bench_classify_compare,classify_io}`, `vendor/lz77index`
  (tries Patricia DFUDS + DAC).

**Eliminado (no citar)**: `src/main.cpp`, `src/wavelet/wt_rmq_min.*` (WtMinRmq
legacy), `tools/bench_special.cpp`, `tools/sr_index_template.cpp`,
`tests/test_wt_rmq_min.cpp`. La RCSA (CSA fwd/rev del sr-index) ya se había
eliminado en el audit 2026-06-29 (el matching es por tries Patricia DFUDS).

## 5. Notas metodológicas

- **Ground truth = FM-index de sdsl** (`sdsl::csa_wt<>` + `sdsl::locate`), el
  estándar del proyecto para correctitud.
- **Costo del ground truth a escala**: como el texto es hiper-repetitivo, cada
  patrón toca ~todas las copias; `locate` enumera cientos de miles de posiciones.
  Por eso `verify_lca_equiv` construye LZ+FM **una sola vez** y barre varios
  largos, y se capa el nº de patrones en los tamaños grandes.
- **Gráficos**: usar el intérprete `../.venv/bin/python` (tiene matplotlib; el
  `python3` del sistema no).
- **Recursos**: la construcción a 1.8 GB usa un SA a 64 bits (~14 GB solo forward);
  se evita `scaling_2GB` por un segfault conocido en construcción.
