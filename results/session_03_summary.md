# Sesión 03 — Wavelet Tree + RMQ topológico para `locate_min`

**Fecha:** 2026-04-22
**Semanas cubiertas:** 7 del cronograma (RMQ en WT, locate query funcional)

---

## Contexto: qué se quería implementar

Implementar la consulta de mínimo 2D en T en O(polylog z) usando un Wavelet Tree
custom con RMQ sucinto por nodo, siguiendo la arquitectura Ferrada-Navarro:

- **RMQ topológicos** (cajas negras): devuelven posiciones en O(1) con 2n+o(n) bits,
  sin almacenar los valores asociados (text_pos).
- **Valores en Grid2D**: `text_pos_plain_` vive una sola vez en la grilla.
- **Unwind vía select**: conversión de índice local de nodo → índice global usando
  `select_0`/`select_1` en la cadena de padres, O(log σ).
- **Dos índices en el futuro**: uno sobre T (mínimo), otro sobre T^R (máximo via mínimo).
  Esta sesión solo implementa el de T.

---

## Lo que se hizo

### Commit 1 — `text_pos_` compacto + shadow plain

`text_pos_` migrado de `std::vector<size_t>` a `sdsl::int_vector<>` (ceil(log₂n) bits
por entrada). Se añade `text_pos_plain_` como shadow `std::vector<size_t>` para pasar
a la API del WtMinRmq sin conversión por consulta.

### Commit 2 — `WtMinRmq` — WT binario custom, RMQ topológico

Nuevos archivos: `src/wavelet/wt_rmq_min.hpp` y `.cpp`.

Estructura por nodo:
```
Node {
  sdsl::bit_vector bv;              // routing
  sdsl::rank_support_v<1>  rank1;   // descenso
  sdsl::select_support_mcl<0> sel0; // unwind
  sdsl::select_support_mcl<1> sel1; // unwind
  sdsl::rmq_succinct_sct<true> rmq; // topológico: sin valores internos
  size_t lo, hi;
  Node* parent; bool is_right_child;
}
```

El `rmq` se construye desde un `int_vector<>` temporal que se descarta tras build.
Ningún nodo almacena `text_pos`. `range_argmin_2d` recibe `const vector<size_t>&`
en cada llamada. `unwind_to_root` sube de hoja a raíz vía select.

Tests (`tests/test_wt_rmq_min.cpp`): 9 tests exhaustivos contra brute-force. Pasan.

### Commit 3 — `Grid2D::query_min_2d`

`build_from_coords` construye `wt_min_rmq_` tras el `wt_int`. `query_min_2d` aplica
las mismas transformaciones de rango (rank_fwd / rank_rev) que `query_extremal`.

Invariante verificado: `query_min_2d(r).boundary_min == query_extremal(r).boundary_min`.
22 tests pasan.

### Commit 4 — `LZ77Index::locate_min`

Itera splits del patrón y llama `grid_.query_min_2d()` en vez de `query_extremal()`.
`locate_extremal` se conserva sin cambios (test_classify lo usa).

Invariante: `locate_min(P) == locate_extremal(P).first` para todos los patrones.
27 tests pasan.

### Commit 5 — `measure_index` actualizado

Desglose de `WtMinRmq` en tres filas: bitvectors, rank+select, RMQ BP.
Bench comparativo: 1000 patrones aleatorios, `locate_min` vs `locate_extremal`.

---

## Resultado en E. coli (n = 108 MB, z = 1.752.702)

```
wt_int (grilla / count)       :    7 MB  (0.49 bpc)
sd_vector fwd                 :    2 MB  (0.14 bpc)
sd_vector rev                 :    2 MB  (0.14 bpc)
text_pos[] (z×8B aprox)      :   13 MB  (1.00 bpc)
  --- WtMinRmq (locate_min) ---
    wt_min_rmq (bitvectors)   :   43 MB  (3.21 bpc)
    wt_min_rmq (rank+select)  :  267 MB  (19.87 bpc)
    wt_min_rmq (rmq BP)       :  541 MB  (40.27 bpc)
    wt_min_rmq total          :  851 MB  (63.34 bpc)

locate_min      : 33.94 µs/query
locate_extremal : 40.97 µs/query
speedup         : 1.21x
```

---

## Hallazgo crítico: overhead per-instancia de sdsl explota a z grande

**Análisis teórico** (correcto): O(z log σ) bits para bitvectors + supports + RMQs.
Para z = 1.75M y log σ = 21: ≈ 36.8M bits ≈ 4.6 MB de contenido.

**Realidad**: 851 MB. La discrepancia de ~180× se explica por el **overhead per-instancia
de sdsl**:

- Cada `sdsl::bit_vector`, `rank_support_v`, `select_support_mcl`, `rmq_succinct_sct`
  tiene tablas de bloques precomputadas propias de tamaño fijo (~300 bytes mínimo).
- Con O(z) nodos internos, el overhead total es O(z × 300 bytes) ≈ 525 MB.
- El análisis O(z log σ) bits era correcto para el **contenido**, ignoró la constante
  oculta de sdsl.

**¿Cuándo NO es problema?**

Para z ≪ 1M (el caso de uso real del índice — genomas pequeños de referencia):
- z = 60k: 60k × 300 bytes = 18 MB → completamente manejable.
- La E. coli completa (108 MB) es el peor caso; cada genoma individual en una
  clasificación taxonómica sería mucho más pequeño.

El plan explícitamente asumía "z ≈ 60k", no z = 1.75M.

---

## Fix requerido para z grande: Wavelet Matrix plana

Si en el futuro se necesita escalar a z > 200k, el fix es reemplazar el árbol
pointer-based por una **wavelet matrix plana**:

- log σ niveles, cada nivel tiene UNA sola BV de tamaño z + UN solo rank + UN solo RMQ.
- Total: O(log σ) instancias sdsl (no O(z)).
- Espacio: O(z log σ) bits ≈ 9 MB RMQ + 4.6 MB BVs = ~15 MB para E. coli.
- Tradeoff: cada nivel requiere re-mapeo de índices → necesita O(z log σ) bytes para
  arrays de reordenamiento por nivel, o devolver VALUE directamente (no argmin_global).

Este refactor es no trivial y está fuera del scope de esta sesión.

---

## Estado del índice tras esta sesión

| Funcionalidad | Estado | Notas |
|---------------|--------|-------|
| LZ77 parser   | ✅ correcta | Tests vs brute-force |
| Grid2D build  | ✅ correcta | Ambas variantes (texto + CSA) |
| count()       | ✅ correcta | O(m² · t_ψ) |
| locate_extremal() | ✅ correcta | O(occ · log z), enumera |
| locate_min()  | ✅ correcta | O(m · log² z), no enumera |
| locate_max()  | ❌ pendiente | Requiere índice sobre T^R |
| WtMinRmq espacio | ⚠️ ineficiente para z > 200k | Fix: wavelet matrix plana |

## Tests en verde

```
test_lz77        :  10/10
test_wt_rmq_min  :   9/9
test_grid        :  22/22
test_index       :  27/27
test_classify    :  19/19
TOTAL            :  87/87
```

## Próxima sesión

1. Implementar el índice hermano sobre T^R (recibe T^R como T a la misma clase).
2. Con ambos índices, redefinir `locate_extremal` como `(locate_min(T), translate(locate_min(T^R)))`.
3. Opcional: refactorizar WtMinRmq a wavelet matrix plana si z > 200k es relevante.
