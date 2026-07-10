# Rediseño a dos índices (T y Tᴿ) — validación y espacio · 2026-07-10

Motivado por la retroalimentación de G. Navarro: la ocurrencia más a la derecha
NO es el máximo de las primarias del texto directo. Solución adoptada: dos
LZ-index (sobre T y sobre Tᴿ), **ambos con RMQ de mínimo**; el max-RMQ se eliminó.
`rightmost(P en T) = leftmost(Pᴿ en Tᴿ)`, remapeado con `s = (n-1) - m - j`.

## Correctitud (el punto serio, p19/p38) — RESUELTO

`tools/verify_lca_equiv` sobre `scaling10kb_10MB` (10 000 patrones):

| Diseño | EQUAL (de los que ocurren) | LZ_DESCENDANT | ANCESTOR/INCOMP |
|--------|---------------------------|---------------|-----------------|
| Antiguo (max de primarias) | ~94–99 % | ~1–6 % | 0 |
| **Nuevo (índice reverso)** | **100 %** (9000/9000) | **0 %** | 0 |

Los casos DESCENDANT del diseño antiguo eran exactamente las ocurrencias más a
la derecha secundarias que no se veían (la conjetura de Navarro en p38). Con el
índice reverso desaparecen. 136/136 tests unitarios verdes.

## Espacio (bpc) — el costo del rediseño

`tools/measure_index` sobre `scaling10kb_10MB` (n=10 MB, z/n=0.0011):

| Componente | bpc |
|------------|-----|
| Índice directo (grilla+min-RMQ+tries) | 0.22 |
| Índice reverso (idem) | 0.22 |
| **TOTAL dos índices** | **0.44** |
| Diseño antiguo (1 índice + max-RMQ), estimado | ~0.26 |
| **sr-index (s=16), mismo dataset** | **0.35** |

**Implicación:** el diseño ingenuo de dos índices independientes **duplica** el
bpc y queda POR ENCIMA del sr-index (0.44 > 0.35) a esta escala — se pierde la
ventaja de espacio que era resultado central de la memoria.

## Consecuencia para el diseño

Para preservar la ventaja de espacio hay que evitar pagar dos veces la maquinaria
de búsqueda. Dirección (idea Navarro + Edgar): reemplazar los tries (SST) por un
**bitmap disperso** de inicios de frase sobre el `SA_fwd`/`SA_rev` que la
construcción del **sr-index ya materializa** — no pagar un FM-index nuevo. Esto es
la Parte B del plan; con estos números pasa de "optimización a evaluar" a
**necesaria**. Requiere decidir si el FM-index/sr-index pasa a ser parte del índice
entregable (cambia la base de comparación).
