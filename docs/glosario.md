# Glosario de nomenclatura

Este documento fija una traducción práctica entre los nombres usados en el proyecto, la literatura de índices comprimidos y la implementación actual. La idea no es reemplazar la teoría, sino dar un mapa mental para leer el código sin perderse.

## 1. Objetos base del proyecto

### Texto de referencia
Texto concatenado sobre el cual se construye el índice. En este proyecto corresponde a la concatenación de genomas de referencia en orden DFS del árbol filogenético.

### Query
Lectura o patrón que se quiere buscar sobre la referencia.

### Patrón
Sinónimo de query cuando se habla desde el punto de vista del índice. En el código suele aparecer como `pattern`.

### Genoma
Secuencia individual dentro del texto concatenado. Cada genoma ocupa un rango `[start, end)` en la referencia.

### Nodo filogenético
Nodo del árbol taxonómico. Un genoma hoja se representa mediante un `node_id`.

### LCA
`Lowest Common Ancestor`. Dado un conjunto de genomas o nodos del árbol, es el ancestro común más bajo. Es la salida final de la clasificación taxonómica.

## 2. Estructuras clásicas de indexación

### SA
`Suffix Array`. Arreglo de sufijos.

Es un arreglo que lista las posiciones iniciales de todos los sufijos del texto en orden lexicográfico. Si `SA[i] = p`, entonces el sufijo `T[p..n-1]` es el i-ésimo sufijo más pequeño.

Sirve para buscar patrones como rangos contiguos de sufijos.

### ISA
`Inverse Suffix Array`. Inversa del suffix array.

Responde la pregunta inversa: dado un sufijo que empieza en la posición `p`, ¿en qué lugar del `SA` aparece? Formalmente, `ISA[p] = i` si `SA[i] = p`.

En este proyecto aparece en la definición geométrica de la grilla:
- `X_k` usa `ISA` del texto forward.
- `Y_k` usa `ISA` del texto reverso.

### BWT
`Burrows-Wheeler Transform`.

Transformación del texto estrechamente ligada al `SA`. Reordena caracteres según el contexto de sus sufijos. No es solo compresión: también permite construir índices de búsqueda eficientes.

### RLBWT
`Run-Length BWT`.

Representación comprimida de la BWT por runs. Es importante cuando el texto es repetitivo. El baseline `sr-index` y herramientas como `ropebwt3` viven cerca de esta idea.

### FM-index
Índice comprimido basado en BWT que soporta búsqueda por conteo y localización. En la práctica, permite hacer `count` y `locate` sin almacenar el texto crudo ni el `SA` completo de forma explícita.

### CSA
`Compressed Suffix Array`.

Familia de estructuras que representan el `SA` de forma comprimida, normalmente apoyadas en la BWT. En `sdsl-lite`, `csa_wt<>` es una implementación de CSA montada sobre wavelet tree.

En el proyecto:
- `csa_fwd_` indexa `T`.
- `csa_rev_` indexa `T^R`.

### Backward search
Técnica estándar sobre FM-index/CSA para buscar un patrón extendiéndolo de derecha a izquierda. El resultado es un intervalo `[sp, ep]` del `SA` que contiene todos los sufijos que tienen el patrón como prefijo.

## 3. Conceptos LZ77

### Parsing LZ77
Descomposición greedy del texto en frases. El proyecto usa esa descomposición como esqueleto del índice.

### Frase
Unidad del parsing LZ77. En el código se representa con `Phrase`.

En [phrase.hpp](/home/edgarmorales/uni/memoria/LCA-lzindex/src/lz77/phrase.hpp:1), una frase tiene:
- `start_pos`: dónde empieza en el texto.
- `length`: largo de la copia.
- `next_char`: carácter explícito final.

En esta versión del proyecto no se guarda `source`, porque no se implementa propagación de ocurrencias secundarias.

### Límite de frase
Separación entre una frase y la siguiente. En la documentación y en el código suele aparecer como `boundary`.

Si la frase `k` termina y luego empieza la frase `k+1`, el boundary relevante para la grilla está en `start_{k+1}`.

### Ocurrencia primaria
Ocurrencia de un patrón que cruza al menos un límite de frase. El proyecto mantiene solo estas ocurrencias porque para clasificación taxonómica bastan los extremos.

### Ocurrencia secundaria
Ocurrencia que no se obtiene directamente por cruce de boundary y que en el LZ-index clásico se propaga desde las primarias. Este proyecto no la construye.

### `z`
Cantidad de frases del parsing LZ77. Es una medida estándar de compresibilidad del texto repetitivo.

## 4. Grilla 2D del índice

### Grilla 2D
Estructura geométrica donde cada boundary entre frases se convierte en un punto. Es el corazón del índice de ocurrencias primarias.

### Punto
Representación geométrica de un boundary. Si hay `z` frases, hay `z-1` puntos.

### Coordenada `X_k`
Rango en el `SA` forward del sufijo que empieza en `start_{k+1}`.

Formalmente:
`X_k = ISA[start_{k+1}]`

### Coordenada `Y_k`
Rango en el `SA` reverso del prefijo que termina justo antes del boundary.

Formalmente:
`Y_k = ISA_rev[n - 2 - end_k]`

### `wt_`
Wavelet tree que almacena la grilla en forma compacta. Conceptualmente contiene un arreglo `R` tal que, para cada punto ordenado por `X`, guarda su coordenada `Y`.

### `bv_fwd_` y `bv_rev_`
Bitvectors sparse que marcan qué posiciones globales del `SA` y del `rev-SA` participan realmente en la grilla. Sirven para convertir rangos globales del índice a rangos relativos sobre el wavelet tree.

### `wt_idx`
Índice interno dentro del arreglo almacenado por el wavelet tree. No es una posición del texto, ni un índice del `SA`, ni un `node_id`.

### `text_pos_`
Vector auxiliar en la grilla. Para cada `wt_idx`, guarda la posición del boundary correspondiente en el texto. Semánticamente sería más claro pensar en esto como "boundary position".

## 5. Rangos e intervalos

### Intervalo `[sp, ep]`
Rango cerrado de posiciones en un `SA` o `CSA` devuelto por `backward_search`. Contiene todos los sufijos compatibles con un patrón o con una mitad del patrón.

### Rango BWT global
Intervalo `[sp, ep]` expresado en coordenadas del índice completo.

### Rango relativo en la grilla
Intervalo convertido con `rank` sobre `bv_fwd_` o `bv_rev_`, para consultar el `wt_` solo sobre los puntos válidos de la grilla.

### Rectángulo de búsqueda
Región `[x1, x2] x [y1, y2]` en la grilla. Corresponde a un split del patrón:
- lado derecho buscado en `csa_fwd`
- lado izquierdo reversado y buscado en `csa_rev`

Los puntos dentro del rectángulo representan ocurrencias primarias compatibles con ese split.

## 6. MEMs y clasificación

### MEM
`Maximal Exact Match`.

Match exacto entre query y referencia que no puede extenderse ni a izquierda ni a derecha sin producir un mismatch. En el código se representa con `MEM` en [mem.hpp](/home/edgarmorales/uni/memoria/LCA-lzindex/src/mem/mem.hpp:1).

### `query_start`
Posición donde parte el MEM dentro de la query.

### `ref_start`
Posición donde parte el MEM dentro de la referencia.

### `GenomeRange`
Rango `[start, end)` que ocupa un genoma dentro del texto concatenado, junto con el `node_id` de la hoja correspondiente.

### `genome_of(pos)`
Operación que mapea una posición del texto concatenado al genoma que la contiene.

### Clasificación extremal
Estrategia que usa solo la posición mínima y máxima de las ocurrencias primarias para computar el LCA. Funciona porque la concatenación de genomas sigue orden DFS.

## 7. Traducción rápida de nombres del código

### `Phrase`
Frase LZ77.

### `LZ77Parsing`
Vector completo de frases del texto.

### `LZ77Parser`
Constructor del parsing LZ77.

### `Grid2D`
Índice geométrico de puntos primarios.

### `count(pattern)`
Cuenta ocurrencias primarias del patrón, no necesariamente todas las ocurrencias textuales.

### `locate_extremal(pattern)`
Devuelve posiciones extremales de ocurrencias primarias del patrón en el texto.

### `boundary_min` / `boundary_max`
Extremos de boundaries encontrados en la grilla antes de traducirlos a posición real de ocurrencia.

### `pos_min` / `pos_max`
Extremos de ocurrencias reales en el texto, después de ajustar por el split del patrón.

### `csa_fwd`
CSA del texto original.

### `csa_rev`
CSA del texto reversado.

### `node_id`
Identificador de nodo en el árbol filogenético.

## 8. Regla práctica para no perderse

Al leer el código, conviene preguntarse siempre en cuál de estas cuatro capas estás:

1. Texto y genomas: posiciones reales en la referencia concatenada.
2. Parsing LZ77: frases y boundaries.
3. Índice geométrico: puntos, rangos `SA/CSA`, wavelet tree.
4. Taxonomía: genomas, nodos y `LCA`.

La mayoría de la confusión aparece cuando un nombre de una capa se interpreta como si perteneciera a otra.

## 9. Renames utiles para codigo

Esta sección no propone cambios obligatorios. Son renombres sugeridos para hacer más explícita la capa conceptual a la que pertenece cada dato.

### En parsing LZ77
- `start_pos` -> `phrase_start`
- `length` -> `phrase_len`
- `phrases` -> `lz_phrases`

### En boundaries y posiciones de texto
- `text_pos_` -> `boundary_pos_`
- `text_pos(wt_idx)` -> `boundary_pos(wt_idx)`
- `boundary_min` -> `min_boundary_pos`
- `boundary_max` -> `max_boundary_pos`
- `pos_min` -> `min_occ_start`
- `pos_max` -> `max_occ_start`

### En rangos del indice
- `sp_right` -> `fwd_sp`
- `ep_right` -> `fwd_ep`
- `sp_left` -> `rev_sp`
- `ep_left` -> `rev_ep`
- `query(...)` -> `query_rectangle(...)`
- `query_extremal(...)` -> `query_rectangle_extremal(...)`

### En la grilla
- `coords` -> `grid_points`
- `boundaries` -> `boundary_positions`
- `wt_idx` -> `grid_idx`

### En taxonomia
- `result` -> `lca_node`
- `genome_node` -> `genome_leaf`
- `genomes_` -> `genome_ranges_`

### Regla general de prefijos
Si en el futuro quieres ordenar más la nomenclatura, una convención razonable sería usar prefijos por capa:
- `phrase_*` para parsing LZ77
- `boundary_*` para límites entre frases
- `occ_*` para ocurrencias reales en el texto
- `fwd_*` y `rev_*` para rangos del índice forward/reverso
- `genome_*` para rangos de genomas
- `node_*` para nodos del árbol filogenético
