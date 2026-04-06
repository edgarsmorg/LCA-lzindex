# Pipeline LCA con sr-index — Bosquejo de arquitectura y pesos

**Dataset:** Datos sintéticos (4 genomas, árbol de profundidad 2)  
**Fecha:** 2026-04-06  
**Propósito:** Baseline de comparación para el LZ77-index

---

## 1. Árbol filogenético y genomas de entrada

```
        root (0)
       /        \
     A (1)      B (2)
     / \        / \
  A1(3) A2(4) B1(5) B2(6)
```

| Genoma | Node ID | Secuencia                    | Largo |
|--------|---------|------------------------------|-------|
| A1     | 3       | `AAAAACGTACGTACGTACGT`        | 20    |
| A2     | 4       | `ACGTACGTACGTGGGGGGGG`        | 20    |
| B1     | 5       | `GGGGGGGGTTTTTTGCATGCA`       | 21    |
| B2     | 6       | `TGCATGCATGCACCCCCCCC`        | 20    |

**Texto de referencia** (concatenación en orden DFS, separador `#`):

```
AAAAACGTACGTACGTACGT#ACGTACGTACGTGGGGGGGG#GGGGGGGGTTTTTTGCATGCA#TGCATGCATGCACCCCCCCC#
^--- A1 [0,20) ---^ ^---- A2 [21,41) ----^ ^------- B1 [42,63) -------^ ^-- B2 [64,84) --^
```

- **n** = 85 caracteres (sin contar el `$` de la BWT)
- **Alfabeto** σ = {A, C, G, T, #}

---

## 2. Construcción de la BWT (ropebwt3)

```
ropebwt3 build -bo index.fmr reference.fa   # índice dinámico
ropebwt3 build -i index.fmr -do index.fmd   # convertir a estático
```

| Archivo        | Tamaño | Descripción                             |
|----------------|--------|-----------------------------------------|
| `reference.fa` | 101 B  | Genomas en formato FASTA                |
| `index.fmr`    | 421 B  | BWT dinámica (B+-tree RLE, formato .fmr)|
| `index.fmd`    | 272 B  | BWT estática comprimida (formato .fmd)  |

> **Nota:** ropebwt3 indexa forward + reverse complement por defecto.
> Para clasificación taxonómica se usa `-R` para solo forward.

**Estadísticas de la BWT:**
- Símbolos indexados: 170 (incluyendo RC)
- Distribución: ($=8, A=33, C=48, G=48, T=33, N=0)
- **r = 32 runs** de BWT (medido por sr-index)

---

## 3. Construcción del sr-index

```
bm_construct_ri -data=reference.txt -min_s=4 -max_s=16
```

### Estructuras internas (del log de memoria sdsl)

| Estructura                | Descripción                                    |
|---------------------------|------------------------------------------------|
| Mark2Sample Links         | Marcas de muestreo en runs de la BWT           |
| Predecessor               | Estructura de predecesor para navigate en BWT  |
| (unknown)                 | Estructura base del R-Index (FM-index RLBWT)   |

### Parámetros de muestreo medidos

| Variante       | s  | r' (runs submuestreados) | Tiempo construcción |
|----------------|----|--------------------------|---------------------|
| R-Index        | —  | — (r=32 runs totales)    | 1.26 ms             |
| SR-Index/4     | 4  | 21                       | 1.27 ms             |
| SR-Index/8     | 8  | 13                       | 1.32 ms             |
| SR-Index/16    | 16 | 8                        | 1.31 ms             |
| SR-Index-VA/4  | 4  | 21                       | 1.08 ms             |
| SR-Index-VA/8  | 8  | 13                       | 0.76 ms             |
| SR-Index-VA/16 | 16 | 8                        | 0.74 ms             |

> **SR-Index-VA** (Valid Area) es la variante recomendada: mejor relación espacio/tiempo.
> A mayor `s`, menos espacio pero locate más lento (factor s).

---

## 4. Pipeline de consulta (query → LCA)

```
┌─────────────────────────────────────────────────────────────────┐
│                     QUERY (read de ADN)                         │
└────────────────────────────┬────────────────────────────────────┘
                             │
                    ┌────────▼────────┐
                    │  Extracción     │
                    │  de MEMs        │  ← ropebwt3 mem -l31
                    │  (longitud ≥31) │
                    └────────┬────────┘
                             │  MEM = {query_start, ref_start, length}
                    ┌────────▼────────┐
                    │   sr-index      │
                    │   .locate(MEM)  │  ← O(r'/s) por ocurrencia
                    └────────┬────────┘
                             │  posiciones en texto ref: {p₁, p₂, ..., pₖ}
                    ┌────────▼────────┐
                    │  Mapeo pos →    │
                    │  genoma         │  ← bitvector rank sobre límites
                    │  (rank query)   │     de genoma [O(1) con sd_vector]
                    └────────┬────────┘
                             │  nodos hoja: {leaf₁, leaf₂, ..., leafₖ}
                    ┌────────▼────────┐
                    │   LCA(leaf_min, │
                    │       leaf_max) │  ← O(log n) naive / O(1) F-C&B
                    └────────┬────────┘
                             │
                    ┌────────▼────────┐
                    │  Clasificación  │
                    │  taxonómica     │
                    └─────────────────┘
```

### Costo por paso (datos sintéticos, n=85, r=32)

| Paso                 | Complejidad       | Estructuras usadas               |
|----------------------|-------------------|----------------------------------|
| Extracción MEMs      | O(n · m)          | BWT .fmd (ropebwt3)              |
| locate(MEM)          | O(m + occ·r'/s)   | SR-Index RLBWT + samples         |
| Mapeo pos→genoma     | O(1)              | `sd_vector` rank sobre B_genomes |
| LCA(leaf_min,leaf_max)| O(log n) naive   | PhyloTree (naive) / RMQ (prod.)  |

---

## 5. Estructuras de datos del pipeline completo

```
Disco (persistentes):
├── reference.txt          85 B    texto concatenado DFS
├── index.fmd             272 B    BWT comprimida (ropebwt3)
└── sr-index (en caché sdsl)
    ├── RLBWT runs          ~r·log(n) bits
    ├── Mark2Sample Links   ~r'·log(n) bits  (depende de s)
    └── Predecessor struct  ~r'·log(n) bits

RAM durante consulta:
├── SR-Index cargado        O(r + r'/s · log n) bits
├── BWT .fmd (mmap)        272 B
├── Bitvector B_genomes     n bits sparse → O(G·log(n/G)) bits Elias-Fano
└── PhyloTree               O(nodos · log n) bits
```

---

## 6. Validación de queries sintéticas (5/5 ✓)

| Query      | Patrón      | Posiciones     | Genomas tocados | LCA esperado | LCA obtenido |
|------------|-------------|----------------|-----------------|--------------|--------------|
| q1         | `ACGTACGT`  | 4,8,12,21,25   | A1, A2          | A            | A ✓          |
| q2         | `TGCATGCA`  | 55,64,68       | B1, B2          | B            | B ✓          |
| q3         | `GGGGGGGG`  | 33,42          | A2, B1          | root         | root ✓       |
| q4         | `AAAAACGT`  | 0              | A1              | A1           | A1 ✓         |
| q5         | `CCCCCCCC`  | 76             | B2              | B2           | B2 ✓         |

**Observación clave:** El LCA extremal (solo `pos_min` y `pos_max`) produce
resultados idénticos al LCA completo (todas las posiciones). Esto valida la
premisa del diseño primarias-only del LZ77-index.

---

## 7. Diferencia con el LZ77-index (nuestro índice)

| Aspecto              | sr-index (baseline)                | LZ77-index (nuestro)               |
|----------------------|------------------------------------|------------------------------------|
| Locate               | Todas las ocurrencias (occ total)  | Solo ocurrencias primarias (≤ z)   |
| Para LCA             | enumerate all → LCA                | pos_min + pos_max directo via RMQ  |
| Espacio              | O(r log n) bits                    | O(z log n) bits (z ≤ r típicamente)|
| Locate time          | O(m + occ · r'/s)                  | O(m log n + log z) amortizado      |
| Correctitud LCA      | Exacto siempre                     | Exacto para MEMs que cruzan límite |

> **Caso patológico LZ77-index:** MEMs completamente dentro de una frase
> no tienen ocurrencia primaria. Mitigación: usar MEMs de longitud ≥ longitud
> mínima de frase LZ77 (~log n).

---

## 8. Próximos pasos

- [ ] Implementar LZ77 parser sobre SA (Semana 3)
- [ ] Construir grilla 2D y Wavelet Tree (Semana 4-6)
- [ ] Conectar RMQ para localización extremal sin enumerate (Semana 7)
- [ ] Benchmark comparativo sr-index vs LZ77-index sobre corpus real (Semana 11-12)
