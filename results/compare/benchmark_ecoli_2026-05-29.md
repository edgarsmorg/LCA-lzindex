# Benchmark comparativo: LZ77-index vs sr-index
**Dataset:** Escherichia_Coli (n = 112,689,515 bytes = 107.47 MB)  
**Fecha:** 2026-05-29  
**Máquina:** Linux x86-64, GCC, -O3 -march=native  
**Patrones:** extraídos aleatoriamente del propio texto (seed=42)

---

## Espacio en disco (serialize en memoria)

| Índice          | Bytes        | MB     | bpc  |
|-----------------|-------------|--------|------|
| sr-index sr=4   | 78,484,583  | 74.85  | 5.57 |
| sr-index sr=8   | 54,688,047  | 52.15  | 3.88 |
| sr-index sr=16  | 41,775,207  | 39.84  | 2.97 |
| sr-index sr=32  | 35,187,791  | 33.56  | 2.50 |
| sr-index sr=64  | 30,896,415  | 29.47  | 2.19 |
| **LZ77-index**  | ~107,000,000 | ~107  | 7.98 |

LZ77-index desglose (7.98 bpc total):
- wt_int (grilla / count): 0.49 bpc
- sd_vector fwd + rev: 0.28 bpc
- text_pos (int_vector): 0.42 bpc
- phrase_total_len: 0.31 bpc
- wm_min_rmq (bv+rank+rmq): 1.37 bpc
- wm_max_rmq (bv+rank+rmq): 1.37 bpc
- CSA forward: 1.88 bpc
- CSA reverse: 1.88 bpc

---

## Tiempo locate_extremal (µs/query) — 1000 patrones mezclados m=20–50

| Índice         | µs/query | Peak RSS |
|----------------|----------|----------|
| sr-index sr=4  | 27.68    | 246 MB   |
| sr-index sr=8  | 25.59    | 143 MB   |
| sr-index sr=16 | 27.09    | 122 MB   |
| sr-index sr=32 | 29.30    | 116 MB   |
| sr-index sr=64 | 34.48    |  79 MB   |
| **LZ77-index** | **649**  | 2307 MB  |

## Tiempo locate_extremal (µs/query) — patrones largos

| m      | sr=16  | sr=64  | LZ77-index |
|--------|--------|--------|------------|
| 100    | 60.0   | 60.7   | —          |
| 200    | 108.6  | 98.7   | —          |
| 500    | 214.3  | 230.2  | —          |
| 1000   | 426.1  | 427.1  | —          |
| 100–1000 mix | 194.6 | 196.2 | **~231,000** |

---

## Análisis

### Espacio
El sr-index es ~3.6× más comprimido que el LZ77-index (2.19 vs 7.98 bpc).
El sr-index es O(r log n) con r = runs de BWT; el LZ77-index es O(z log n) con
z = frases LZ77. Para DNA altamente repetitivo, r << z. En E. coli, z = 1,752,702
frases. Con sr=64, el sr-index ocupa sólo 29 MB.

### Tiempo
El sr-index escala O(m · r/sr) en backward search + O(occ) en enumerate.
Para locate_extremal se enumeran todas las ocurrencias (O(occ · sr) steps de Phi).
El tiempo crece ~linealmente con m: ~0.42 µs/carácter.

El LZ77-index hace m−1 splits, cada uno con dos búsquedas RCSA (una fwd de
longitud m−i y una rev de longitud i). Costo total: O(m²·t_psi + m·log²z).
Para m=50: ~649 µs. Para m≈450 (mix): ~231,000 µs. Ratio 356× vs predicción
cuadrática 81×: el costo del RCSA también crece con m.

### Conclusión
El LZ77-index no es competitivo contra el sr-index en E. coli en ninguna métrica.
El trade-off teórico (no enumerar secundarias, RMQ directo) no compensa el costo
O(m²) de las m−1 búsquedas RCSA. El sr-index con sr=16–64 domina en espacio
(2.19–2.97 bpc) y en tiempo (27–427 µs/query para m=50–1000).
