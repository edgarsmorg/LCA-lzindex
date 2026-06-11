#pragma once

/**
 * test_helpers.hpp — utilidades compartidas entre los tests de lz77tax.
 *
 * Proporciona:
 *   1. make_toy_tree()                  — árbol filogenético de juguete de 7 nodos.
 *   2. make_toy_genome_ranges()         — rangos de genoma para el árbol toy
 *                                         con texto de 85 caracteres (semántica
 *                                         usada en test_classify, test_mem y
 *                                         test_baseline_sr).
 *   3. make_toy_genome_ranges_exclusive() — mismos genomas con rangos con un
 *                                           offset diferente (semántica usada en
 *                                           test_compare_indexes).
 *   4. BfResult / brute_force_range_2d() — argmin brute-force para verificación
 *                                          de WtMinRmq y WmMinRmq.
 *
 * Topología del árbol toy:
 *
 *        root(0)
 *       /       \
 *     A(1)      B(2)
 *    /    \    /    \
 *  A1(3) A2(4) B1(5) B2(6)
 *
 * Texto DFS (separador '#'):
 *   A1: "AAAAACGTACGTACGTACGT"  [0..19]
 *   A2: "ACGTACGTACGTGGGGGGGG"  [21..40]
 *   B1: "GGGGGGGGTTTTTTGCATGCA" [42..62]
 *   B2: "TGCATGCATGCACCCCCCCC"  [64..83]
 */

#include "../src/taxonomy/lca.hpp"
#include "../src/taxonomy/classifier.hpp"

#include <climits>
#include <vector>

namespace lz77tax::test {

// ─────────────────────────────────────────────────────────────────────────────
// 1. Árbol filogenético toy
// ─────────────────────────────────────────────────────────────────────────────

/// Construye y retorna el árbol toy de 7 nodos usado en la mayoría de tests.
inline PhyloTree make_toy_tree() {
    PhyloTree t;
    t.build(
        {0, 1, 2, 3, 4, 5, 6},
        {"root", "A", "B", "A1", "A2", "B1", "B2"},
        {PhyloTree::NO_PARENT, 0, 0, 1, 1, 2, 2}
    );
    return t;
}

// ─────────────────────────────────────────────────────────────────────────────
// 2. Rangos de genoma — variante usada en test_classify / test_mem /
//    test_baseline_sr (end ligeramente diferente al de test_compare_indexes).
// ─────────────────────────────────────────────────────────────────────────────

/**
 * Rangos correspondientes al árbol toy con el texto de 85 chars.
 *
 * Semántica: GenomeRange::end es el primer índice FUERA del genoma.
 *
 *   A1: [0, 20)   → posiciones 0..19
 *   A2: [21, 41)  → posiciones 21..40
 *   B1: [42, 63)  → posiciones 42..62
 *   B2: [64, 84)  → posiciones 64..83
 *
 * Usados en: test_classify.cpp, test_mem.cpp, test_baseline_sr.cpp
 */
inline std::vector<Classifier::GenomeRange> make_toy_genome_ranges() {
    return {
        {  0, 20, 3},   // A1: [0..19]
        { 21, 41, 4},   // A2: [21..40]
        { 42, 63, 5},   // B1: [42..62]
        { 64, 84, 6},   // B2: [64..83]
    };
}

/**
 * Rangos con límites extendidos en uno (end = start_siguiente_genoma).
 *
 *   A1: [0, 21)   → inclusive del separador '#' en pos 20
 *   A2: [21, 42)  → inclusive del separador '#' en pos 41
 *   B1: [42, 64)  → inclusive del separador '#' en pos 63
 *   B2: [64, 85)  → inclusive del separador '#' en pos 84
 *
 * Usados en: test_compare_indexes.cpp
 */
inline std::vector<Classifier::GenomeRange> make_toy_genome_ranges_exclusive() {
    return {
        {  0, 21, 3},   // A1: [0..20] exclusive
        { 21, 42, 4},   // A2: [21..41] exclusive
        { 42, 64, 5},   // B1: [42..63] exclusive
        { 64, 85, 6},   // B2: [64..84] exclusive
    };
}

// ─────────────────────────────────────────────────────────────────────────────
// 3. Brute-force argmin 2D
// ─────────────────────────────────────────────────────────────────────────────

/**
 * Resultado del argmin/argmax brute-force sobre un rectángulo 2D.
 *
 * - count     : número de puntos (x, ys[x]) con x in [x_lo,x_hi] y ys[x] in [y_lo,y_hi].
 * - argmin    : índice x en [0,n) del punto con vals[x] mínimo (SIZE_MAX si count==0).
 * - argmin_val: valor vals[argmin] (SIZE_MAX si count==0).
 */
struct BfResult {
    size_t count;
    size_t argmin;      ///< índice global en [x_lo, x_hi]
    size_t argmin_val;  ///< vals[argmin]
};

/**
 * Escanea todos los puntos (x, ys[x]) con x in [x_lo, x_hi]
 * que satisfacen ys[x] in [y_lo, y_hi] y retorna el argmin de vals[x].
 *
 * Equivalente a las funciones locales brute_force() / brute_force_argmin()
 * que existían en test_wm_rmq_min.cpp y test_wt_rmq_min.cpp.
 *
 * @param ys    coordenada Y de cada punto, indexado por X.
 * @param vals  peso/valor de cada punto, indexado por X.
 * @param x_lo  límite inferior del rango X (inclusivo).
 * @param x_hi  límite superior del rango X (inclusivo).
 * @param y_lo  límite inferior del rango Y (inclusivo).
 * @param y_hi  límite superior del rango Y (inclusivo).
 */
inline BfResult brute_force_range_2d(
        const std::vector<size_t>& ys,
        const std::vector<size_t>& vals,
        size_t x_lo, size_t x_hi,
        size_t y_lo, size_t y_hi)
{
    size_t count = 0;
    size_t best_v = SIZE_MAX;
    size_t best_i = SIZE_MAX;

    for (size_t x = x_lo; x <= x_hi; ++x) {
        if (ys[x] < y_lo || ys[x] > y_hi) continue;
        ++count;
        if (vals[x] < best_v) {
            best_v = vals[x];
            best_i = x;
        }
    }
    return {count, best_i, best_v};
}

} // namespace lz77tax::test
