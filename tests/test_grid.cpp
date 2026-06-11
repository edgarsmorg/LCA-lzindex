#include <gtest/gtest.h>

#include "lz77/parser.hpp"
#include "lz77/grid.hpp"

#include <algorithm>
#include <set>
#include <string>
#include <vector>

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/csa_wt.hpp>

using namespace lz77tax;

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Añade terminador nulo al texto para que el SA sea correcto.
/// El '\0' es el menor carácter del alfabeto → lexicográficamente correcto.
static std::string with_sentinel(const std::string& text) {
    std::string s = text;
    s += '\0';
    return s;
}

/// Devuelve las posiciones en el texto donde `pattern` ocurre (solapamientos incluidos).
static std::vector<size_t> find_all(const std::string& text, const std::string& pattern) {
    std::vector<size_t> positions;
    size_t pos = 0;
    while ((pos = text.find(pattern, pos)) != std::string::npos) {
        positions.push_back(pos);
        ++pos;
    }
    return positions;
}

/**
 * Para un split concreto P[0..i-1] | P[i..m-1] (i = 1..m-1):
 * retorna el número de boundaries k cruzados que tienen
 *   start_{k+1} - occ_pos == i   (la ocurrencia cae en ese split exacto)
 */
static size_t brute_force_split_count(const std::string& text,
                                       const std::string& pattern,
                                       const LZ77Parsing& phrases,
                                       size_t split_i) {
    const auto occurrences = find_all(text, pattern);
    size_t count = 0;
    for (size_t occ : occurrences) {
        for (size_t k = 0; k + 1 < phrases.size(); ++k) {
            const size_t boundary = phrases[k + 1].start_pos;
            if (boundary == occ + split_i) count++;
        }
    }
    return count;
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests de construcción
// ─────────────────────────────────────────────────────────────────────────────

TEST(Grid2D_Build, SingleChar) {
    const std::string text = with_sentinel("a");
    LZ77Parser parser;
    const auto phrases = parser.parse(text);
    // "a\0" → 2 frases (literal 'a', luego literal '\0')
    // z=2 → z-1=1 punto si hay dos frases, o z=1 → 0 puntos si solo una
    // El parser ve "a\0": depende de la implementación; en cualquier caso z >= 1.
    Grid2D grid;
    grid.build(phrases, text);
    // Si z == 1: 0 puntos. Si z == 2: 1 punto.
    EXPECT_EQ(grid.point_count(), phrases.size() <= 1 ? 0u : phrases.size() - 1);
}

TEST(Grid2D_Build, Abracadabra) {
    const std::string text = with_sentinel("abracadabra");
    LZ77Parser parser;
    const auto phrases = parser.parse(text);
    ASSERT_GE(phrases.size(), 2u) << "Se esperan al menos 2 frases";

    Grid2D grid;
    grid.build(phrases, text);

    EXPECT_EQ(grid.point_count(), phrases.size() - 1);
}

TEST(Grid2D_Build, Repetitive) {
    const std::string text = with_sentinel("ACGTACGTACGTACGT");
    LZ77Parser parser;
    const auto phrases = parser.parse(text);
    ASSERT_GE(phrases.size(), 2u);

    Grid2D grid;
    grid.build(phrases, text);

    EXPECT_EQ(grid.point_count(), phrases.size() - 1);
}

TEST(Grid2D_Build, AllLiterals) {
    // Cada carácter es único → z = n frases (todas literales)
    const std::string text = with_sentinel("abcdefgh");
    LZ77Parser parser;
    const auto phrases = parser.parse(text);

    Grid2D grid;
    grid.build(phrases, text);

    EXPECT_EQ(grid.point_count(), phrases.size() - 1);
}

// ─────────────────────────────────────────────────────────────────────────────
// Invariantes de las coordenadas
// ─────────────────────────────────────────────────────────────────────────────

TEST(Grid2D_Invariants, CoordsDistinct_Abracadabra) {
    const std::string text = with_sentinel("abracadabra");
    LZ77Parser parser;
    const auto phrases = parser.parse(text);

    Grid2D grid;
    grid.build(phrases, text);

    // bv_fwd_ debe tener exactamente point_count() unos
    size_t ones_fwd = 0;
    for (size_t i = 0; i < grid.bv_fwd().size(); ++i)
        ones_fwd += grid.bv_fwd()[i];
    EXPECT_EQ(ones_fwd, grid.point_count());

    size_t ones_rev = 0;
    for (size_t i = 0; i < grid.bv_rev().size(); ++i)
        ones_rev += grid.bv_rev()[i];
    EXPECT_EQ(ones_rev, grid.point_count());
}

TEST(Grid2D_Invariants, RisPermutation) {
    // R debe ser una permutación de {0..z-2}
    const std::string text = with_sentinel("abracadabra");
    LZ77Parser parser;
    const auto phrases = parser.parse(text);
    if (phrases.size() <= 1) GTEST_SKIP();

    Grid2D grid;
    grid.build(phrases, text);

    const size_t z1 = grid.point_count();
    std::vector<bool> seen(z1, false);
    for (size_t j = 0; j < z1; ++j) {
        const size_t val = grid.wt()[j];
        ASSERT_LT(val, z1) << "Valor fuera de rango en R[" << j << "]";
        ASSERT_FALSE(seen[val]) << "Valor duplicado en R: " << val;
        seen[val] = true;
    }
}

TEST(Grid2D_Invariants, RisPermutation_Repetitive) {
    const std::string text = with_sentinel("ACGTACGTACGTACGTACGT");
    LZ77Parser parser;
    const auto phrases = parser.parse(text);
    if (phrases.size() <= 1) GTEST_SKIP();

    Grid2D grid;
    grid.build(phrases, text);

    const size_t z1 = grid.point_count();
    std::vector<bool> seen(z1, false);
    for (size_t j = 0; j < z1; ++j) {
        const size_t val = grid.wt()[j];
        ASSERT_LT(val, z1);
        ASSERT_FALSE(seen[val]);
        seen[val] = true;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests de consulta: brute-force vs query()
// ─────────────────────────────────────────────────────────────────────────────

/**
 * Para un texto y un patrón, itera sobre todos los splits i (1..m-1) y compara
 * el resultado de grid.query() con brute_force_split_count().
 *
 * Para poder llamar a grid.query() necesitamos rangos BWT: los obtenemos
 * construyendo un csa_wt de sdsl como ground truth.
 */
static void check_query_vs_brute_force(const std::string& text_raw,
                                        const std::string& pattern) {
    const std::string text = with_sentinel(text_raw);
    const size_t n = text.size();
    const size_t m = pattern.size();
    if (m < 2) return;  // necesitamos al menos split en i=1

    LZ77Parser parser;
    const auto phrases = parser.parse(text);
    if (phrases.size() <= 1) return;

    Grid2D grid;
    grid.build(phrases, text);

    // Construir csa_wt sobre el texto SIN sentinel para obtener rangos BWT.
    // sdsl::construct_im añade '\0' internamente → el SA resultante es idéntico
    // al SA de "text_raw\0" que usa la grilla (construido con divsufsort sobre
    // text_with_sentinel). Los rangos BWT son directamente comparables.
    sdsl::csa_wt<> csa_fwd;
    sdsl::construct_im(csa_fwd, text_raw, 1);  // 1 = byte alphabet

    const std::string text_raw_rev(text_raw.rbegin(), text_raw.rend());
    sdsl::csa_wt<> csa_rev;
    sdsl::construct_im(csa_rev, text_raw_rev, 1);

    // Para cada split i = 1..m-1:
    //   right_part = pattern[i..m-1]
    //   left_part_rev = reverse(pattern[0..i-1])
    for (size_t i = 1; i < m; ++i) {
        const std::string right_part(pattern.begin() + i, pattern.end());
        const std::string left_rev(pattern.rbegin() + (m - i), pattern.rend());

        // backward_search sobre el SA forward para right_part
        size_t sp_r = 0, ep_r = n - 1;
        sdsl::backward_search(csa_fwd, 0, n - 1,
                              right_part.begin(), right_part.end(),
                              sp_r, ep_r);

        // backward_search sobre el SA reverso para left_rev
        size_t sp_l = 0, ep_l = n - 1;
        sdsl::backward_search(csa_rev, 0, n - 1,
                              left_rev.begin(), left_rev.end(),
                              sp_l, ep_l);

        const size_t expected = brute_force_split_count(text_raw, pattern, phrases, i);

        if (sp_r > ep_r || sp_l > ep_l) {
            // El patrón no existe en ese split → 0 ocurrencias
            EXPECT_EQ(0u, expected)
                << "split=" << i << " pattern=" << pattern;
            continue;
        }

        const auto [count, pts] = grid.query(sp_r, ep_r, sp_l, ep_l);
        EXPECT_EQ(count, expected)
            << "split=" << i
            << " pattern=" << pattern
            << " right=[" << sp_r << "," << ep_r << "]"
            << " left=[" << sp_l << "," << ep_l << "]";
    }
}

TEST(Grid2D_Query, Abracadabra_Patterns) {
    const std::string text = "abracadabra";
    for (const std::string p : {"abr", "ra", "bra", "cada", "abra", "a"}) {
        SCOPED_TRACE("pattern=" + p);
        check_query_vs_brute_force(text, p);
    }
}

TEST(Grid2D_Query, Repetitive_DNA) {
    const std::string text = "ACGTACGTACGTACGT";
    for (const std::string p : {"ACGT", "CGTA", "GT", "ACGTACGT"}) {
        SCOPED_TRACE("pattern=" + p);
        check_query_vs_brute_force(text, p);
    }
}

TEST(Grid2D_Query, DNA_Synthetic) {
    const std::string text = "AAACCCGGGTTTTAAACCC";
    for (const std::string p : {"AAA", "CCC", "AAACCC", "GGG", "TTT"}) {
        SCOPED_TRACE("pattern=" + p);
        check_query_vs_brute_force(text, p);
    }
}

TEST(Grid2D_Query, NoMatch) {
    const std::string text = "abracadabra";
    // Patrón que no existe en el texto
    const std::string pattern = "xyz";
    const std::string text_s = with_sentinel(text);
    LZ77Parser parser;
    const auto phrases = parser.parse(text_s);
    if (phrases.size() <= 1) GTEST_SKIP();

    Grid2D grid;
    grid.build(phrases, text_s);

    // Rangos vacíos (sp > ep) → query debe devolver 0
    const auto [count, pts] = grid.query(5, 3, 5, 3);
    EXPECT_EQ(count, 0u);
    EXPECT_TRUE(pts.empty());
}

TEST(Grid2D_Query, AllSameChar) {
    // Texto muy compresible: pocas frases, muchos boundaries cruzados
    const std::string text = "aaaaaaaaaaaaaaaa";  // 16 'a'
    for (const std::string p : {"aa", "aaa", "aaaa"}) {
        SCOPED_TRACE("pattern=" + p);
        check_query_vs_brute_force(text, p);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests de text_pos: invariante de consistencia
// ─────────────────────────────────────────────────────────────────────────────

TEST(Grid2D_Extremal, TextPosConsistency) {
    // Para cada punto del WT, text_pos(j) debe ser start_{k+1} de alguna frase
    const std::string text_raw = "abracadabra";
    const std::string text = with_sentinel(text_raw);
    LZ77Parser parser;
    const auto phrases = parser.parse(text);
    if (phrases.size() <= 1) GTEST_SKIP();

    Grid2D grid;
    grid.build(phrases, text);

    // Recopilar el conjunto de start_{k+1} reales
    std::set<size_t> expected_boundaries;
    for (size_t k = 0; k + 1 < phrases.size(); ++k)
        expected_boundaries.insert(phrases[k + 1].start_pos);

    // Cada text_pos(j) debe estar en ese conjunto
    for (size_t j = 0; j < grid.point_count(); ++j) {
        EXPECT_TRUE(expected_boundaries.count(grid.text_pos(j)) > 0)
            << "text_pos(" << j << ") = " << grid.text_pos(j)
            << " no es start de ninguna frase";
    }
    // Y todos los boundaries intermedios (no el último) deben aparecer
    EXPECT_EQ(grid.point_count(), expected_boundaries.size());
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests de query_min_2d: invariante query_min_2d.boundary_min == brute-force min
// ─────────────────────────────────────────────────────────────────────────────

/**
 * Para cada split de cada patrón, verifica que query_min_2d retorna el mismo
 * count y el mismo boundary_min que la enumeración brute-force vía query().
 */
static void check_min_matches_enumeration(const std::string& text_raw,
                                           const std::string& pattern) {
    const std::string text = with_sentinel(text_raw);
    const size_t m = pattern.size();
    if (m < 2) return;

    LZ77Parser parser;
    const auto phrases = parser.parse(text);
    if (phrases.size() <= 1) return;

    Grid2D grid;
    grid.build(phrases, text);

    sdsl::csa_wt<> csa_fwd;
    sdsl::construct_im(csa_fwd, text_raw, 1);
    const std::string text_rev(text_raw.rbegin(), text_raw.rend());
    sdsl::csa_wt<> csa_rev;
    sdsl::construct_im(csa_rev, text_rev, 1);

    const size_t n = csa_fwd.size();

    for (size_t i = 1; i < m; ++i) {
        const std::string right_part(pattern.begin() + i, pattern.end());
        const std::string left_rev(pattern.rbegin() + (m - i), pattern.rend());

        size_t sp_r = 0, ep_r = n - 1;
        sdsl::backward_search(csa_fwd, 0, n - 1, right_part.begin(), right_part.end(), sp_r, ep_r);
        size_t sp_l = 0, ep_l = n - 1;
        sdsl::backward_search(csa_rev, 0, n - 1, left_rev.begin(), left_rev.end(), sp_l, ep_l);

        if (sp_r > ep_r || sp_l > ep_l) continue;

        // Brute-force via query() enumeration
        const auto [cnt, pts] = grid.query(sp_r, ep_r, sp_l, ep_l);
        size_t exp_bmin = SIZE_MAX;
        for (const auto& [wt_idx, y_rel] : pts) {
            const size_t bp = grid.text_pos(wt_idx);
            if (bp < exp_bmin) exp_bmin = bp;
        }

        const auto mn = grid.query_min_2d(sp_r, ep_r, sp_l, ep_l);

        EXPECT_EQ(mn.count, cnt)
            << "count mismatch split=" << i << " pattern=" << pattern;
        if (cnt > 0) {
            EXPECT_EQ(mn.boundary_min, exp_bmin)
                << "boundary_min mismatch split=" << i << " pattern=" << pattern;
        }
    }
}

TEST(Grid2D_MinRmq, MinMatchesEnumeration_Abracadabra) {
    const std::string text = "abracadabra";
    for (const std::string p : {"abr", "ra", "bra", "cada", "abra", "ab"}) {
        SCOPED_TRACE("pattern=" + p);
        check_min_matches_enumeration(text, p);
    }
}

TEST(Grid2D_MinRmq, MinMatchesEnumeration_Repetitive) {
    const std::string text = "ACGTACGTACGTACGT";
    for (const std::string p : {"ACGT", "CGTA", "GT", "ACGTACGT"}) {
        SCOPED_TRACE("pattern=" + p);
        check_min_matches_enumeration(text, p);
    }
}

TEST(Grid2D_MinRmq, MinMatchesEnumeration_AllSameChar) {
    const std::string text = "aaaaaaaaaaaaaaaa";
    for (const std::string p : {"aa", "aaa", "aaaa"}) {
        SCOPED_TRACE("pattern=" + p);
        check_min_matches_enumeration(text, p);
    }
}

TEST(Grid2D_MinRmq, MinMatchesEnumeration_DNASynthetic) {
    const std::string text = "AAACCCGGGTTTTAAACCC";
    for (const std::string p : {"AAA", "CCC", "AAACCC", "GGG", "TTT"}) {
        SCOPED_TRACE("pattern=" + p);
        check_min_matches_enumeration(text, p);
    }
}

TEST(Grid2D_MinRmq, EmptyGrid) {
    // Texto de 1 frase → grilla vacía
    const std::string text = with_sentinel("a");
    LZ77Parser parser;
    const auto phrases = parser.parse(text);
    Grid2D grid;
    grid.build(phrases, text);

    const auto mn = grid.query_min_2d(0, 10, 0, 10);
    EXPECT_EQ(mn.count, 0u);
    EXPECT_EQ(mn.boundary_min, SIZE_MAX);
}

// ── Grid2D_MaxRmq: query_max_2d == brute-force max via query() ───────────────

static void check_max_matches_enumeration(const std::string& text_raw,
                                           const std::string& pattern) {
    const std::string text = with_sentinel(text_raw);
    const size_t m = pattern.size();
    if (m < 2) return;

    LZ77Parser parser;
    const auto phrases = parser.parse(text);
    if (phrases.size() <= 1) return;

    Grid2D grid;
    grid.build(phrases, text);

    sdsl::csa_wt<> csa_fwd;
    sdsl::construct_im(csa_fwd, text_raw, 1);
    const std::string text_rev(text_raw.rbegin(), text_raw.rend());
    sdsl::csa_wt<> csa_rev;
    sdsl::construct_im(csa_rev, text_rev, 1);

    const size_t n = csa_fwd.size();

    for (size_t i = 1; i < m; ++i) {
        const std::string right_part(pattern.begin() + i, pattern.end());
        const std::string left_rev(pattern.rbegin() + (m - i), pattern.rend());

        size_t sp_r = 0, ep_r = n - 1;
        sdsl::backward_search(csa_fwd, 0, n - 1, right_part.begin(), right_part.end(), sp_r, ep_r);
        size_t sp_l = 0, ep_l = n - 1;
        sdsl::backward_search(csa_rev, 0, n - 1, left_rev.begin(), left_rev.end(), sp_l, ep_l);

        if (sp_r > ep_r || sp_l > ep_l) continue;

        // Brute-force via query() enumeration
        const auto [cnt, pts] = grid.query(sp_r, ep_r, sp_l, ep_l);
        size_t exp_bmax = 0;
        for (const auto& [wt_idx, y_rel] : pts) {
            const size_t bp = grid.text_pos(wt_idx);
            if (bp > exp_bmax) exp_bmax = bp;
        }

        const auto mx = grid.query_max_2d(sp_r, ep_r, sp_l, ep_l);

        EXPECT_EQ(mx.count, cnt)
            << "count mismatch split=" << i << " pattern=" << pattern;
        if (cnt > 0) {
            EXPECT_EQ(mx.boundary_max, exp_bmax)
                << "boundary_max mismatch split=" << i << " pattern=" << pattern;
        }
    }
}

TEST(Grid2D_MaxRmq, MaxMatchesEnumeration_Abracadabra) {
    const std::string text = "abracadabra";
    for (const std::string p : {"abr", "ra", "bra", "cada", "abra", "ab"}) {
        SCOPED_TRACE("pattern=" + p);
        check_max_matches_enumeration(text, p);
    }
}

TEST(Grid2D_MaxRmq, MaxMatchesEnumeration_Repetitive) {
    const std::string text = "ACGTACGTACGTACGT";
    for (const std::string p : {"ACGT", "CGTA", "GT", "ACGTACGT"}) {
        SCOPED_TRACE("pattern=" + p);
        check_max_matches_enumeration(text, p);
    }
}

TEST(Grid2D_MaxRmq, MaxMatchesEnumeration_AllSameChar) {
    const std::string text = "aaaaaaaaaaaaaaaa";
    for (const std::string p : {"aa", "aaa", "aaaa"}) {
        SCOPED_TRACE("pattern=" + p);
        check_max_matches_enumeration(text, p);
    }
}

TEST(Grid2D_MaxRmq, MaxMatchesEnumeration_DNASynthetic) {
    const std::string text = "AAACCCGGGTTTTAAACCC";
    for (const std::string p : {"AAA", "CCC", "AAACCC", "GGG", "TTT"}) {
        SCOPED_TRACE("pattern=" + p);
        check_max_matches_enumeration(text, p);
    }
}

TEST(Grid2D_MaxRmq, EmptyGrid) {
    const std::string text = with_sentinel("a");
    LZ77Parser parser;
    const auto phrases = parser.parse(text);
    Grid2D grid;
    grid.build(phrases, text);

    const auto mx = grid.query_max_2d(0, 10, 0, 10);
    EXPECT_EQ(mx.count, 0u);
    EXPECT_EQ(mx.boundary_max, 0u);
}

TEST(Grid2D_MaxRmq, MinLeqMax) {
    // Verificar que boundary_min <= boundary_max para todos los splits de cada patrón
    const std::string text_raw = "ACGTACGTACGTACGT";
    const std::string text = with_sentinel(text_raw);
    LZ77Parser parser;
    const auto phrases = parser.parse(text);
    Grid2D grid;
    grid.build(phrases, text);

    sdsl::csa_wt<> csa_fwd;
    sdsl::construct_im(csa_fwd, text_raw, 1);
    const std::string text_rev(text_raw.rbegin(), text_raw.rend());
    sdsl::csa_wt<> csa_rev;
    sdsl::construct_im(csa_rev, text_rev, 1);

    const size_t n = csa_fwd.size();

    for (const std::string pattern : {"ACGT", "CGTA", "GT", "ACGTACGT"}) {
        const size_t m = pattern.size();
        for (size_t i = 1; i < m; ++i) {
            size_t sp_r = 0, ep_r = n - 1;
            const std::string right_part(pattern.begin() + i, pattern.end());
            sdsl::backward_search(csa_fwd, 0, n-1, right_part.begin(), right_part.end(), sp_r, ep_r);
            size_t sp_l = 0, ep_l = n - 1;
            const std::string left_rev(pattern.rbegin() + (m - i), pattern.rend());
            sdsl::backward_search(csa_rev, 0, n-1, left_rev.begin(), left_rev.end(), sp_l, ep_l);
            if (sp_r > ep_r || sp_l > ep_l) continue;

            const auto mn = grid.query_min_2d(sp_r, ep_r, sp_l, ep_l);
            const auto mx = grid.query_max_2d(sp_r, ep_r, sp_l, ep_l);
            ASSERT_EQ(mn.count, mx.count);
            if (mn.count > 0) {
                EXPECT_LE(mn.boundary_min, mx.boundary_max)
                    << "min > max para split=" << i << " pattern=" << pattern;
            }
        }
    }
}
