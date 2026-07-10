#include <gtest/gtest.h>

#include "index.hpp"
#include "lz77/parser.hpp"

#include <algorithm>
#include <climits>
#include <random>
#include <string>
#include <vector>

using namespace lz77tax;

// ─────────────────────────────────────────────────────────────────────────────
// Brute-force helper
// ─────────────────────────────────────────────────────────────────────────────

/**
 * Cuenta pares (ocurrencia p, boundary k) primarios: crossings + end-aligned.
 *
 * Crossing: b = phrases[k+1].start_pos cae estrictamente dentro del patrón
 *   → p < b && b <= p+m-1
 *
 * End-aligned (special): P termina en el trailing_char de la frase k
 *   → p + m - 1 == end_k  (= phrases[k].start_pos + phrases[k].length)
 *   → phrase_total_len_k >= m  (P cabe en la frase)
 *
 * Cada condición genera exactamente un par; no hay solapamiento entre ellas
 * para el mismo boundary (crossings requieren b dentro de P; special requiere
 * b = p+m, que está justo después).
 */
static size_t bf_primary_count(const std::string& text,
                                const std::string& pattern,
                                const LZ77Parsing& phrases) {
    const size_t m = pattern.size();
    // El índice retorna 0 para m < 2 (mismo contrato que count()).
    if (m < 2) return 0;

    size_t count = 0;
    for (size_t p = 0; p + m <= text.size(); ++p) {
        if (text.compare(p, m, pattern) != 0) continue;

        // Crossings: p starts within phrase k (i = b-p <= phrase_total_len_k)
        for (size_t k = 0; k + 1 < phrases.size(); ++k) {
            const size_t b = phrases[k + 1].start_pos;
            if (b >= text.size()) continue;  // boundary de centinela
            if (p < b && b <= p + m - 1 && b - p <= phrases[k].length + 1) count++;
        }

        // End-aligned (special): P ends at end_k
        for (size_t k = 0; k + 1 < phrases.size(); ++k) {
            const size_t end_k = phrases[k].start_pos + phrases[k].length;
            if (end_k >= text.size()) continue;  // frase centinela
            if (phrases[k].length + 1 >= m && p + m - 1 == end_k) count++;
        }
    }
    return count;
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests de build + count
// ─────────────────────────────────────────────────────────────────────────────

TEST(LZ77Index_Build, Abracadabra) {
    LZ77Index idx;
    idx.build("abracadabra");

    EXPECT_GE(idx.phrase_count(), 2u);
    EXPECT_EQ(idx.grid_points(), idx.phrase_count() - 1);
    // text_size incluye el centinela '\0'
    EXPECT_EQ(idx.text_size(), std::string("abracadabra").size() + 1);
}

TEST(LZ77Index_Build, Repetitive) {
    LZ77Index idx;
    idx.build("ACGTACGTACGTACGT");

    EXPECT_GE(idx.phrase_count(), 2u);
    EXPECT_EQ(idx.grid_points(), idx.phrase_count() - 1);
}

TEST(LZ77Index_Build, AllLiterals) {
    LZ77Index idx;
    idx.build("abcdefgh");

    EXPECT_GE(idx.phrase_count(), 2u);
    EXPECT_EQ(idx.grid_points(), idx.phrase_count() - 1);
}

TEST(LZ77Index_Build, SingleChar) {
    LZ77Index idx;
    idx.build("a");
    // Puede haber 1 o 2 frases según el parser; en cualquier caso la grilla es consistente
    EXPECT_EQ(idx.grid_points(), idx.phrase_count() <= 1u ? 0u : idx.phrase_count() - 1);
}

// ─────────────────────────────────────────────────────────────────────────────
// count() == brute-force
// ─────────────────────────────────────────────────────────────────────────────

static void check_count(const std::string& text,
                        const std::string& pattern) {
    LZ77Index idx;
    idx.build(text);

    const auto phrases_bf = LZ77Parser().parse(text + '\0');
    const size_t expected = bf_primary_count(text, pattern, phrases_bf);
    const size_t got      = idx.count(pattern);

    EXPECT_EQ(got, expected)
        << "text=" << text << " pattern=" << pattern;

    // Sanidad: primarias ≤ totales (no verificamos contra sdsl::count aquí
    // para no añadir una dependencia extra, pero la desigualdad debe cumplirse)
    EXPECT_LE(got, text.size())  // cota grosera
        << "count mayor que el largo del texto";
}

TEST(LZ77Index_Count, Abracadabra) {
    const std::string text = "abracadabra";
    for (const std::string p : {"abr", "ra", "bra", "cada", "abra", "a", "ab"}) {
        SCOPED_TRACE("pattern=" + p);
        check_count(text, p);
    }
}

TEST(LZ77Index_Count, Repetitive_DNA) {
    const std::string text = "ACGTACGTACGTACGT";
    for (const std::string p : {"ACGT", "CGTA", "GT", "ACGTACGT", "AC"}) {
        SCOPED_TRACE("pattern=" + p);
        check_count(text, p);
    }
}

TEST(LZ77Index_Count, DNA_Synthetic) {
    const std::string text = "AAACCCGGGTTTTAAACCC";
    for (const std::string p : {"AAA", "CCC", "AAACCC", "GGG", "TTT", "CCCGGG"}) {
        SCOPED_TRACE("pattern=" + p);
        check_count(text, p);
    }
}

TEST(LZ77Index_Count, AllSameChar) {
    const std::string text = "aaaaaaaaaaaaaaaa";  // 16 'a'
    for (const std::string p : {"aa", "aaa", "aaaa", "aaaaa"}) {
        SCOPED_TRACE("pattern=" + p);
        check_count(text, p);
    }
}

TEST(LZ77Index_Count, AllLiterals) {
    const std::string text = "abcdefghijklmno";
    for (const std::string p : {"ab", "bc", "cd", "abc", "bcd"}) {
        SCOPED_TRACE("pattern=" + p);
        check_count(text, p);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Patrones que no existen → count debe devolver 0
// ─────────────────────────────────────────────────────────────────────────────

TEST(LZ77Index_Count, NoMatch) {
    LZ77Index idx;
    idx.build("abracadabra");

    EXPECT_EQ(idx.count("xyz"), 0u);
    EXPECT_EQ(idx.count("zzz"), 0u);
    EXPECT_EQ(idx.count("abracadabraz"), 0u);
}

// ─────────────────────────────────────────────────────────────────────────────
// Patrones de longitud 0 o 1 → count siempre 0
// ─────────────────────────────────────────────────────────────────────────────

TEST(LZ77Index_Count, ShortPatterns) {
    LZ77Index idx;
    idx.build("abracadabra");

    EXPECT_EQ(idx.count(""), 0u);
    EXPECT_EQ(idx.count("a"), 0u);
    EXPECT_EQ(idx.count("b"), 0u);
}

// ─────────────────────────────────────────────────────────────────────────────
// Consistencia de accesores
// ─────────────────────────────────────────────────────────────────────────────

TEST(LZ77Index_Accessors, PhraseCountMatchesParser) {
    const std::string text = "abracadabra";
    const std::string text_s = text + '\0';

    LZ77Index idx;
    idx.build(text);

    LZ77Parser parser;
    const auto phrases_direct = parser.parse(text_s);

    EXPECT_EQ(idx.phrase_count(), phrases_direct.size());
}

TEST(LZ77Index_Accessors, GridPointsIsPhrasesMinusOne) {
    for (const std::string text : {"abracadabra", "ACGTACGTACGT", "aaaaaaaaa"}) {
        SCOPED_TRACE("text=" + text);
        LZ77Index idx;
        idx.build(text);
        if (idx.phrase_count() >= 2) {
            EXPECT_EQ(idx.grid_points(), idx.phrase_count() - 1);
        } else {
            EXPECT_EQ(idx.grid_points(), 0u);
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// locate_extremal()
// ─────────────────────────────────────────────────────────────────────────────

/**
 * Brute-force para locate_extremal: retorna {pos_min, pos_max} de las posiciones
 * de inicio de TODAS las ocurrencias del patrón en el texto.
 *
 * Con el diseño de dos índices (directo + reverso, ambos con RMQ de mínimo),
 * locate_extremal devuelve el mínimo y el máximo exactos sobre todas las
 * ocurrencias: la más a la izquierda siempre es primaria (índice directo) y la
 * más a la derecha se obtiene como la más a la izquierda de P^R (índice reverso).
 */
static std::pair<size_t, size_t> bf_locate_extremal(
        const std::string& text,
        const std::string& pattern) {
    const size_t m = pattern.size();
    size_t pos_min = SIZE_MAX, pos_max = 0;
    bool found = false;
    for (size_t p = 0; p + m <= text.size(); ++p) {
        if (text.compare(p, m, pattern) != 0) continue;
        if (p < pos_min) pos_min = p;
        if (p > pos_max) pos_max = p;
        found = true;
    }
    return found ? std::make_pair(pos_min, pos_max)
                 : std::make_pair(SIZE_MAX, size_t(0));
}

static void check_locate_extremal(const std::string& text,
                                   const std::string& pattern) {
    LZ77Index idx;
    idx.build(text);

    const auto [got_min, got_max] = idx.locate_extremal(pattern);
    const auto [exp_min, exp_max] = bf_locate_extremal(text, pattern);

    EXPECT_EQ(got_min, exp_min) << "pos_min  text=" << text << " pattern=" << pattern;
    EXPECT_EQ(got_max, exp_max) << "pos_max  text=" << text << " pattern=" << pattern;

    // Si se encontró algo, verificar que pos_min y pos_max son ocurrencias reales
    if (got_min != SIZE_MAX) {
        EXPECT_LE(got_min, got_max) << "pos_min > pos_max";
        EXPECT_EQ(text.compare(got_min, pattern.size(), pattern), 0)
            << "pos_min no es ocurrencia del patrón";
        EXPECT_EQ(text.compare(got_max, pattern.size(), pattern), 0)
            << "pos_max no es ocurrencia del patrón";
    }
}

TEST(LZ77Index_LocateExtremal, Abracadabra) {
    const std::string text = "abracadabra";
    for (const std::string p : {"abr", "ra", "bra", "cada", "abra", "ab"}) {
        SCOPED_TRACE("pattern=" + p);
        check_locate_extremal(text, p);
    }
}

TEST(LZ77Index_LocateExtremal, Repetitive_DNA) {
    const std::string text = "ACGTACGTACGTACGT";
    for (const std::string p : {"ACGT", "CGTA", "GT", "ACGTACGT"}) {
        SCOPED_TRACE("pattern=" + p);
        check_locate_extremal(text, p);
    }
}

TEST(LZ77Index_LocateExtremal, AllSameChar) {
    const std::string text = "aaaaaaaaaaaaaaaa";
    for (const std::string p : {"aa", "aaa", "aaaa"}) {
        SCOPED_TRACE("pattern=" + p);
        check_locate_extremal(text, p);
    }
}

TEST(LZ77Index_LocateExtremal, DNA_Synthetic) {
    const std::string text = "AAACCCGGGTTTTAAACCC";
    for (const std::string p : {"AAA", "CCC", "AAACCC", "GGG"}) {
        SCOPED_TRACE("pattern=" + p);
        check_locate_extremal(text, p);
    }
}

TEST(LZ77Index_LocateExtremal, NoMatch) {
    LZ77Index idx;
    idx.build("abracadabra");

    const auto [mn, mx] = idx.locate_extremal("xyz");
    EXPECT_EQ(mn, SIZE_MAX);
    EXPECT_EQ(mx, 0u);
}

TEST(LZ77Index_LocateExtremal, ShortPatterns) {
    LZ77Index idx;
    idx.build("abracadabra");

    // Patrones de longitud < 2 nunca tienen ocurrencias primarias
    EXPECT_EQ(idx.locate_extremal("").first,  SIZE_MAX);
    EXPECT_EQ(idx.locate_extremal("a").first,  SIZE_MAX);
}

TEST(LZ77Index_LocateExtremal, ConsistencyWithCount) {
    // Si count(P) > 0 → locate_extremal(P) debe devolver algo válido
    // Si count(P) == 0 → locate_extremal(P) debe devolver {SIZE_MAX, 0}
    const std::string text = "ACGTACGTACGTACGT";
    LZ77Index idx;
    idx.build(text);

    for (const std::string p : {"ACGT", "GT", "xyz", "ACGTACGT", "ZZ"}) {
        SCOPED_TRACE("pattern=" + p);
        const size_t cnt = idx.count(p);
        const auto [mn, mx] = idx.locate_extremal(p);
        if (cnt > 0) {
            EXPECT_NE(mn, SIZE_MAX) << "count>0 pero locate_extremal retornó vacío";
        } else {
            EXPECT_EQ(mn, SIZE_MAX) << "count=0 pero locate_extremal retornó algo";
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Task 5: Tests de correctitud con special occurrences
// ─────────────────────────────────────────────────────────────────────────────

// "aabaab" / "ab": la ocurrencia en p=1 es end-aligned (especial); la de p=4
// queda dentro de la frase k=2 (secundaria). count solo cuenta la primaria (1),
// pero locate_extremal reporta el máximo REAL (p=4) vía el índice reverso.
TEST(LZ77Index_SpecialCorrectness, SpecialOnlyText) {
    const std::string text    = "aabaab";
    const std::string pattern = "ab";

    LZ77Index idx;
    idx.build(text);

    const auto phrases_bf = LZ77Parser().parse(text + '\0');
    // Brute-force: 1 primaria (end-aligned en p=1); p=4 es secundaria.
    ASSERT_EQ(bf_primary_count(text, pattern, phrases_bf), 1u)
        << "bf_primary_count inesperado";

    EXPECT_EQ(idx.count(pattern), 1u);

    const auto [mn, mx] = idx.locate_extremal(pattern);
    EXPECT_EQ(mn, 1u);   // más a la izquierda (primaria)
    EXPECT_EQ(mx, 4u);   // más a la derecha (secundaria, vía índice reverso)
}

// Para textos con crossings Y specials, count == bf sin duplicados.
TEST(LZ77Index_SpecialCorrectness, NoDoubleCount) {
    for (const std::string text : {
            std::string("aabaab"),
            std::string("ACGTACGTACGT"),
            std::string("abracadabra"),
            std::string("AAACCCGGGTTTTAAACCC")}) {
        SCOPED_TRACE("text=" + text);
        LZ77Index idx;
        idx.build(text);
        const auto phrases_bf = LZ77Parser().parse(text + '\0');

        for (size_t plen = 2; plen <= std::min(text.size(), size_t(6)); ++plen) {
            for (size_t s = 0; s + plen <= text.size(); s += 3) {
                const std::string p = text.substr(s, plen);
                EXPECT_EQ(idx.count(p),
                          bf_primary_count(text, p, phrases_bf))
                    << "plen=" << plen << " s=" << s;
            }
        }
    }
}

// Si idx.count(P) > 0 entonces P existe en T (sdsl::count > 0).
// idx.count mide pares (ocurrencia, split_i), no posiciones únicas,
// por lo que puede superar sdsl::count — pero implica existencia real.
TEST(LZ77Index_SpecialCorrectness, CountImpliesExistence) {
    const std::string text = "abracadabra";
    sdsl::csa_wt<> csa;
    sdsl::construct_im(csa, text, 1);

    LZ77Index idx;
    idx.build(text);

    for (const std::string p : {"ab", "abr", "ra", "bra", "cada",
                                  "abra", "a", "x", "aa", "abracadabra"}) {
        SCOPED_TRACE("pattern=" + p);
        if (idx.count(p) > 0) {
            EXPECT_GT(sdsl::count(csa, p.begin(), p.end()), 0u)
                << "idx.count>0 pero sdsl::count=0 para pattern=" << p;
        }
    }
}

// 100 substrings aleatorios de texto repetitivo: count == bf_primary_count.
TEST(LZ77Index_SpecialCorrectness, GroundTruth100Patterns) {
    const std::string unit = "ACGTTTACGTAAACGT";
    std::string text;
    for (int i = 0; i < 8; ++i) text += unit;  // 128 chars

    LZ77Index idx;
    idx.build(text);
    const auto phrases_bf = LZ77Parser().parse(text + '\0');

    std::mt19937 rng(42);
    int mismatches = 0;
    for (int trial = 0; trial < 100; ++trial) {
        const size_t plen = 2 + rng() % 7;
        if (plen > text.size()) continue;
        const size_t pos = rng() % (text.size() - plen + 1);
        const std::string p = text.substr(pos, plen);

        const size_t got = idx.count(p);
        const size_t exp = bf_primary_count(text, p, phrases_bf);
        if (got != exp) {
            ++mismatches;
            ADD_FAILURE()
                << "trial=" << trial << " pos=" << pos
                << " plen=" << plen << " pattern=" << p
                << "  got=" << got << " exp=" << exp;
        }
    }
    EXPECT_EQ(mismatches, 0);
}