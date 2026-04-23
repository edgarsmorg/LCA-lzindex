#include <gtest/gtest.h>

#include "index.hpp"
#include "lz77/parser.hpp"

#include <algorithm>
#include <climits>
#include <string>
#include <vector>

using namespace lz77tax;

// ─────────────────────────────────────────────────────────────────────────────
// Brute-force helper
// ─────────────────────────────────────────────────────────────────────────────

/**
 * Cuenta pares (ocurrencia p, boundary k) donde P cruza un límite de frase.
 * Una ocurrencia en p cruza boundary b = phrases[k+1].start_pos si:
 *   p < b  &&  b <= p + m - 1
 * (equivalentemente: b ∈ (p, p+m-1], i.e., el boundary cae dentro del patrón,
 *  no en el primer ni el último carácter porque la posición b es el inicio de la
 *  siguiente frase — el carácter text[p] es el prefijo de la frase, no el inicio).
 *
 * NOTA: las frases están parseadas sobre text_s = text + '\0'. Las boundaries
 * con start_pos >= text.size() (la frase del centinela) no se cuentan.
 */
static size_t bf_primary_count(const std::string& text,
                                const std::string& pattern,
                                const LZ77Parsing& phrases) {
    const size_t m = pattern.size();
    if (m == 0) return 0;

    size_t count = 0;
    for (size_t p = 0; p + m <= text.size(); ++p) {
        if (text.compare(p, m, pattern) != 0) continue;
        for (size_t k = 0; k + 1 < phrases.size(); ++k) {
            const size_t b = phrases[k + 1].start_pos;
            if (b >= text.size()) continue;  // boundary de centinela
            if (p < b && b <= p + m - 1) count++;
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

    // Brute-force sobre el mismo parsing interno
    const size_t expected = bf_primary_count(text, pattern, idx.phrases());
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
 * Brute-force para locate_extremal: retorna {pos_min, pos_max} de las
 * posiciones de inicio de las ocurrencias primarias (las que cruzan al menos
 * un boundary de frase) del patrón en el texto.
 */
static std::pair<size_t, size_t> bf_locate_extremal(
        const std::string& text,
        const std::string& pattern,
        const LZ77Parsing& phrases) {
    const size_t m = pattern.size();
    size_t pos_min = SIZE_MAX, pos_max = 0;

    for (size_t p = 0; p + m <= text.size(); ++p) {
        if (text.compare(p, m, pattern) != 0) continue;
        for (size_t k = 0; k + 1 < phrases.size(); ++k) {
            const size_t b = phrases[k + 1].start_pos;
            if (b >= text.size()) continue;  // boundary del centinela
            if (p < b && b <= p + m - 1) {
                // Ocurrencia primaria en p
                if (p < pos_min) pos_min = p;
                if (p > pos_max) pos_max = p;
                break;  // basta un boundary cruzado para que sea primaria
            }
        }
    }
    return {pos_min, pos_max};
}

static void check_locate_extremal(const std::string& text,
                                   const std::string& pattern) {
    LZ77Index idx;
    idx.build(text);

    const auto [got_min, got_max] = idx.locate_extremal(pattern);
    const auto [exp_min, exp_max] = bf_locate_extremal(text, pattern, idx.phrases());

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
// locate_min() == locate_extremal().first
// ─────────────────────────────────────────────────────────────────────────────

static void check_locate_min(const std::string& text, const std::string& pattern) {
    LZ77Index idx;
    idx.build(text);

    const size_t got_min       = idx.locate_min(pattern);
    const auto [exp_min, exp_max] = idx.locate_extremal(pattern);

    EXPECT_EQ(got_min, exp_min)
        << "text=" << text << " pattern=" << pattern;
}

TEST(LZ77Index_LocateMin, Abracadabra) {
    const std::string text = "abracadabra";
    for (const std::string p : {"abr", "ra", "bra", "cada", "abra", "ab"}) {
        SCOPED_TRACE("pattern=" + p);
        check_locate_min(text, p);
    }
}

TEST(LZ77Index_LocateMin, Repetitive_DNA) {
    const std::string text = "ACGTACGTACGTACGT";
    for (const std::string p : {"ACGT", "CGTA", "GT", "ACGTACGT"}) {
        SCOPED_TRACE("pattern=" + p);
        check_locate_min(text, p);
    }
}

TEST(LZ77Index_LocateMin, AllSameChar) {
    const std::string text = "aaaaaaaaaaaaaaaa";
    for (const std::string p : {"aa", "aaa", "aaaa"}) {
        SCOPED_TRACE("pattern=" + p);
        check_locate_min(text, p);
    }
}

TEST(LZ77Index_LocateMin, DNA_Synthetic) {
    const std::string text = "AAACCCGGGTTTTAAACCC";
    for (const std::string p : {"AAA", "CCC", "AAACCC", "GGG"}) {
        SCOPED_TRACE("pattern=" + p);
        check_locate_min(text, p);
    }
}

TEST(LZ77Index_LocateMin, NoMatch) {
    LZ77Index idx;
    idx.build("abracadabra");
    EXPECT_EQ(idx.locate_min("xyz"), SIZE_MAX);
}

TEST(LZ77Index_LocateMin, ShortPatterns) {
    LZ77Index idx;
    idx.build("abracadabra");
    EXPECT_EQ(idx.locate_min(""),  SIZE_MAX);
    EXPECT_EQ(idx.locate_min("a"), SIZE_MAX);
}

TEST(LZ77Index_LocateMin, ConsistencyWithCount) {
    const std::string text = "ACGTACGTACGTACGT";
    LZ77Index idx;
    idx.build(text);

    for (const std::string p : {"ACGT", "GT", "xyz", "ACGTACGT", "ZZ"}) {
        SCOPED_TRACE("pattern=" + p);
        const size_t cnt = idx.count(p);
        const size_t mn  = idx.locate_min(p);
        if (cnt > 0) {
            EXPECT_NE(mn, SIZE_MAX) << "count>0 pero locate_min retornó vacío";
        } else {
            EXPECT_EQ(mn, SIZE_MAX) << "count=0 pero locate_min retornó algo";
        }
    }
}
