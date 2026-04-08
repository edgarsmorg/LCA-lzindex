#include <gtest/gtest.h>
#include "lz77/parser.hpp"
#include "lz77/phrase.hpp"

#include <string>
#include <numeric>

using namespace lz77tax;

// ── Helpers ───────────────────────────────────────────────────────────────────

/// Verifica que las frases parten el texto correctamente:
///   - phrases[0].start_pos == 0
///   - phrases[i+1].start_pos == phrases[i].start_pos + phrases[i].length + 1
///   - La suma total de (length + 1) == n
///   - next_char de cada frase coincide con text[start_pos + length]
static void check_partition(const std::string& text, const LZ77Parsing& phrases) {
    ASSERT_FALSE(phrases.empty()) << "Texto no vacío debe producir al menos una frase";
    EXPECT_EQ(phrases[0].start_pos, 0u);

    size_t pos = 0;
    for (size_t k = 0; k < phrases.size(); ++k) {
        EXPECT_EQ(phrases[k].start_pos, pos)
            << "Frase " << k << " empieza en " << phrases[k].start_pos
            << " pero esperábamos " << pos;

        EXPECT_EQ(phrases[k].next_char, (uint8_t)text[pos + phrases[k].length])
            << "next_char incorrecto en frase " << k;

        pos += phrases[k].length + 1;
    }
    EXPECT_EQ(pos, text.size()) << "Las frases no cubren el texto completo";
}

/// Verifica la propiedad LZ77: cada frase con length > 0 tiene una fuente en
/// j < start_pos tal que T[j..j+length-1] == T[start_pos..start_pos+length-1].
/// Las copias solapadas (j + length > start_pos) son VÁLIDAS en LZ77.
static void check_lz77_property(const std::string& text, const LZ77Parsing& phrases) {
    for (size_t k = 0; k < phrases.size(); ++k) {
        const auto& ph = phrases[k];
        if (ph.length == 0) continue;

        bool found = false;
        for (size_t j = 0; j < ph.start_pos && !found; ++j) {
            bool match = true;
            for (size_t l = 0; l < ph.length && match; ++l)
                if (text[j + l] != text[ph.start_pos + l]) match = false;
            if (match) found = true;
        }

        EXPECT_TRUE(found)
            << "Frase " << k << " en pos " << ph.start_pos
            << " con length=" << ph.length
            << ": no hay fuente previa válida (las copias solapadas son ok)";
    }
}

// ── Tests básicos ─────────────────────────────────────────────────────────────

TEST(LZ77Parser, EmptyString) {
    LZ77Parser parser;
    auto phrases = parser.parse("");
    EXPECT_TRUE(phrases.empty());
    EXPECT_EQ(parser.phrase_count(), 0u);
}

TEST(LZ77Parser, SingleChar) {
    LZ77Parser parser;
    auto phrases = parser.parse("a");
    ASSERT_EQ(phrases.size(), 1u);
    EXPECT_EQ(phrases[0].start_pos, 0u);
    EXPECT_EQ(phrases[0].length, 0u);
    EXPECT_EQ(phrases[0].next_char, (uint8_t)'a');
    check_partition("a", phrases);
}

TEST(LZ77Parser, AllUnique) {
    // Sin repeticiones → cada carácter es una frase literal → z = n
    LZ77Parser parser;
    const std::string text = "abcdefgh";
    auto phrases = parser.parse(text);
    EXPECT_EQ(phrases.size(), text.size());
    for (auto& ph : phrases) EXPECT_EQ(ph.length, 0u);
    check_partition(text, phrases);
    check_lz77_property(text, phrases);
}

TEST(LZ77Parser, AllSameChar) {
    // "aaaa" (n=4): frase 1 = literal 'a' (i=0),
    // frase 2 = copia de largo 2 + next_char 'a' → cubre i=1..3 → z=2
    LZ77Parser parser;
    const std::string text = "aaaa";
    auto phrases = parser.parse(text);
    EXPECT_EQ(phrases.size(), 2u);
    check_partition(text, phrases);
    check_lz77_property(text, phrases);
}

TEST(LZ77Parser, AllSameCharLonger) {
    LZ77Parser parser;
    const std::string text(20, 'a');
    auto phrases = parser.parse(text);
    EXPECT_LT(phrases.size(), text.size());  // debe haber compresión
    check_partition(text, phrases);
    check_lz77_property(text, phrases);
}

// ── Tests canónicos ───────────────────────────────────────────────────────────

TEST(LZ77Parser, Abracadabra) {
    // "abracadabra" (n=11)
    // Parsing esperado: a | b | r | (a,c) | (a,d) | (abr,a) → z=6
    LZ77Parser parser;
    const std::string text = "abracadabra";
    auto phrases = parser.parse(text);
    EXPECT_EQ(phrases.size(), 6u);
    check_partition(text, phrases);
    check_lz77_property(text, phrases);
}

TEST(LZ77Parser, AcgtRepeated) {
    // "ACGTACGTACGT" — patrón de 4 repetido 3 veces (ADN sintético)
    LZ77Parser parser;
    const std::string text = "ACGTACGTACGT";
    auto phrases = parser.parse(text);
    // Primera ocurrencia: 4 literales (A,C,G,T), luego copia de 8 chars
    // z esperado ≤ 6 (hay compresión significativa)
    EXPECT_LT(phrases.size(), text.size());
    check_partition(text, phrases);
    check_lz77_property(text, phrases);
}

// ── Test sobre el texto sintético real ────────────────────────────────────────

TEST(LZ77Parser, SyntheticReference) {
    // El texto sintético de 85 chars tiene repeticiones (ADN con 4 genomas)
    const std::string text =
        "AAAAACGTACGTACGTACGT#ACGTACGTACGTGGGGGGGG#GGGGGGGGTTTTTTGCATGCA"
        "#TGCATGCATGCACCCCCCCC#";

    LZ77Parser parser;
    auto phrases = parser.parse(text);

    EXPECT_GT(phrases.size(), 0u);
    EXPECT_LT(phrases.size(), text.size());  // hay compresión
    check_partition(text, phrases);
    check_lz77_property(text, phrases);
}

// ── Test: phrase_count() es coherente ────────────────────────────────────────

TEST(LZ77Parser, PhraseCountConsistency) {
    LZ77Parser parser;
    const std::string text = "banana";
    auto phrases = parser.parse(text);
    EXPECT_EQ(parser.phrase_count(), phrases.size());
}

// ── Test: textos con separador '#' (como los datos reales) ───────────────────

TEST(LZ77Parser, WithSeparator) {
    const std::string text = "ACGT#TGCA#AAAA#";
    LZ77Parser parser;
    auto phrases = parser.parse(text);
    check_partition(text, phrases);
    check_lz77_property(text, phrases);
}
