#include <gtest/gtest.h>

#include "index.hpp"

#include <climits>
#include <filesystem>
#include <random>
#include <string>
#include <vector>

using namespace lz77tax;
namespace fs = std::filesystem;

// ─── helpers ─────────────────────────────────────────────────────────────────

static std::string build_text(const std::string& unit, int reps) {
    std::string t;
    for (int i = 0; i < reps; ++i) t += unit;
    return t;
}

static fs::path tmp_prefix(const std::string& label) {
    return fs::temp_directory_path() / ("lz77rt_" + label);
}

// Borra archivos generados por save() tras el test.
static void cleanup(const fs::path& prefix) {
    for (const auto& ext : {".meta", ".grid"})
        fs::remove(prefix.string() + ext);
}

// ─── roundtrip genérico ───────────────────────────────────────────────────────

static void check_roundtrip(const std::string& text, const std::string& label,
                              const std::vector<std::string>& patterns) {
    const fs::path pre = tmp_prefix(label);
    // Construir y guardar
    LZ77Index orig;
    orig.build(text);
    orig.save(pre);

    // Cargar en instancia nueva
    LZ77Index loaded;
    loaded.load(pre, text);

    cleanup(pre);

    // Comparar metadatos
    EXPECT_EQ(loaded.text_size(),   orig.text_size());
    EXPECT_EQ(loaded.phrase_count(), orig.phrase_count());
    EXPECT_EQ(loaded.grid_points(),  orig.grid_points());

    // Comparar queries para cada patrón
    for (const auto& p : patterns) {
        SCOPED_TRACE("pattern=" + p);

        EXPECT_EQ(loaded.count(p), orig.count(p));

        const auto [lo_o, hi_o] = orig.locate_extremal(p);
        const auto [lo_l, hi_l] = loaded.locate_extremal(p);
        EXPECT_EQ(lo_l, lo_o);
        EXPECT_EQ(hi_l, hi_o);
    }
}

// ─── tests ───────────────────────────────────────────────────────────────────

TEST(IndexRoundtrip, Abracadabra) {
    const std::string text = "abracadabra";
    check_roundtrip(text, "abra",
        {"abr", "ra", "bra", "cada", "abra", "ab", "xyz"});
}

TEST(IndexRoundtrip, RepetitiveDNA) {
    const std::string text = build_text("ACGTACGTTTACGT", 10);
    check_roundtrip(text, "dna_rep",
        {"ACGT", "TTAC", "CGT", "ACGTA", "ZZZZZ"});
}

TEST(IndexRoundtrip, AllLiterals) {
    const std::string text = "abcdefghijklm";
    check_roundtrip(text, "literals",
        {"ab", "cd", "ef", "abc", "xyz"});
}

TEST(IndexRoundtrip, Random100Patterns) {
    const std::string unit = "ACGTTTACGTAAACGTGG";
    const std::string text = build_text(unit, 8);
    const fs::path pre = tmp_prefix("rnd100");

    LZ77Index orig;
    orig.build(text);
    orig.save(pre);

    LZ77Index loaded;
    loaded.load(pre, text);
    cleanup(pre);

    std::mt19937 rng(99);
    int mismatches = 0;
    for (int t = 0; t < 100; ++t) {
        const size_t plen = 2 + rng() % 8;
        if (plen >= text.size()) continue;
        const size_t pos = rng() % (text.size() - plen + 1);
        const std::string p = text.substr(pos, plen);

        if (loaded.count(p) != orig.count(p)) ++mismatches;
        if (loaded.locate_extremal(p) != orig.locate_extremal(p)) ++mismatches;
    }
    EXPECT_EQ(mismatches, 0);
}

TEST(IndexRoundtrip, EmptyIndex) {
    // Texto de una sola frase: grilla vacía. save/load debe funcionar igual.
    const std::string text = "a";
    const fs::path pre = tmp_prefix("empty");

    LZ77Index orig;
    orig.build(text);
    orig.save(pre);

    LZ77Index loaded;
    loaded.load(pre, text);
    cleanup(pre);

    EXPECT_EQ(loaded.text_size(),    orig.text_size());
    EXPECT_EQ(loaded.phrase_count(), orig.phrase_count());
    EXPECT_EQ(loaded.grid_points(),  orig.grid_points());
    EXPECT_EQ(loaded.count("ab"),    0u);
}
