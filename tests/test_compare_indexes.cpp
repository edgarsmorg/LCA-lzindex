/**
 * test_compare_indexes — verifica la equivalencia entre LZ77Index y SrIndexLocator
 * en el pipeline classify_read → LCA.
 *
 * Propiedad central verificada:
 *   Para cualquier patrón P, classify_read(P, lz_idx) está SIEMPRE en el
 *   subárbol de classify_read(P, sr_loc). El sr-index ve todas las ocurrencias
 *   → da el LCA real. El LZ77-index solo ve primarias → puede dar un descendiente.
 *   Nunca puede dar un nodo incomparable (eso sería un bug).
 *
 * Categorías de resultado para cada patrón:
 *   EQUAL              — ambos retornan el mismo nodo
 *   LZ_DESCENDANT_OF_SR — lz_result en subárbol de sr_result (limitación expected)
 *   INCOMPARABLE       — lz_result y sr_result en ramas distintas (bug)
 *
 * Dataset sintético (85 chars, 4 genomas, árbol de altura 2):
 *        root(0)
 *       /       \
 *     A(1)      B(2)
 *    /    \    /    \
 *  A1(3) A2(4) B1(5) B2(6)
 */

#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <vector>
#include <unistd.h>

#include "index.hpp"
#include "taxonomy/lca.hpp"
#include "taxonomy/classifier.hpp"
#include "mem/extractor.hpp"
#include "baseline/sr_index_locator.hpp"

#include <sdsl/suffix_arrays.hpp>

using namespace lz77tax;
namespace fs = std::filesystem;

static const std::string REFERENCE =
    "AAAAACGTACGTACGTACGT#"
    "ACGTACGTACGTGGGGGGGG#"
    "GGGGGGGGTTTTTTGCATGCA#"
    "TGCATGCATGCACCCCCCCC#";

class CompareIndexesTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Directorio temporal para archivos sdsl del sr-index
        tmp_dir_ = fs::temp_directory_path()
                 / ("compare_idx_test_" + std::to_string(getpid()));
        fs::create_directories(tmp_dir_);

        // Escribir referencia a disco para el sr-index
        const fs::path data_file = tmp_dir_ / "reference.txt";
        {
            std::ofstream f(data_file);
            f << REFERENCE;
        }

        // Construir LZ77Index (en memoria)
        lz_idx_.build(REFERENCE);

        // Construir SrIndexLocator (desde archivo)
        sr_loc_.build(data_file.string(), tmp_dir_.string(), /*sr=*/4);

        // Árbol filogenético
        tree_.build(
            {0, 1, 2, 3, 4, 5, 6},
            {"root", "A", "B", "A1", "A2", "B1", "B2"},
            {PhyloTree::NO_PARENT, 0, 0, 1, 1, 2, 2}
        );

        // Classifier
        classifier_.setup(tree_, {
            {  0, 21, 3},   // A1: [0..20] exclusive
            { 21, 42, 4},   // A2: [21..41] exclusive
            { 42, 64, 5},   // B1: [42..63] exclusive
            { 64, 85, 6},   // B2: [64..84] exclusive
        });

        // FM-index + MEMExtractor para ground-truth y classify_read
        sdsl::construct_im(fm_, REFERENCE, 1);
        ext_ = std::make_unique<MEMExtractor>(fm_);
    }

    void TearDown() override {
        fs::remove_all(tmp_dir_);
    }

    // Categoría de equivalencia entre LZ77 y sr-index
    enum class Category { EQUAL, LZ_DESCENDANT_OF_SR, INCOMPARABLE, EITHER_UNCLASSIFIED };

    Category compare(const std::string& pattern, int min_mem = 8) const {
        const int lz = classifier_.classify_read(pattern, lz_idx_, *ext_, min_mem);
        const int sr = classifier_.classify_read(pattern, sr_loc_, *ext_, min_mem);
        if (lz == -1 || sr == -1) return Category::EITHER_UNCLASSIFIED;
        if (lz == sr)             return Category::EQUAL;
        // lz está en subárbol de sr ↔ lca(lz,sr) == sr
        if (tree_.lca(lz, sr) == sr) return Category::LZ_DESCENDANT_OF_SR;
        return Category::INCOMPARABLE;
    }

    fs::path                       tmp_dir_;
    LZ77Index                      lz_idx_;
    SrIndexLocator                 sr_loc_;
    PhyloTree                      tree_;
    Classifier                     classifier_;
    sdsl::csa_wt<>                 fm_;
    std::unique_ptr<MEMExtractor>  ext_;
};

// ─────────────────────────────────────────────────────────────────────────────
// Propiedad fundamental: LZ77 nunca es INCOMPARABLE con sr-index.
// Puede ser EQUAL o LZ_DESCENDANT_OF_SR, pero nunca en una rama distinta.
// ─────────────────────────────────────────────────────────────────────────────

TEST_F(CompareIndexesTest, NeverIncomparable_UniqueA1) {
    EXPECT_NE(compare("AAAAACGT"), Category::INCOMPARABLE);
}

TEST_F(CompareIndexesTest, NeverIncomparable_UniqueB2) {
    EXPECT_NE(compare("CCCCCCCC"), Category::INCOMPARABLE);
}

TEST_F(CompareIndexesTest, NeverIncomparable_SharedB1B2) {
    EXPECT_NE(compare("TGCATGCA"), Category::INCOMPARABLE);
}

TEST_F(CompareIndexesTest, NeverIncomparable_SharedA1A2) {
    EXPECT_NE(compare("ACGTACGT"), Category::INCOMPARABLE);
}

TEST_F(CompareIndexesTest, NeverIncomparable_SharedA2B1) {
    EXPECT_NE(compare("GGGGGGGG"), Category::INCOMPARABLE);
}

// ─────────────────────────────────────────────────────────────────────────────
// Casos con resultado EXACT: ambos índices coinciden
// ─────────────────────────────────────────────────────────────────────────────

TEST_F(CompareIndexesTest, Equal_UniqueA1) {
    // "AAAAACGT" solo en A1 → ambos deben dar A1(3)
    EXPECT_EQ(compare("AAAAACGT"), Category::EQUAL);
}

TEST_F(CompareIndexesTest, Equal_UniqueB2) {
    EXPECT_EQ(compare("CCCCCCCC"), Category::EQUAL);
}

TEST_F(CompareIndexesTest, Equal_SharedB1B2) {
    // "TGCATGCA" cruza boundary B1/B2 → primarias visibles para LZ77
    EXPECT_EQ(compare("TGCATGCA"), Category::EQUAL);
}

// ─────────────────────────────────────────────────────────────────────────────
// Casos con limitación conocida: LZ77 da descendiente del real (primary-only)
// ─────────────────────────────────────────────────────────────────────────────

TEST_F(CompareIndexesTest, LzDescendant_SharedA1A2) {
    // "ACGTACGT" en A1 y A2. sr-index da A(1). LZ77 solo ve ocurrencias
    // que cruzan boundaries → puede dar A1(3) que es descendiente de A(1).
    const auto cat = compare("ACGTACGT");
    EXPECT_NE(cat, Category::INCOMPARABLE)
        << "LZ77 no puede estar en rama distinta a sr-index";
    // Documentar el resultado actual (puede ser EQUAL si hay primarias en A2)
    if (cat == Category::LZ_DESCENDANT_OF_SR) {
        SUCCEED() << "Limitación primary-only confirmada: LZ77 da descendiente de A(1)";
    } else if (cat == Category::EQUAL) {
        SUCCEED() << "LZ77 coincide con sr-index para ACGTACGT";
    }
}

TEST_F(CompareIndexesTest, LzDescendant_SharedA2B1) {
    const auto cat = compare("GGGGGGGG");
    EXPECT_NE(cat, Category::INCOMPARABLE)
        << "LZ77 no puede estar en rama distinta a sr-index";
}

// ─────────────────────────────────────────────────────────────────────────────
// Barrido exhaustivo: ningún patrón del vocabulario es INCOMPARABLE
// ─────────────────────────────────────────────────────────────────────────────

TEST_F(CompareIndexesTest, NeverIncomparable_AllPatterns) {
    const std::vector<std::string> patterns = {
        "AAAAACGT", "ACGTACGT", "TGCATGCA", "GGGGGGGG", "CCCCCCCC",
        "ACGTACGTACGT", "TGCATGCATGCA", "GCATGCATGCA",
    };
    for (const auto& p : patterns) {
        SCOPED_TRACE("pattern=" + p);
        EXPECT_NE(compare(p), Category::INCOMPARABLE)
            << "LZ77 en rama distinta a sr-index: bug en la implementación";
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Verificar que sr-index concuerda con FM-index (exactitud del baseline)
// ─────────────────────────────────────────────────────────────────────────────

TEST_F(CompareIndexesTest, SrIndex_ExactEqualsGroundTruth) {
    const std::vector<std::pair<std::string,int>> cases = {
        {"ACGTACGT",  1},   // A1+A2 → A
        {"TGCATGCA",  2},   // B1+B2 → B
        {"GGGGGGGG",  0},   // A2+B1 → root
        {"AAAAACGT",  3},   // solo A1
        {"CCCCCCCC",  6},   // solo B2
    };
    for (const auto& [pat, true_lca] : cases) {
        SCOPED_TRACE("pattern=" + pat);
        const int sr = classifier_.classify_read(pat, sr_loc_, *ext_, 8);
        ASSERT_NE(sr, -1);
        EXPECT_EQ(sr, true_lca)
            << "sr-index debe dar LCA exacto (sin limitación primary-only)";
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Verificar que cuando ambos clasifican, LZ77 es ancestro-o-igual del sr-result
// ─────────────────────────────────────────────────────────────────────────────

TEST_F(CompareIndexesTest, LzAlwaysInSubtreeOfSr) {
    const std::vector<std::string> patterns = {
        "AAAAACGT", "ACGTACGT", "TGCATGCA", "GGGGGGGG", "CCCCCCCC",
    };
    for (const auto& p : patterns) {
        SCOPED_TRACE("pattern=" + p);
        const int lz = classifier_.classify_read(p, lz_idx_, *ext_, 8);
        const int sr = classifier_.classify_read(p, sr_loc_, *ext_, 8);
        if (lz == -1 || sr == -1) continue;
        EXPECT_EQ(tree_.lca(lz, sr), sr)
            << "lz=" << lz << " no está en subárbol de sr=" << sr;
    }
}
