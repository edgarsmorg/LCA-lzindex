/**
 * test_baseline_sr — tests del pipeline LCA usando SrIndexLocator como baseline.
 *
 * Dataset sintético idéntico al de test_classify.cpp:
 *
 *        root(0)
 *       /       \
 *     A(1)      B(2)
 *    /    \    /    \
 *  A1(3) A2(4) B1(5) B2(6)
 *
 *  Texto (DFS order, separador '#'):
 *    A1: "AAAAACGTACGTACGTACGT"  [0..19]
 *    A2: "ACGTACGTACGTGGGGGGGG"  [21..40]
 *    B1: "GGGGGGGGTTTTTTGCATGCA" [42..62]
 *    B2: "TGCATGCATGCACCCCCCCC"  [64..83]
 *
 * Diferencia clave con test_classify.cpp:
 *   SrIndexLocator ve TODAS las ocurrencias (primarias + secundarias),
 *   por lo que el LCA obtenido debe ser EXACTAMENTE igual al ground truth
 *   del FM-index, sin la limitación primary-only del LZ77-index.
 */

#include <gtest/gtest.h>

#include "baseline/sr_index_locator.hpp"
#include "taxonomy/lca.hpp"
#include "taxonomy/classifier.hpp"
#include "mem/extractor.hpp"
#include "test_helpers.hpp"

#include <sdsl/suffix_arrays.hpp>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <climits>
#include <memory>
#include <string>
#include <vector>

using namespace lz77tax;
namespace fs = std::filesystem;

static const std::string REFERENCE =
    "AAAAACGTACGTACGTACGT#"
    "ACGTACGTACGTGGGGGGGG#"
    "GGGGGGGGTTTTTTGCATGCA#"
    "TGCATGCATGCACCCCCCCC#";

// ─────────────────────────────────────────────────────────────────────────────
// Fixture
// ─────────────────────────────────────────────────────────────────────────────

class BaselineSrTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Construir sr-index sobre el texto sintético.
        // Escribimos a un archivo temporal, construimos, y limpiamos al final.
        tmp_dir_ = fs::temp_directory_path() / ("sr_baseline_test_" + std::to_string(getpid()));
        fs::create_directories(tmp_dir_);

        const fs::path data_file = tmp_dir_ / "reference.txt";
        {
            std::ofstream f(data_file);
            f << REFERENCE;
        }

        sr_loc_.build(data_file.string(), tmp_dir_.string(), /*sr=*/4);

        // FM-index para ground truth
        sdsl::construct_im(fm_, REFERENCE, 1);
        fm_ext_ = std::make_unique<MEMExtractor>(fm_);

        // Árbol filogenético
        tree_ = lz77tax::test::make_toy_tree();
        classifier_.setup(tree_, lz77tax::test::make_toy_genome_ranges());
    }

    void TearDown() override {
        fs::remove_all(tmp_dir_);
    }

    // Ground truth: FM-index enumera todas las ocurrencias, toma min/max → LCA.
    int bf_classify(const std::string& pattern) const {
        const auto occ = sdsl::locate(fm_, pattern);
        if (occ.empty()) return -1;
        std::vector<size_t> positions(occ.begin(), occ.end());
        return classifier_.classify(positions);
    }

    // Ground truth: min/max de todas las ocurrencias via FM-index.
    std::pair<size_t,size_t> bf_extremal(const std::string& pattern) const {
        const auto occ = sdsl::locate(fm_, pattern);
        if (occ.empty()) return {SIZE_MAX, 0};
        return {*std::min_element(occ.begin(), occ.end()),
                *std::max_element(occ.begin(), occ.end())};
    }

    bool is_in_subtree(int descendant, int ancestor) const {
        return tree_.lca(descendant, ancestor) == ancestor;
    }

    fs::path                       tmp_dir_;
    SrIndexLocator                 sr_loc_;
    PhyloTree                      tree_;
    Classifier                     classifier_;
    sdsl::csa_wt<>                 fm_;
    std::unique_ptr<MEMExtractor>  fm_ext_;
};

// ─────────────────────────────────────────────────────────────────────────────
// Tests de locate_extremal: sr-index debe ser exacto (no approximate)
// ─────────────────────────────────────────────────────────────────────────────

TEST_F(BaselineSrTest, LocateExtremal_ExactMatchesFMIndex) {
    const std::vector<std::string> patterns = {
        "AAAAACGT", "ACGTACGT", "TGCATGCA", "GGGGGGGG", "CCCCCCCC"
    };
    for (const auto& p : patterns) {
        SCOPED_TRACE("pattern=" + p);
        const auto [sr_min, sr_max] = sr_loc_.locate_extremal(p);
        const auto [bf_min, bf_max] = bf_extremal(p);
        EXPECT_EQ(sr_min, bf_min) << "pos_min difiere del ground truth";
        EXPECT_EQ(sr_max, bf_max) << "pos_max difiere del ground truth";
    }
}

TEST_F(BaselineSrTest, LocateExtremal_NoMatch_ReturnsMaxSize) {
    const auto [mn, mx] = sr_loc_.locate_extremal("NNNNNNNN");
    EXPECT_EQ(mn, SIZE_MAX);
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests de LCA exacto (propiedad FUERTE vs la débil del LZ77-index)
// ─────────────────────────────────────────────────────────────────────────────

TEST_F(BaselineSrTest, LCA_Exact_UniqueA1) {
    const auto [mn, mx] = sr_loc_.locate_extremal("AAAAACGT");
    ASSERT_NE(mn, SIZE_MAX);
    EXPECT_EQ(classifier_.classify_extremal(mn, mx), 3) << "AAAAACGT solo en A1";
}

TEST_F(BaselineSrTest, LCA_Exact_UniqueB2) {
    const auto [mn, mx] = sr_loc_.locate_extremal("CCCCCCCC");
    ASSERT_NE(mn, SIZE_MAX);
    EXPECT_EQ(classifier_.classify_extremal(mn, mx), 6) << "CCCCCCCC solo en B2";
}

TEST_F(BaselineSrTest, LCA_Exact_SharedA1A2) {
    // "ACGTACGT" aparece en A1 y A2 → LCA exacto = A(1).
    // (El LZ77-index tiene limitación primary-only aquí y retorna A1(3)).
    const auto [mn, mx] = sr_loc_.locate_extremal("ACGTACGT");
    ASSERT_NE(mn, SIZE_MAX);
    EXPECT_EQ(classifier_.classify_extremal(mn, mx), bf_classify("ACGTACGT"));
}

TEST_F(BaselineSrTest, LCA_Exact_SharedB1B2) {
    const auto [mn, mx] = sr_loc_.locate_extremal("TGCATGCA");
    ASSERT_NE(mn, SIZE_MAX);
    EXPECT_EQ(classifier_.classify_extremal(mn, mx), 2) << "TGCATGCA en B1+B2 → B";
}

TEST_F(BaselineSrTest, LCA_Exact_SharedA2B1) {
    // "GGGGGGGG" aparece en A2 y B1 → LCA exacto = root(0).
    // (El LZ77-index tiene limitación primary-only y puede retornar A2(4)).
    const auto [mn, mx] = sr_loc_.locate_extremal("GGGGGGGG");
    ASSERT_NE(mn, SIZE_MAX);
    EXPECT_EQ(classifier_.classify_extremal(mn, mx), 0) << "GGGGGGGG en A2+B1 → root";
}

// Verificación exhaustiva: sr-index == FM-index para todos los patrones de prueba.
TEST_F(BaselineSrTest, LCA_ExactEqualsGroundTruth_AllPatterns) {
    const std::vector<std::pair<std::string,int>> cases = {
        {"ACGTACGT", 1},   // A1+A2 → A
        {"TGCATGCA", 2},   // B1+B2 → B
        {"GGGGGGGG", 0},   // A2+B1 → root
        {"AAAAACGT", 3},   // solo A1
        {"CCCCCCCC", 6},   // solo B2
    };
    for (const auto& [pat, true_lca] : cases) {
        SCOPED_TRACE("pattern=" + pat);
        const int bf = bf_classify(pat);
        ASSERT_EQ(bf, true_lca) << "bf_classify incorrecto";

        const auto [mn, mx] = sr_loc_.locate_extremal(pat);
        ASSERT_NE(mn, SIZE_MAX);
        const int sr_lca = classifier_.classify_extremal(mn, mx);
        EXPECT_EQ(sr_lca, true_lca) << "sr-index debería dar LCA exacto";
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Pipeline completa: classify_read con SrIndexLocator (template)
// ─────────────────────────────────────────────────────────────────────────────

TEST_F(BaselineSrTest, ClassifyRead_UniqueA1_Exact) {
    const int result = classifier_.classify_read("AAAAACGT", sr_loc_, *fm_ext_, 8);
    ASSERT_NE(result, -1);
    EXPECT_EQ(result, 3);
}

TEST_F(BaselineSrTest, ClassifyRead_UniqueB2_Exact) {
    const int result = classifier_.classify_read("CCCCCCCC", sr_loc_, *fm_ext_, 8);
    ASSERT_NE(result, -1);
    EXPECT_EQ(result, 6);
}

TEST_F(BaselineSrTest, ClassifyRead_SharedA1A2_Exact) {
    // sr-index debe retornar A(1) exacto, a diferencia del LZ77-index que retorna A1(3).
    const int result = classifier_.classify_read("ACGTACGT", sr_loc_, *fm_ext_, 8);
    ASSERT_NE(result, -1);
    EXPECT_EQ(result, bf_classify("ACGTACGT"));
}

TEST_F(BaselineSrTest, ClassifyRead_SharedA2B1_Exact) {
    const int result = classifier_.classify_read("GGGGGGGG", sr_loc_, *fm_ext_, 8);
    ASSERT_NE(result, -1);
    EXPECT_EQ(result, 0) << "GGGGGGGG en A2+B1 → root(0)";
}

TEST_F(BaselineSrTest, ClassifyRead_NoMatch) {
    EXPECT_EQ(classifier_.classify_read("NNNNNNNN", sr_loc_, *fm_ext_, 8), -1);
}

// ─────────────────────────────────────────────────────────────────────────────
// locate_count: debe coincidir con el count real del FM-index
// ─────────────────────────────────────────────────────────────────────────────

TEST_F(BaselineSrTest, LocateCount_MatchesFMIndex) {
    const std::vector<std::string> patterns = {
        "AAAAACGT", "ACGTACGT", "TGCATGCA", "GGGGGGGG", "CCCCCCCC"
    };
    for (const auto& p : patterns) {
        SCOPED_TRACE("pattern=" + p);
        const size_t sr_cnt = sr_loc_.locate_count(p);
        const size_t fm_cnt = sdsl::count(fm_, p);
        EXPECT_EQ(sr_cnt, fm_cnt);
    }
}
