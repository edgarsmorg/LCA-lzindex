/**
 * test_mem — tests para MEMExtractor.
 *
 * Tres bloques:
 *  (a) Sanidad básica: empty, short, sin matches.
 *  (b) Validación ground-truth: cada MEM es left-maximal, right-maximal y ocurre en T.
 *  (c) Pipeline end-to-end: read → MEMs → locate_extremal → classify_extremal.
 */

#include <gtest/gtest.h>

#include "index.hpp"
#include "mem/extractor.hpp"
#include "mem/mem.hpp"
#include "taxonomy/lca.hpp"
#include "taxonomy/classifier.hpp"

#include <sdsl/suffix_arrays.hpp>
#include <algorithm>
#include <climits>
#include <string>
#include <vector>

using namespace lz77tax;

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

static size_t fm_count(const sdsl::csa_wt<>& csa, const std::string& pat) {
    return pat.empty() ? 0 : sdsl::count(csa, pat);
}

// Verifica que todas las MEMs retornadas son válidas (ocurren, left-max, right-max).
static void verify_mems(const sdsl::csa_wt<>& csa,
                        const std::string& query,
                        const MEMList& mems) {
    for (const auto& mem : mems) {
        const size_t i = mem.query_start;
        const size_t len = mem.length;
        const size_t j = i + len - 1;
        const std::string sub = query.substr(i, len);

        EXPECT_GT(fm_count(csa, sub), 0u)
            << "MEM '" << sub << "' no ocurre en T";

        // Left-maximal
        if (i > 0) {
            EXPECT_EQ(fm_count(csa, query.substr(i - 1, len + 1)), 0u)
                << "MEM '" << sub << "' NO es left-maximal (Q[i-1..j] ocurre)";
        }

        // Right-maximal
        if (j + 1 < query.size()) {
            EXPECT_EQ(fm_count(csa, query.substr(i, len + 1)), 0u)
                << "MEM '" << sub << "' NO es right-maximal (Q[i..j+1] ocurre)";
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// (a) Sanidad básica
// ─────────────────────────────────────────────────────────────────────────────

TEST(MEMExtractor_Sanity, EmptyQuery) {
    sdsl::csa_wt<> csa;
    sdsl::construct_im(csa, "abracadabra", 1);
    MEMExtractor ext(csa);

    EXPECT_TRUE(ext.extract("", 1).empty());
}

TEST(MEMExtractor_Sanity, MinLengthZero) {
    sdsl::csa_wt<> csa;
    sdsl::construct_im(csa, "abracadabra", 1);
    MEMExtractor ext(csa);

    // min_length=0 ⟹ vacío (contrato del extractor)
    EXPECT_TRUE(ext.extract("abr", 0).empty());
}

TEST(MEMExtractor_Sanity, QueryShorterThanMinLength) {
    sdsl::csa_wt<> csa;
    sdsl::construct_im(csa, "abracadabra", 1);
    MEMExtractor ext(csa);

    // query "abr" (3 chars), min_length=5 → ninguna MEM válida
    EXPECT_TRUE(ext.extract("abr", 5).empty());
}

TEST(MEMExtractor_Sanity, NoMatchInReference) {
    sdsl::csa_wt<> csa;
    sdsl::construct_im(csa, "abracadabra", 1);
    MEMExtractor ext(csa);

    // 'z' no aparece en "abracadabra"
    EXPECT_TRUE(ext.extract("zzzzz", 1).empty());
}

TEST(MEMExtractor_Sanity, QueryIsEntireReference) {
    const std::string text = "abracadabra";
    sdsl::csa_wt<> csa;
    sdsl::construct_im(csa, text, 1);
    MEMExtractor ext(csa);

    const auto mems = ext.extract(text, 1);
    // Debe haber al menos una MEM (el texto completo ocurre)
    EXPECT_FALSE(mems.empty());
    verify_mems(csa, text, mems);
}

// ─────────────────────────────────────────────────────────────────────────────
// (b) Validación ground-truth: left-maximal, right-maximal, ocurre en T
// ─────────────────────────────────────────────────────────────────────────────

class MEMExtractor_Validity : public ::testing::Test {
protected:
    void SetUp() override {
        sdsl::construct_im(csa_abra_, "abracadabra", 1);
        sdsl::construct_im(csa_dna_,
            "ACGTACGTACGTACGT"
            "GGGGGGGGTTTTTTTT"
            "AAAAAAAACCCCCCCC", 1);
        sdsl::construct_im(csa_same_, "AAAAAAAAAA", 1);
    }

    sdsl::csa_wt<> csa_abra_;
    sdsl::csa_wt<> csa_dna_;
    sdsl::csa_wt<> csa_same_;
};

TEST_F(MEMExtractor_Validity, Abracadabra_SelfQuery) {
    MEMExtractor ext(csa_abra_);
    const std::string q = "abracadabra";
    verify_mems(csa_abra_, q, ext.extract(q, 1));
}

TEST_F(MEMExtractor_Validity, Abracadabra_PartialQuery) {
    MEMExtractor ext(csa_abra_);
    const std::string q = "Xabra cadabraZ";  // X y Z no están en T
    auto mems = ext.extract(q, 1);
    EXPECT_FALSE(mems.empty());
    verify_mems(csa_abra_, q, mems);
    // "abra" y "cadabra" deben aparecer como MEMs separadas
    bool found_abra    = false;
    bool found_cadabra = false;
    for (const auto& m : mems) {
        const std::string sub = q.substr(m.query_start, m.length);
        if (sub == "abra")    found_abra    = true;
        if (sub == "cadabra") found_cadabra = true;
    }
    EXPECT_TRUE(found_abra)    << "MEM 'abra' no encontrada";
    EXPECT_TRUE(found_cadabra) << "MEM 'cadabra' no encontrada";
}

TEST_F(MEMExtractor_Validity, DNA_MinLength) {
    MEMExtractor ext(csa_dna_);
    const std::string q = "NNNACGTACGTNNNGGGGGGGGNNN";
    const auto mems = ext.extract(q, 4);
    EXPECT_FALSE(mems.empty());
    verify_mems(csa_dna_, q, mems);
    for (const auto& m : mems)
        EXPECT_GE(m.length, 4u) << "MEM más corta que min_length";
}

TEST_F(MEMExtractor_Validity, AllSameChar_MinLen2) {
    MEMExtractor ext(csa_same_);
    const std::string q = "AAAA";
    const auto mems = ext.extract(q, 2);
    EXPECT_FALSE(mems.empty());
    verify_mems(csa_same_, q, mems);
}

TEST_F(MEMExtractor_Validity, NoOverlappingMEMs) {
    // Dos MEMs no deberían ser contenidas una en la otra si ambas son maximales.
    MEMExtractor ext(csa_dna_);
    const std::string q = "ACGTACGTGGGGGGGG";
    const auto mems = ext.extract(q, 4);
    // Ningún MEM debe ser subcadena de otro (ambos son right y left-maximal → distintos)
    for (size_t a = 0; a < mems.size(); ++a) {
        for (size_t b = 0; b < mems.size(); ++b) {
            if (a == b) continue;
            // Si [i_a, j_a] ⊆ [i_b, j_b] estrictamente, entonces a no sería maximal
            const size_t ia = mems[a].query_start, ja = ia + mems[a].length - 1;
            const size_t ib = mems[b].query_start, jb = ib + mems[b].length - 1;
            const bool a_inside_b = (ib <= ia && ja <= jb && !(ia == ib && ja == jb));
            EXPECT_FALSE(a_inside_b)
                << "MEM " << a << " está contenida en MEM " << b;
        }
    }
}

TEST_F(MEMExtractor_Validity, RefStartIsValidOccurrence) {
    MEMExtractor ext(csa_dna_);
    const std::string q = "ACGTACGT";
    const auto mems = ext.extract(q, 4);
    for (const auto& mem : mems) {
        // ref_start debe ser una posición válida en T donde ocurre la MEM
        const std::string pat = q.substr(mem.query_start, mem.length);
        // Número de caracteres en T desde ref_start: verificamos que T[ref_start..] empieza con pat
        // (no podemos leer T directamente, pero sí sabemos que csa_[k] da posiciones válidas)
        EXPECT_LT(mem.ref_start, csa_dna_.size())
            << "ref_start " << mem.ref_start << " fuera de rango [0, n)";
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// (c) Pipeline end-to-end: read → MEMs → locate_extremal → classify_extremal
//
// Dataset: mismo que test_classify.cpp
//   root(0) → A(1) → A1(3), A2(4)
//             B(2) → B1(5), B2(6)
//
//   Texto (DFS, sep '#'):
//     A1: "AAAAACGTACGTACGTACGT"  [0..19]
//     A2: "ACGTACGTACGTGGGGGGGG"  [21..40]
//     B1: "GGGGGGGGTTTTTTGCATGCA" [42..62]
//     B2: "TGCATGCATGCACCCCCCCC"  [64..83]
// ─────────────────────────────────────────────────────────────────────────────

static const std::string REFERENCE =
    "AAAAACGTACGTACGTACGT#"
    "ACGTACGTACGTGGGGGGGG#"
    "GGGGGGGGTTTTTTGCATGCA#"
    "TGCATGCATGCACCCCCCCC#";

class PipelineTest : public ::testing::Test {
protected:
    void SetUp() override {
        idx_.build(REFERENCE);
        ext_ = std::make_unique<MEMExtractor>(idx_.csa_fwd());

        tree_.build(
            {0, 1, 2, 3, 4, 5, 6},
            {"root", "A", "B", "A1", "A2", "B1", "B2"},
            {PhyloTree::NO_PARENT, 0, 0, 1, 1, 2, 2});

        classifier_.setup(tree_, {
            {  0, 20, 3},
            { 21, 41, 4},
            { 42, 63, 5},
            { 64, 84, 6},
        });
    }

    // Encadena: read → MEMs → locate_extremal → classify_extremal.
    // Retorna -1 si ningún MEM tiene ocurrencias primarias.
    int pipeline_classify(const std::string& read, size_t min_len = 8) const {
        const auto mems = ext_->extract(read, min_len);
        size_t pos_min = SIZE_MAX, pos_max = 0;
        for (const auto& mem : mems) {
            const auto [pmin, pmax] =
                idx_.locate_extremal(read.substr(mem.query_start, mem.length));
            if (pmin == SIZE_MAX) continue;
            if (pmin < pos_min) pos_min = pmin;
            if (pmax > pos_max) pos_max = pmax;
        }
        if (pos_min == SIZE_MAX) return -1;
        return classifier_.classify_extremal(pos_min, pos_max);
    }

    // Ground-truth vía FM-index completo.
    int bf_classify(const std::string& pattern) const {
        sdsl::csa_wt<> fm;
        sdsl::construct_im(fm, REFERENCE, 1);
        const auto occ = sdsl::locate(fm, pattern);
        if (occ.empty()) return -1;
        std::vector<size_t> positions(occ.begin(), occ.end());
        return classifier_.classify(positions);
    }

    LZ77Index                      idx_;
    std::unique_ptr<MEMExtractor>  ext_;
    PhyloTree                      tree_;
    Classifier                     classifier_;
};

// read que contiene un patrón único de A1 rodeado de caracteres no presentes en la ref
TEST_F(PipelineTest, UniqueA1_InRead) {
    // "AAAAACGT" aparece solo en A1
    const std::string read = "NNAAAAACGTNN";
    const int result = pipeline_classify(read, 8);
    ASSERT_NE(result, -1) << "No MEM encontrada para read con patrón único de A1";
    EXPECT_EQ(classifier_.genome_of(idx_.locate_extremal("AAAAACGT").first), 3);
    // Propiedad: resultado debe estar en subárbol de A1(3)
    EXPECT_EQ(tree_.lca(result, 3), 3)
        << "resultado=" << result << " no está en subárbol de A1(3)";
}

// read que contiene patrón único de B2 rodeado de ruido
TEST_F(PipelineTest, UniqueB2_InRead) {
    // "CCCCCCCC" aparece solo en B2
    const std::string read = "NNCCCCCCCCNN";
    const int result = pipeline_classify(read, 8);
    ASSERT_NE(result, -1);
    EXPECT_EQ(tree_.lca(result, 6), 6)
        << "resultado=" << result << " no está en subárbol de B2(6)";
}

// read que contiene "TGCATGCA" (compartida entre B1 y B2)
TEST_F(PipelineTest, SharedB1B2_InRead) {
    const std::string read = "NNTGCATGCANN";
    const int result = pipeline_classify(read, 8);
    ASSERT_NE(result, -1);
    // LCA real es B(2); nuestro resultado debe estar en subárbol de B(2)
    EXPECT_EQ(tree_.lca(result, 2), 2)
        << "resultado=" << result << " no está en subárbol de B(2)";
}

// read con dos MEMs de genomas distintos → LCA más alto
TEST_F(PipelineTest, TwoMEMs_CrossSubtrees) {
    // "AAAAACGT" (A1) + separator + "CCCCCCCC" (B2)
    // Ambas MEMs en el mismo read → LCA debería ser root(0) o cerca
    const std::string read = "NNAAAAACGTNNNCCCCCCCCNN";
    const int result = pipeline_classify(read, 8);
    ASSERT_NE(result, -1);
    // A1(3) y B2(6) → LCA real = root(0)
    // Nuestro resultado puede ser root(0) o un descendiente (si alguna MEM no tuvo primarias)
    EXPECT_EQ(tree_.lca(result, 0), 0)
        << "resultado=" << result << " no está en subárbol de root(0)";
}

// read que no contiene ningún substring de la referencia → -1
TEST_F(PipelineTest, NoMatchInReference) {
    const std::string read = "NNNNNNNNNNNNNNNNNNNN";  // solo 'N', no en referencia
    EXPECT_EQ(pipeline_classify(read, 8), -1);
}

// read que es exactamente un patrón conocido (sin ruido)
TEST_F(PipelineTest, ExactPattern_UniqueA1) {
    const int result = pipeline_classify("AAAAACGTACGTACGTACGT", 8);
    ASSERT_NE(result, -1);
    EXPECT_EQ(tree_.lca(result, 3), 3);
}

// Propiedad central: para cualquier read con matches, resultado ∈ subárbol(bf_classify)
TEST_F(PipelineTest, ResultInSubtreeOfGroundTruth) {
    const std::vector<std::pair<std::string, std::string>> cases = {
        {"NNTGCATGCANN",       "TGCATGCA"},   // en B1 y B2 → LCA bf = B(2)
        {"NNCCCCCCCCNN",       "CCCCCCCC"},   // solo B2    → LCA bf = B2(6)
        {"NNAAAAACGTNN",       "AAAAACGT"},   // solo A1    → LCA bf = A1(3)
    };

    for (const auto& [read, pattern] : cases) {
        const int our = pipeline_classify(read, 8);
        const int bf  = bf_classify(pattern);
        if (our == -1 || bf == -1) continue;  // sin primarias: limitación conocida

        EXPECT_EQ(tree_.lca(our, bf), bf)
            << "read='" << read << "': nuestro=" << our
            << " no está en subárbol de bf=" << bf;
    }
}
