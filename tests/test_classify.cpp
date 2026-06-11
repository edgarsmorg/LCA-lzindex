/**
 * test_classify — tests end-to-end del pipeline LZ77Index → locate_extremal → LCA.
 *
 * Dataset sintético (85 chars, 4 genomas, árbol de altura 2):
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
 * LIMITACIÓN CONOCIDA (primary-only approach):
 *   locate_extremal solo encuentra ocurrencias que CRUZAN un boundary de frase.
 *   Si un patrón aparece en un genoma X solo dentro de frases copiadas desde otro
 *   genoma Y (ocurrencias secundarias), nuestra implementación no lo detecta en X.
 *   Resultado: el LCA puede ser más específico que el real (descendiente del correcto).
 *
 *   Propiedad garantizada: si find_extremal retorna (pos_min, pos_max) y classify
 *   retorna L, entonces L es el nodo real O un descendiente de él. Nunca retorna un
 *   nodo incorrecto en el sentido de "ancestro de la respuesta real".
 *
 *   Para mitigar: usar MEMs largos (≥31 bp) donde es improbable que toda ocurrencia
 *   en un genoma quede dentro de una sola frase LZ77.
 */

#include <gtest/gtest.h>

#include "index.hpp"
#include "taxonomy/lca.hpp"
#include "taxonomy/classifier.hpp"
#include "mem/extractor.hpp"
#include "test_helpers.hpp"

#include <sdsl/suffix_arrays.hpp>
#include <climits>
#include <memory>
#include <string>
#include <vector>

using namespace lz77tax;

static const std::string REFERENCE =
    "AAAAACGTACGTACGTACGT#"
    "ACGTACGTACGTGGGGGGGG#"
    "GGGGGGGGTTTTTTGCATGCA#"
    "TGCATGCATGCACCCCCCCC#";

// ─────────────────────────────────────────────────────────────────────────────
// Fixture compartido
// ─────────────────────────────────────────────────────────────────────────────

class ClassifyTest : public ::testing::Test {
protected:
    void SetUp() override {
        idx_.build(REFERENCE);

        tree_ = lz77tax::test::make_toy_tree();
        classifier_.setup(tree_, lz77tax::test::make_toy_genome_ranges());

        sdsl::construct_im(fm_, REFERENCE, 1);
        ext_ = std::make_unique<MEMExtractor>(fm_);
    }

    // Ground truth: localiza todas las ocurrencias via FM-index y clasifica.
    int bf_classify(const std::string& pattern) const {
        const auto occ = sdsl::locate(fm_, pattern);
        if (occ.empty()) return -1;
        std::vector<size_t> positions(occ.begin(), occ.end());
        return classifier_.classify(positions);
    }

    // Verifica que 'descendant' está en el subárbol de 'ancestor' (o es igual).
    bool is_in_subtree(int descendant, int ancestor) const {
        return tree_.lca(descendant, ancestor) == ancestor;
    }

    LZ77Index                     idx_;
    PhyloTree                     tree_;
    Classifier                    classifier_;
    sdsl::csa_wt<>                fm_;
    std::unique_ptr<MEMExtractor> ext_;
};

// ─────────────────────────────────────────────────────────────────────────────
// Tests unitarios del árbol
// ─────────────────────────────────────────────────────────────────────────────

TEST(PhyloTreeTest, Depths) {
    const PhyloTree t = lz77tax::test::make_toy_tree();

    EXPECT_EQ(t.depth(0), 0);
    EXPECT_EQ(t.depth(1), 1);
    EXPECT_EQ(t.depth(2), 1);
    EXPECT_EQ(t.depth(3), 2);
    EXPECT_EQ(t.depth(4), 2);
    EXPECT_EQ(t.depth(5), 2);
    EXPECT_EQ(t.depth(6), 2);
}

TEST(PhyloTreeTest, LCA_SameNode) {
    const PhyloTree t = lz77tax::test::make_toy_tree();

    EXPECT_EQ(t.lca(0, 0), 0);
    EXPECT_EQ(t.lca(3, 3), 3);
    EXPECT_EQ(t.lca(6, 6), 6);
}

TEST(PhyloTreeTest, LCA_ParentChild) {
    const PhyloTree t = lz77tax::test::make_toy_tree();

    EXPECT_EQ(t.lca(1, 3), 1);   // A y A1 → A
    EXPECT_EQ(t.lca(2, 5), 2);   // B y B1 → B
    EXPECT_EQ(t.lca(0, 6), 0);   // root y B2 → root
}

TEST(PhyloTreeTest, LCA_Siblings) {
    const PhyloTree t = lz77tax::test::make_toy_tree();

    EXPECT_EQ(t.lca(3, 4), 1);   // A1, A2 → A
    EXPECT_EQ(t.lca(5, 6), 2);   // B1, B2 → B
}

TEST(PhyloTreeTest, LCA_CrossSubtrees) {
    const PhyloTree t = lz77tax::test::make_toy_tree();

    EXPECT_EQ(t.lca(3, 5), 0);   // A1, B1 → root
    EXPECT_EQ(t.lca(4, 6), 0);   // A2, B2 → root
    EXPECT_EQ(t.lca(3, 6), 0);   // A1, B2 → root
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests del clasificador con posiciones conocidas
// ─────────────────────────────────────────────────────────────────────────────

TEST_F(ClassifyTest, GenomeOf) {
    EXPECT_EQ(classifier_.genome_of( 0), 3);   // A1 inicio
    EXPECT_EQ(classifier_.genome_of(19), 3);   // A1 fin
    EXPECT_EQ(classifier_.genome_of(20), -1);  // '#' A1/A2
    EXPECT_EQ(classifier_.genome_of(21), 4);   // A2 inicio
    EXPECT_EQ(classifier_.genome_of(40), 4);   // A2 fin
    EXPECT_EQ(classifier_.genome_of(41), -1);  // '#' A2/B1
    EXPECT_EQ(classifier_.genome_of(42), 5);   // B1 inicio
    EXPECT_EQ(classifier_.genome_of(62), 5);   // B1 fin
    EXPECT_EQ(classifier_.genome_of(63), -1);  // '#' B1/B2
    EXPECT_EQ(classifier_.genome_of(64), 6);   // B2 inicio
    EXPECT_EQ(classifier_.genome_of(83), 6);   // B2 fin
    EXPECT_EQ(classifier_.genome_of(84), -1);  // '#' final
}

TEST_F(ClassifyTest, ClassifyDirectPositions) {
    EXPECT_EQ(classifier_.classify({ 5, 25}), 1);   // A1+A2 → A
    EXPECT_EQ(classifier_.classify({55, 64}), 2);   // B1+B2 → B
    EXPECT_EQ(classifier_.classify({33, 42}), 0);   // A2+B1 → root
    EXPECT_EQ(classifier_.classify({ 0,  7}), 3);   // ambos A1 → A1
}

TEST_F(ClassifyTest, ClassifyExtremal_Equivalence) {
    EXPECT_EQ(classifier_.classify_extremal( 5, 25), 1);   // A1+A2 → A
    EXPECT_EQ(classifier_.classify_extremal(55, 64), 2);   // B1+B2 → B
    EXPECT_EQ(classifier_.classify_extremal(33, 42), 0);   // A2+B1 → root
    EXPECT_EQ(classifier_.classify_extremal( 0,  7), 3);   // ambos A1 → A1
}

// ─────────────────────────────────────────────────────────────────────────────
// Ground truth: FM-index brute-force classify
// ─────────────────────────────────────────────────────────────────────────────

TEST_F(ClassifyTest, BruteForce_SharedA1A2)   { EXPECT_EQ(bf_classify("ACGTACGT"), 1); }
TEST_F(ClassifyTest, BruteForce_SharedB1B2)   { EXPECT_EQ(bf_classify("TGCATGCA"), 2); }
TEST_F(ClassifyTest, BruteForce_SharedA2B1)   { EXPECT_EQ(bf_classify("GGGGGGGG"), 0); }
TEST_F(ClassifyTest, BruteForce_UniqueA1)     { EXPECT_EQ(bf_classify("AAAAACGT"), 3); }
TEST_F(ClassifyTest, BruteForce_UniqueB2)     { EXPECT_EQ(bf_classify("CCCCCCCC"), 6); }

// ─────────────────────────────────────────────────────────────────────────────
// Pipeline completa: locate_extremal → classify_extremal
// ─────────────────────────────────────────────────────────────────────────────

// Propiedad central garantizada: nuestro resultado es siempre un nodo dentro
// del subárbol del LCA real. Es decir, podemos ser MÁS ESPECÍFICOS pero nunca
// retornamos un nodo que esté "por encima" del LCA correcto.
//
// Formalmente: lca(nuestro_resultado, bf_result) == bf_result
//              ⟺ nuestro_resultado está en el subárbol de bf_result.
TEST_F(ClassifyTest, ExtremalIsInSubtreeOfTrueLCA) {
    const std::vector<std::pair<std::string,int>> cases = {
        {"ACGTACGT", 1},   // verdadero LCA = A
        {"TGCATGCA", 2},   // verdadero LCA = B
        {"GGGGGGGG", 0},   // verdadero LCA = root
        {"AAAAACGT", 3},   // verdadero LCA = A1
        {"CCCCCCCC", 6},   // verdadero LCA = B2
    };
    for (const auto& [pat, true_lca] : cases) {
        const int bf = bf_classify(pat);
        ASSERT_EQ(bf, true_lca) << "bf_classify incorrecto: " << pat;

        auto [pmin, pmax] = idx_.locate_extremal(pat);
        if (pmin == SIZE_MAX) {
            ADD_FAILURE() << "Sin ocurrencias primarias para: " << pat;
            continue;
        }
        const int our_lca = classifier_.classify_extremal(pmin, pmax);
        EXPECT_TRUE(is_in_subtree(our_lca, true_lca))
            << "patron='" << pat << "'  nuestro=" << our_lca
            << "  verdadero=" << true_lca
            << "  (nuestro debería estar en subárbol del verdadero)";
    }
}

// Casos donde locate_extremal encuentra primarias en ambos genomas relevantes
// → el LCA es EXACTAMENTE correcto (no solo un descendiente).
TEST_F(ClassifyTest, LCA_Exact_UniqueA1) {
    // "AAAAACGT" aparece solo en A1 → debe ser exacto
    auto [pmin, pmax] = idx_.locate_extremal("AAAAACGT");
    ASSERT_NE(pmin, SIZE_MAX);
    EXPECT_EQ(classifier_.classify_extremal(pmin, pmax), 3);
}

TEST_F(ClassifyTest, LCA_Exact_UniqueB2) {
    // "CCCCCCCC" aparece solo en B2 → debe ser exacto
    auto [pmin, pmax] = idx_.locate_extremal("CCCCCCCC");
    ASSERT_NE(pmin, SIZE_MAX);
    EXPECT_EQ(classifier_.classify_extremal(pmin, pmax), 6);
}

TEST_F(ClassifyTest, LCA_Exact_SharedB1B2) {
    // "TGCATGCA": el '#' entre B1 y B2 genera boundary que expone ambas como primarias
    auto [pmin, pmax] = idx_.locate_extremal("TGCATGCA");
    ASSERT_NE(pmin, SIZE_MAX);
    EXPECT_EQ(classifier_.classify_extremal(pmin, pmax), 2);
}

// Casos con limitación conocida: ocurrencias secundarias (copias dentro de frases)
// no son detectadas. El LCA retornado es un descendiente del real.
TEST_F(ClassifyTest, LCA_KnownLimitation_SharedA1A2) {
    // "ACGTACGT" aparece en A1 y A2, pero todas las ocurrencias en A2 son copias
    // de A1 (frases LZ77 secundarias). locate_extremal solo ve A1.
    // LCA real = A(1), nuestro resultado = A1(3) o similar descendiente de A.
    auto [pmin, pmax] = idx_.locate_extremal("ACGTACGT");
    ASSERT_NE(pmin, SIZE_MAX);
    const int our_lca = classifier_.classify_extremal(pmin, pmax);
    // Verificamos solo la propiedad débil: nuestro resultado ⊆ subárbol(A=1)
    EXPECT_TRUE(is_in_subtree(our_lca, 1))
        << "LIMITACIÓN: LCA=" << our_lca << " debería estar en subárbol de A(1)";
    // Documentamos que NO obtenemos el LCA exacto en este caso:
    EXPECT_NE(our_lca, 1) << "INESPERADO: si pasa, el índice ahora encuentra A2 como primaria";
}

TEST_F(ClassifyTest, LCA_KnownLimitation_SharedA2B1) {
    // "GGGGGGGG" aparece en A2 y B1. Si las ocurrencias en B1 son todas copias
    // de A2 (frases secundarias), locate_extremal solo ve A2.
    // LCA real = root(0), nuestro resultado = A2(4) o similar descendiente de root.
    auto [pmin, pmax] = idx_.locate_extremal("GGGGGGGG");
    ASSERT_NE(pmin, SIZE_MAX);
    const int our_lca = classifier_.classify_extremal(pmin, pmax);
    EXPECT_TRUE(is_in_subtree(our_lca, 0))
        << "LIMITACIÓN: LCA=" << our_lca << " debería estar en subárbol de root(0)";
    EXPECT_NE(our_lca, 0) << "INESPERADO: si pasa, el índice ahora encuentra B1 como primaria";
}

// ─────────────────────────────────────────────────────────────────────────────
// classify_read() — pipeline completa MEMExtractor → locate_extremal → LCA
// ─────────────────────────────────────────────────────────────────────────────

TEST_F(ClassifyTest, ClassifyRead_UniqueA1) {
    // Patrón único de A1 → classify_read debe retornar nodo en subárbol de A1(3)
    const int result = classifier_.classify_read("AAAAACGT", idx_, *ext_, 8);
    ASSERT_NE(result, -1);
    EXPECT_EQ(tree_.lca(result, 3), 3)
        << "classify_read retornó " << result << ", no en subárbol de A1(3)";
}

TEST_F(ClassifyTest, ClassifyRead_UniqueB2) {
    const int result = classifier_.classify_read("CCCCCCCC", idx_, *ext_, 8);
    ASSERT_NE(result, -1);
    EXPECT_EQ(tree_.lca(result, 6), 6)
        << "classify_read retornó " << result << ", no en subárbol de B2(6)";
}

TEST_F(ClassifyTest, ClassifyRead_SharedB1B2) {
    // "TGCATGCA" aparece en B1 y B2 → LCA debe estar en subárbol de B(2)
    const int result = classifier_.classify_read("TGCATGCA", idx_, *ext_, 8);
    ASSERT_NE(result, -1);
    EXPECT_EQ(tree_.lca(result, 2), 2)
        << "classify_read retornó " << result << ", no en subárbol de B(2)";
}

TEST_F(ClassifyTest, ClassifyRead_NoMatch) {
    EXPECT_EQ(classifier_.classify_read("NNNNNNNN", idx_, *ext_, 8), -1);
}

TEST_F(ClassifyTest, ClassifyRead_InSubtreeOfGroundTruth) {
    // Propiedad: classify_read(P) ∈ subárbol(bf_classify(P))
    const std::vector<std::pair<std::string, std::string>> cases = {
        {"TGCATGCA",           "TGCATGCA"},   // B1+B2 → bf=B(2)
        {"CCCCCCCC",           "CCCCCCCC"},   // B2    → bf=B2(6)
        {"AAAAACGT",           "AAAAACGT"},   // A1    → bf=A1(3)
    };
    for (const auto& [read, pattern] : cases) {
        SCOPED_TRACE("read=" + read);
        const int our = classifier_.classify_read(read, idx_, *ext_, 8);
        const int bf  = bf_classify(pattern);
        if (our == -1 || bf == -1) continue;
        EXPECT_EQ(tree_.lca(our, bf), bf)
            << "our=" << our << " no está en subárbol de bf=" << bf;
    }
}
