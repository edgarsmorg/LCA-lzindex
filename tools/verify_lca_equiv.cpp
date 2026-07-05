/**
 * verify_lca_equiv — verificación de correctitud del diseño primarias-only.
 *
 * Contrasta, para cada patrón, el LCA taxonómico obtenido por dos vías:
 *   - LZ77Index::locate_extremal → (pos_min, pos_max) de ocurrencias PRIMARIAS
 *     → LCA(genoma(pos_min), genoma(pos_max)).
 *   - FM-index de sdsl (ground truth): enumera TODAS las ocurrencias, mapea
 *     cada una a su genoma y calcula el LCA del conjunto completo.
 *
 * El argumento a validar (CLAUDE.md): como la referencia sigue orden DFS, el
 * LCA de todas las ocurrencias == LCA de sus extremos. El LZ77 solo ve primarias,
 * así que su LCA debe ser EQUAL al full, o un DESCENDIENTE (más específico) en el
 * caso patológico de MEMs internos a una frase. Nunca un ancestro ni incomparable.
 *
 * Uso:
 *   verify_lca_equiv <ref_text> <tree_tsv> <genomes_tsv> <patterns_file> [--limit=N]
 */
#include "bench_common.hpp"
#include "classify_io.hpp"
#include "taxonomy/lca.hpp"
#include "index.hpp"

#include <sdsl/suffix_arrays.hpp>

#include <algorithm>
#include <climits>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using lz77tax::LZ77Index;
using lz77tax::PhyloTree;
namespace fs = std::filesystem;

// Mapeo posición → nodo-hoja por búsqueda binaria sobre rangos DFS (O(log g)).
struct GenomeMap {
    std::vector<size_t> starts, ends;
    std::vector<int>    node_ids;
    void build(std::vector<lz77tax::Classifier::GenomeRange> r) {
        std::sort(r.begin(), r.end(),
                  [](const auto& a, const auto& b) { return a.start < b.start; });
        for (const auto& g : r) {
            starts.push_back(g.start); ends.push_back(g.end); node_ids.push_back(g.node_id);
        }
    }
    int genome_of(size_t pos) const {
        auto it = std::upper_bound(starts.begin(), starts.end(), pos);
        if (it == starts.begin()) return -1;
        const size_t i = static_cast<size_t>(it - starts.begin() - 1);
        return (pos < ends[i]) ? node_ids[i] : -1;
    }
};

int main(int argc, char** argv) {
    if (argc < 5) {
        std::cerr << "usage: " << argv[0]
                  << " <ref_text> <tree_tsv> <genomes_tsv> <patterns_file> [--limit=N]\n";
        return 1;
    }
    const fs::path ref_path = argv[1];
    const std::string tree_path = argv[2];
    const std::string genomes_path = argv[3];
    const fs::path patterns_path = argv[4];
    size_t limit = SIZE_MAX;
    std::string label;
    fs::path csv_path;
    for (int i = 5; i < argc; ++i) {
        const std::string arg = argv[i];
        if (auto v = bench::option_value(arg, "limit"); !v.empty()) limit = std::stoull(v);
        else if (auto v = bench::option_value(arg, "label"); !v.empty()) label = v;
        else if (auto v = bench::option_value(arg, "csv"); !v.empty()) csv_path = v;
    }

    // Texto de referencia.
    std::ifstream tf(ref_path, std::ios::binary);
    std::ostringstream tbuf; tbuf << tf.rdbuf();
    const std::string text = tbuf.str();

    // LZ77Index (in-memory) + FM-index ground truth.
    std::cout << "Construyendo LZ77Index..." << std::flush;
    LZ77Index lz; lz.build(text);
    std::cout << " ok\nConstruyendo FM-index (ground truth)..." << std::flush;
    sdsl::csa_wt<> fm;
    sdsl::construct_im(fm, text, 1);
    std::cout << " ok\n";

    // Árbol + rangos de genoma.
    PhyloTree tree;
    std::vector<lz77tax::Classifier::GenomeRange> ranges;
    if (!classify_io::load_tree(tree_path, tree)) return 4;
    if (!classify_io::load_genomes(genomes_path, ranges)) return 5;
    GenomeMap gmap; gmap.build(std::move(ranges));

    const auto patterns = bench::read_patterns(patterns_path);

    enum Cat { EQUAL, LZ_DESCENDANT, LZ_ANCESTOR, INCOMPARABLE, LZ_MISS, NO_OCC, N_CAT };
    const char* names[N_CAT] =
        {"EQUAL", "LZ_DESCENDANT", "LZ_ANCESTOR(BUG)", "INCOMPARABLE(BUG)",
         "LZ_MISS", "NO_OCC"};
    std::map<int, int> counts;
    int examples_bug = 0;

    size_t done = 0;
    for (const auto& p : patterns) {
        if (done >= limit) break;
        ++done;

        // Ground truth: todas las ocurrencias → LCA del conjunto (== LCA de extremos).
        auto occ = sdsl::locate(fm, p.begin(), p.end());
        if (occ.empty()) { counts[NO_OCC]++; continue; }
        int lca_full = -1;
        for (auto pos : occ) {
            const int g = gmap.genome_of(static_cast<size_t>(pos));
            if (g < 0) continue;
            lca_full = (lca_full < 0) ? g : tree.lca(lca_full, g);
        }

        // LZ77 primarias: extremos → LCA.
        const auto [mn, mx] = lz.locate_extremal(p);
        if (mn == SIZE_MAX) { counts[LZ_MISS]++; continue; }
        const int gmin = gmap.genome_of(mn), gmax = gmap.genome_of(mx);
        const int lca_lz = (gmin < 0 || gmax < 0) ? -1 : tree.lca(gmin, gmax);

        int cat;
        if (lca_lz == lca_full)                  cat = EQUAL;
        else if (tree.lca(lca_lz, lca_full) == lca_full) cat = LZ_DESCENDANT;   // lz más específico
        else if (tree.lca(lca_lz, lca_full) == lca_lz)   cat = LZ_ANCESTOR;     // BUG
        else                                             cat = INCOMPARABLE;    // BUG
        counts[cat]++;

        if ((cat == LZ_ANCESTOR || cat == INCOMPARABLE) && examples_bug < 5) {
            ++examples_bug;
            std::cerr << "  [" << names[cat] << "] patrón='" << p
                      << "' lz_node=" << lca_lz << " full_node=" << lca_full
                      << " occ=" << occ.size() << "\n";
        }
    }

    std::cout << "\n=== Verificación LCA: LZ77 (primarias) vs FM-index (todas) ===\n";
    std::cout << "Patrones evaluados: " << done << "\n\n";
    for (int c = 0; c < N_CAT; ++c) {
        const int n = counts[c];
        const double pct = done > 0 ? 100.0 * n / static_cast<double>(done) : 0.0;
        std::cout << "  " << names[c] << std::string(std::max(1, 20 - (int)std::string(names[c]).size()), ' ')
                  << n << "  (" << pct << "%)\n";
    }
    const int bugs = counts[LZ_ANCESTOR] + counts[INCOMPARABLE];
    std::cout << "\n" << (bugs == 0
        ? "OK: ningún caso ancestro/incomparable — el invariante primarias-only se cumple."
        : "FALLO: hay casos ancestro/incomparable (bug de correctitud).") << "\n";

    if (!csv_path.empty()) {
        const std::string header =
            "label,evaluated,equal,lz_descendant,lz_ancestor,incomparable,lz_miss,no_occ";
        std::string row = bench::csv_quote(label.empty() ? patterns_path.string() : label);
        row += "," + std::to_string(done);
        for (int c = 0; c < N_CAT; ++c) row += "," + std::to_string(counts[c]);
        bench::append_csv_row(csv_path, header, row);
    }
    return bugs == 0 ? 0 : 1;
}
