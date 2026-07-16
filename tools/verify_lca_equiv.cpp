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
 * LCA de todas las ocurrencias == LCA de sus extremos. El LZ77 solo ve primarias
 * (todo MEM que ocurre tiene >=1 primaria: su ocurrencia más a la izquierda), así que
 * su LCA debe ser EQUAL al full, o un DESCENDIENTE (más específico) cuando el extremo
 * derecho lo alcanza una secundaria interior a una frase. Nunca un ancestro ni incomparable.
 *
 * Construye LZ77Index + FM-index UNA sola vez y puede barrer varios conjuntos de
 * patrones (p. ej. distintos largos) en el mismo proceso.
 *
 * Uso:
 *   verify_lca_equiv <ref_text> <tree_tsv> <genomes_tsv> [<patterns_file>]
 *       [--label=<L>] [--limit=N] [--csv=<path>] [--pat=<label>:<file> ...]
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
#include <utility>
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

enum Cat { EQUAL, LZ_DESCENDANT, LZ_ANCESTOR, INCOMPARABLE, LZ_MISS, NO_OCC, N_CAT };
static const char* kCatNames[N_CAT] =
    {"EQUAL", "LZ_DESCENDANT", "LZ_ANCESTOR(BUG)", "INCOMPARABLE(BUG)", "LZ_MISS", "NO_OCC"};

// Evalúa un conjunto de patrones; imprime resumen y (opcional) agrega fila CSV.
// Devuelve el número de casos-bug (ancestro + incomparable).
template <class FM>
static int run_set(const std::string& label, const fs::path& patterns_path, size_t limit,
                   const LZ77Index& lz, const FM& fm, const PhyloTree& tree,
                   const GenomeMap& gmap, const fs::path& csv_path) {
    const auto patterns = bench::read_patterns(patterns_path);
    std::map<int, int> counts;
    int examples_bug = 0;
    size_t done = 0;

    for (const auto& p : patterns) {
        if (done >= limit) break;
        ++done;

        auto occ = sdsl::locate(fm, p.begin(), p.end());
        if (occ.empty()) { counts[NO_OCC]++; continue; }
        int lca_full = -1;
        for (auto pos : occ) {
            const int g = gmap.genome_of(static_cast<size_t>(pos));
            if (g < 0) continue;
            lca_full = (lca_full < 0) ? g : tree.lca(lca_full, g);
        }

        const auto [mn, mx] = lz.locate_extremal(p);
        if (mn == SIZE_MAX) { counts[LZ_MISS]++; continue; }
        const int gmin = gmap.genome_of(mn), gmax = gmap.genome_of(mx);
        const int lca_lz = (gmin < 0 || gmax < 0) ? -1 : tree.lca(gmin, gmax);

        int cat;
        if (lca_lz == lca_full)                          cat = EQUAL;
        else if (tree.lca(lca_lz, lca_full) == lca_full) cat = LZ_DESCENDANT;  // lz más específico
        else if (tree.lca(lca_lz, lca_full) == lca_lz)   cat = LZ_ANCESTOR;    // BUG
        else                                             cat = INCOMPARABLE;   // BUG
        counts[cat]++;

        if ((cat == LZ_ANCESTOR || cat == INCOMPARABLE) && examples_bug < 5) {
            ++examples_bug;
            std::cerr << "  [" << kCatNames[cat] << "] patrón='" << p
                      << "' lz_node=" << lca_lz << " full_node=" << lca_full
                      << " occ=" << occ.size() << "\n";
        }
    }

    std::cout << "── [" << label << "]  patrones evaluados: " << done << "\n";
    for (int c = 0; c < N_CAT; ++c) {
        const int n = counts[c];
        if (n == 0 && c != EQUAL && c != LZ_DESCENDANT) continue;
        const double pct = done > 0 ? 100.0 * n / static_cast<double>(done) : 0.0;
        std::cout << "     " << kCatNames[c]
                  << std::string(std::max(1, 20 - (int)std::string(kCatNames[c]).size()), ' ')
                  << n << "  (" << pct << "%)\n";
    }

    const int bugs = counts[LZ_ANCESTOR] + counts[INCOMPARABLE];
    if (!csv_path.empty()) {
        const std::string header =
            "label,evaluated,equal,lz_descendant,lz_ancestor,incomparable,lz_miss,no_occ";
        std::string row = bench::csv_quote(label) + "," + std::to_string(done);
        for (int c = 0; c < N_CAT; ++c) row += "," + std::to_string(counts[c]);
        bench::append_csv_row(csv_path, header, row);
    }
    return bugs;
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "usage: " << argv[0]
                  << " <ref_text> <tree_tsv> <genomes_tsv> [<patterns_file>]"
                  << " [--label=<L>] [--limit=N] [--csv=<path>] [--pat=<label>:<file> ...]\n";
        return 1;
    }
    const fs::path ref_path = argv[1];
    const std::string tree_path = argv[2];
    const std::string genomes_path = argv[3];

    size_t limit = SIZE_MAX;
    std::string label = "patterns";
    fs::path csv_path;
    fs::path positional_pat;
    std::vector<std::pair<std::string, fs::path>> pat_specs;  // (label, file)

    for (int i = 4; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg.rfind("--", 0) != 0) { positional_pat = arg; continue; }
        if (auto v = bench::option_value(arg, "limit"); !v.empty()) limit = std::stoull(v);
        else if (auto v = bench::option_value(arg, "label"); !v.empty()) label = v;
        else if (auto v = bench::option_value(arg, "csv"); !v.empty()) csv_path = v;
        else if (auto v = bench::option_value(arg, "pat"); !v.empty()) {
            const auto colon = v.find(':');
            if (colon == std::string::npos) {
                std::cerr << "formato --pat esperado label:archivo, recibido: " << v << "\n";
                return 1;
            }
            pat_specs.emplace_back(v.substr(0, colon), fs::path(v.substr(colon + 1)));
        }
    }
    if (pat_specs.empty()) {
        if (positional_pat.empty()) {
            std::cerr << "ERROR: falta patterns_file o al menos un --pat=label:file\n";
            return 1;
        }
        pat_specs.emplace_back(label, positional_pat);
    }

    // Texto de referencia.
    std::ifstream tf(ref_path, std::ios::binary);
    std::ostringstream tbuf; tbuf << tf.rdbuf();
    const std::string text = tbuf.str();

    std::cout << "Construyendo LZ77Index..." << std::flush;
    LZ77Index lz; lz.build(text);
    std::cout << " ok\nConstruyendo FM-index (ground truth)..." << std::flush;
    sdsl::csa_wt<> fm;
    sdsl::construct_im(fm, text, 1);
    std::cout << " ok\n";

    PhyloTree tree;
    std::vector<lz77tax::Classifier::GenomeRange> ranges;
    if (!classify_io::load_tree(tree_path, tree)) return 4;
    if (!classify_io::load_genomes(genomes_path, ranges)) return 5;
    GenomeMap gmap; gmap.build(std::move(ranges));

    std::cout << "\n=== Verificación LCA: LZ77 (primarias) vs FM-index (todas) ===\n";
    int total_bugs = 0;
    for (const auto& [lab, file] : pat_specs)
        total_bugs += run_set(lab, file, limit, lz, fm, tree, gmap, csv_path);

    std::cout << "\n" << (total_bugs == 0
        ? "OK: ningún caso ancestro/incomparable — el invariante primarias-only se cumple."
        : "FALLO: hay casos ancestro/incomparable (bug de correctitud).") << "\n";
    return total_bugs == 0 ? 0 : 1;
}
