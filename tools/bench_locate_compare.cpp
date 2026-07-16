#include "bench_common.hpp"
#include "baseline/sr_index_locator.hpp"
#include "classify_io.hpp"
#include "taxonomy/lca.hpp"

#include <algorithm>
#include <climits>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

using lz77tax::LZ77Index;
using lz77tax::SrIndexLocator;
using lz77tax::PhyloTree;
namespace fs = std::filesystem;

static void usage(const char* argv0) {
    std::cerr << "usage: " << argv0
              << " <text_file> <patterns_file> <lz_prefix> <sr_dir>"
              << " [--s=16] [--name=<data_name>] [--reps=3] [--csv=<path>] [--oracle]"
              << " [--tree=<tree.tsv>] [--genomes=<genomes.tsv>]\n";
}

// Mapeo posición-de-texto → nodo-hoja del árbol, por búsqueda binaria.
// Los rangos de genoma no se solapan y siguen el orden DFS (start creciente),
// así que basta upper_bound sobre los starts. O(log g) en vez del O(g) del
// Classifier::genome_of naive (crítico: g llega a cientos de miles de genomas).
struct GenomeMap {
    std::vector<size_t> starts;
    std::vector<size_t> ends;
    std::vector<int>    node_ids;
    bool valid = false;

    void build(std::vector<lz77tax::Classifier::GenomeRange> ranges) {
        std::sort(ranges.begin(), ranges.end(),
                  [](const auto& a, const auto& b) { return a.start < b.start; });
        starts.reserve(ranges.size());
        ends.reserve(ranges.size());
        node_ids.reserve(ranges.size());
        for (const auto& g : ranges) {
            starts.push_back(g.start);
            ends.push_back(g.end);
            node_ids.push_back(g.node_id);
        }
        valid = !starts.empty();
    }

    int genome_of(size_t pos) const {
        auto it = std::upper_bound(starts.begin(), starts.end(), pos);
        if (it == starts.begin()) return -1;
        const size_t idx = static_cast<size_t>(it - starts.begin() - 1);
        return (pos < ends[idx]) ? node_ids[idx] : -1;  // -1 si cae en separador
    }
};

template <class Index>
static std::string run_one(const std::string& dataset, const std::string& index_name,
                           std::size_t s, std::size_t rep,
                           const std::vector<std::string>& patterns,
                           const Index& idx,
                           const PhyloTree* tree, const GenomeMap* gmap) {
    std::size_t hits = 0;
    std::size_t no_hits = 0;
    std::size_t checksum = 1469598103934665603ULL;
    const bool do_lca = (tree != nullptr && gmap != nullptr && gmap->valid);

    double locate_s = 0.0;
    double lca_s = 0.0;

    for (const auto& p : patterns) {
        const auto tl = bench::Clock::now();
        const auto [mn, mx] = idx.locate_extremal(p);
        locate_s += std::chrono::duration<double>(bench::Clock::now() - tl).count();

        if (mn == SIZE_MAX) {
            ++no_hits;
            checksum ^= p.size() + 0x9e3779b97f4a7c15ULL;
            continue;
        }
        ++hits;
        checksum ^= mn + 0x9e3779b97f4a7c15ULL + (checksum << 6) + (checksum >> 2);
        checksum ^= mx + 0x9e3779b97f4a7c15ULL + (checksum << 6) + (checksum >> 2);

        if (do_lca) {
            const auto tc = bench::Clock::now();
            const int nmin = gmap->genome_of(mn);
            const int nmax = gmap->genome_of(mx);
            const int node = (nmin < 0 || nmax < 0) ? -1 : tree->lca(nmin, nmax);
            lca_s += std::chrono::duration<double>(bench::Clock::now() - tc).count();
            checksum ^= static_cast<std::size_t>(node + 1) + 0x9e3779b97f4a7c15ULL
                        + (checksum << 6) + (checksum >> 2);
        }
    }

    const double np = static_cast<double>(patterns.size());
    const double locate_us = locate_s * 1e6 / np;
    const double lca_us    = do_lca ? lca_s * 1e6 / np : 0.0;
    const double total_us  = locate_us + lca_us;

    // us_per_query se mantiene = locate_us por retro-compatibilidad del merge.
    return bench::csv_quote(dataset) + "," + index_name + "," + std::to_string(s) + "," +
        std::to_string(rep) + "," + std::to_string(patterns.size()) + "," +
        std::to_string(hits) + "," + std::to_string(no_hits) + "," +
        bench::seconds(locate_s) + "," + bench::seconds(locate_us) + "," +
        bench::seconds(lca_us) + "," + bench::seconds(total_us) + "," +
        std::to_string(checksum);
}

int main(int argc, char** argv) {
    if (argc < 5) { usage(argv[0]); return 1; }

    const fs::path text_path = argv[1];
    const fs::path patterns_path = argv[2];
    const fs::path lz_prefix = argv[3];
    const fs::path sr_dir = argv[4];
    std::size_t s = 16;
    std::size_t reps = 3;
    bool use_oracle = false;   // --oracle: consultar el texto vía accesor LZ77
    std::string dataset = bench::basename(text_path);
    fs::path csv_path;
    fs::path tree_path;
    fs::path genomes_path;

    for (int i = 5; i < argc; ++i) {
        const std::string arg = argv[i];
        if (auto v = bench::option_value(arg, "s"); !v.empty()) s = std::stoull(v);
        else if (auto v = bench::option_value(arg, "name"); !v.empty()) dataset = v;
        else if (auto v = bench::option_value(arg, "reps"); !v.empty()) reps = std::stoull(v);
        else if (auto v = bench::option_value(arg, "csv"); !v.empty()) csv_path = v;
        else if (auto v = bench::option_value(arg, "tree"); !v.empty()) tree_path = v;
        else if (auto v = bench::option_value(arg, "genomes"); !v.empty()) genomes_path = v;
        else if (arg == "--oracle") use_oracle = true;
        else { std::cerr << "unknown argument: " << arg << "\n"; usage(argv[0]); return 1; }
    }

    // tree.tsv / genomes.tsv viven junto al reference.txt salvo override.
    if (tree_path.empty())    tree_path    = text_path.parent_path() / "tree.tsv";
    if (genomes_path.empty()) genomes_path = text_path.parent_path() / "genomes.tsv";

    const auto patterns = bench::read_patterns(patterns_path);

    std::ifstream tf(text_path, std::ios::binary);
    std::ostringstream tbuf;
    tbuf << tf.rdbuf();
    const std::string text = tbuf.str();

    LZ77Index lz;
    lz.load(lz_prefix, text);
    if (use_oracle) {
        // El texto deja de consultarse: las comparaciones se resuelven extrayendo
        // desde el parsing LZ77 (accesor), como haria un self-index.
        lz.build_oracle();
        lz.use_oracle(true);
        std::cerr << "accesor LZ77 activo: "
                  << lz.oracle_bytes() / 1048576.0 << " MB\n";
    }
    SrIndexLocator sr;
    sr.load(dataset, sr_dir.string(), s);

    // Árbol filogenético + mapeo genoma para el paso de LCA (opcional).
    PhyloTree tree;
    GenomeMap gmap;
    bool lca_ready = false;
    if (fs::exists(tree_path) && fs::exists(genomes_path)) {
        std::vector<lz77tax::Classifier::GenomeRange> ranges;
        if (classify_io::load_tree(tree_path.string(), tree) &&
            classify_io::load_genomes(genomes_path.string(), ranges)) {
            gmap.build(std::move(ranges));
            lca_ready = gmap.valid;
        }
    }
    if (!lca_ready) {
        std::cerr << "AVISO: sin tree.tsv/genomes.tsv (" << tree_path << ", "
                  << genomes_path << ") — se mide locate sin LCA\n";
    }
    const PhyloTree* tree_ptr = lca_ready ? &tree : nullptr;
    const GenomeMap* gmap_ptr = lca_ready ? &gmap : nullptr;

    for (std::size_t i = 0; i < std::min<std::size_t>(patterns.size(), 10); ++i) {
        (void)lz.locate_extremal(patterns[i]);
        (void)sr.locate_extremal(patterns[i]);
    }

    const std::string header =
        "dataset,index,s,rep,n_patterns,hits,no_hits,"
        "total_seconds,us_per_query,lca_us,total_us,checksum";
    for (std::size_t rep = 0; rep < reps; ++rep) {
        const auto lz_row = run_one(dataset, "lz77", 0, rep, patterns, lz, tree_ptr, gmap_ptr);
        const auto sr_row = run_one(dataset, "sr", s, rep, patterns, sr, tree_ptr, gmap_ptr);
        bench::append_csv_row(csv_path, header, lz_row);
        bench::append_csv_row(csv_path, header, sr_row);
        std::cout << lz_row << "\n" << sr_row << "\n";
    }
    return 0;
}
