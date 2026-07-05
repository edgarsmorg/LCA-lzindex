/**
 * bench_classify_compare — benchmark de la pipeline de clasificación CON MEMs:
 * read → extracción de MEMs → locate_extremal por MEM → LCA, comparando
 * LZ77Index vs SrIndexLocator.
 *
 * A diferencia de bench_locate_compare (locate puro de patrones de largo fijo),
 * aquí cada query es una READ: se extraen sus MEMs (min_len) y cada MEM se
 * localiza para computar el LCA taxonómico. Se reporta throughput por índice,
 * desglose por fase (MEM / locate / LCA), precisión vs ground-truth y stats de MEMs.
 *
 * La extracción de MEMs es independiente del índice (usa un FM-index sdsl),
 * por eso se cronometra una sola vez y se suma a ambos para el total por read.
 *
 * Uso:
 *   bench_classify_compare <ref_text> <tree_tsv> <genomes_tsv> <reads_tsv>
 *                          <lz_prefix> <sr_dir>
 *                          [--s=16] [--name=<data>] [--reps=3] [--min-mem=31]
 *                          [--csv=<path>]
 */

#include "bench_common.hpp"
#include "classify_io.hpp"
#include "index.hpp"
#include "taxonomy/lca.hpp"
#include "taxonomy/classifier.hpp"
#include "mem/extractor.hpp"
#include "baseline/sr_index_locator.hpp"

#include <sdsl/suffix_arrays.hpp>

#include <climits>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

using namespace lz77tax;
using classify_io::ReadEntry;
namespace fs = std::filesystem;

static void usage(const char* argv0) {
    std::cerr << "usage: " << argv0
              << " <ref_text> <tree_tsv> <genomes_tsv> <reads_tsv> <lz_prefix> <sr_dir>"
                 " [--s=16] [--name=<data>] [--reps=3] [--min-mem=31] [--csv=<path>]\n";
}

// Precisión de un nodo clasificado respecto al verdadero.
enum class Acc { CORRECT, OVER, INCORRECT, UNCLASSIFIED };

static Acc accuracy(int result, int truth, const PhyloTree& tree) {
    if (result < 0) return Acc::UNCLASSIFIED;
    if (result == truth) return Acc::CORRECT;
    // result ancestro de truth ⇒ clasificación menos específica (esperado en primary-only)
    if (tree.lca(result, truth) == result) return Acc::OVER;
    return Acc::INCORRECT;  // incomparable o descendiente inesperado
}

struct PhaseTimes {
    double locate_s = 0.0;
    double lca_s = 0.0;
};

// Clasifica todas las reads con un índice, reusando los MEMs precomputados.
// Cronometra locate y LCA por separado; devuelve el nodo clasificado por read.
template <class Index>
static PhaseTimes classify_all(const std::vector<ReadEntry>& reads,
                               const std::vector<MEMList>& mems_per_read,
                               const Index& index,
                               const Classifier& classifier,
                               const PhyloTree& tree,
                               std::vector<int>& results_out,
                               std::uint64_t& checksum) {
    PhaseTimes pt;
    results_out.assign(reads.size(), -1);

    for (std::size_t r = 0; r < reads.size(); ++r) {
        const std::string& seq = reads[r].sequence;
        int result = -1;
        for (const auto& mem : mems_per_read[r]) {
            const std::string pattern = seq.substr(mem.query_start, mem.length);

            const auto tl = bench::Clock::now();
            const auto [mn, mx] = index.locate_extremal(pattern);
            pt.locate_s += std::chrono::duration<double>(bench::Clock::now() - tl).count();

            if (mn == SIZE_MAX) continue;

            const auto tc = bench::Clock::now();
            const int node = classifier.classify_extremal(mn, mx);
            if (node >= 0) result = (result < 0) ? node : tree.lca(result, node);
            pt.lca_s += std::chrono::duration<double>(bench::Clock::now() - tc).count();
        }
        results_out[r] = result;
        checksum ^= static_cast<std::uint64_t>(result + 1) * 0x9e3779b97f4a7c15ULL
                  + (checksum << 6) + (checksum >> 2);
    }
    return pt;
}

static std::string make_row(const std::string& dataset, std::uintmax_t n_bytes,
                            const std::string& index_name, std::size_t s,
                            std::size_t rep, const std::vector<ReadEntry>& reads,
                            std::size_t read_len, std::size_t min_mem,
                            double mem_s, const PhaseTimes& pt,
                            double avg_mems, double avg_mem_len, double pct_no_mem,
                            const std::vector<int>& results, const PhyloTree& tree,
                            std::uint64_t checksum) {
    const std::size_t n = reads.size();
    const double total_s = mem_s + pt.locate_s + pt.lca_s;
    const double us_per_read = total_s * 1e6 / static_cast<double>(n);
    const double reads_per_s = total_s > 0 ? static_cast<double>(n) / total_s : 0.0;

    std::size_t correct = 0, over = 0, incorrect = 0, unclassified = 0;
    for (std::size_t r = 0; r < n; ++r) {
        switch (accuracy(results[r], reads[r].true_node_id, tree)) {
            case Acc::CORRECT:      ++correct; break;
            case Acc::OVER:         ++over; break;
            case Acc::INCORRECT:    ++incorrect; break;
            case Acc::UNCLASSIFIED: ++unclassified; break;
        }
    }
    auto pct = [&](std::size_t c) { return 100.0 * static_cast<double>(c) / static_cast<double>(n); };

    return bench::csv_quote(dataset) + "," + std::to_string(n_bytes) + "," + index_name + "," +
        std::to_string(s) + "," + std::to_string(rep) + "," + std::to_string(n) + "," +
        std::to_string(read_len) + "," + std::to_string(min_mem) + "," +
        bench::seconds(reads_per_s) + "," + bench::seconds(us_per_read) + "," +
        bench::seconds(mem_s * 1e6 / n) + "," + bench::seconds(pt.locate_s * 1e6 / n) + "," +
        bench::seconds(pt.lca_s * 1e6 / n) + "," + bench::seconds(avg_mems) + "," +
        bench::seconds(avg_mem_len) + "," + bench::seconds(pct_no_mem) + "," +
        bench::seconds(pct(correct)) + "," + bench::seconds(pct(over)) + "," +
        bench::seconds(pct(incorrect)) + "," + bench::seconds(pct(unclassified)) + "," +
        std::to_string(checksum);
}

int main(int argc, char** argv) {
    if (argc < 7) { usage(argv[0]); return 1; }

    const fs::path ref_path = argv[1];
    const std::string tree_path = argv[2];
    const std::string genomes_path = argv[3];
    const std::string reads_path = argv[4];
    const fs::path lz_prefix = argv[5];
    const fs::path sr_dir = argv[6];

    std::size_t s = 16, reps = 3, min_mem = 31;
    std::string dataset = bench::basename(ref_path);
    fs::path csv_path;

    for (int i = 7; i < argc; ++i) {
        const std::string arg = argv[i];
        if (auto v = bench::option_value(arg, "s"); !v.empty()) s = std::stoull(v);
        else if (auto v = bench::option_value(arg, "name"); !v.empty()) dataset = v;
        else if (auto v = bench::option_value(arg, "reps"); !v.empty()) reps = std::stoull(v);
        else if (auto v = bench::option_value(arg, "min-mem"); !v.empty()) min_mem = std::stoull(v);
        else if (auto v = bench::option_value(arg, "csv"); !v.empty()) csv_path = v;
        else { std::cerr << "unknown argument: " << arg << "\n"; usage(argv[0]); return 1; }
    }

    // ── Carga (no cronometrada) ──────────────────────────────────────────────
    std::string text;
    try { text = bench::read_text_file(ref_path); }
    catch (const std::exception& e) { std::cerr << e.what() << "\n"; return 2; }
    const std::uintmax_t n_bytes = text.size();

    LZ77Index lz;
    lz.load(lz_prefix, text);

    SrIndexLocator sr;
    try { sr.load(dataset, sr_dir.string(), s); }
    catch (const std::exception& e) {
        std::cerr << "Error al cargar sr-index: " << e.what() << "\n"; return 3;
    }

    sdsl::csa_wt<> fm;
    sdsl::construct_im(fm, text, 1);
    MEMExtractor ext(fm);

    PhyloTree tree;
    if (!classify_io::load_tree(tree_path, tree)) return 4;
    std::vector<Classifier::GenomeRange> genome_ranges;
    if (!classify_io::load_genomes(genomes_path, genome_ranges)) return 5;
    Classifier classifier;
    classifier.setup(tree, genome_ranges);

    std::vector<ReadEntry> reads;
    if (!classify_io::load_reads(reads_path, reads)) return 6;
    if (reads.empty()) { std::cerr << "reads vacíos\n"; return 7; }
    const std::size_t read_len = reads.front().sequence.size();

    // warmup
    for (std::size_t i = 0; i < std::min<std::size_t>(reads.size(), 10); ++i) {
        const auto m = ext.extract(reads[i].sequence, min_mem);
        for (const auto& mem : m) {
            const std::string p = reads[i].sequence.substr(mem.query_start, mem.length);
            (void)lz.locate_extremal(p);
            (void)sr.locate_extremal(p);
        }
    }

    const std::string header =
        "dataset,n_bytes,index,s,rep,n_reads,read_len,min_mem,reads_per_s,us_per_read,"
        "mem_us_per_read,locate_us_per_read,lca_us_per_read,avg_mems_per_read,avg_mem_len,"
        "pct_no_mem,pct_correct,pct_over,pct_incorrect,pct_unclassified,checksum";

    for (std::size_t rep = 0; rep < reps; ++rep) {
        // ── Fase MEM (compartida) ────────────────────────────────────────────
        std::vector<MEMList> mems_per_read(reads.size());
        const auto tm = bench::Clock::now();
        for (std::size_t r = 0; r < reads.size(); ++r)
            mems_per_read[r] = ext.extract(reads[r].sequence, min_mem);
        const double mem_s = std::chrono::duration<double>(bench::Clock::now() - tm).count();

        // stats de MEMs
        std::size_t total_mems = 0, total_mem_len = 0, reads_no_mem = 0;
        for (const auto& m : mems_per_read) {
            total_mems += m.size();
            if (m.empty()) ++reads_no_mem;
            for (const auto& mem : m) total_mem_len += mem.length;
        }
        const double avg_mems = static_cast<double>(total_mems) / reads.size();
        const double avg_mem_len = total_mems ? static_cast<double>(total_mem_len) / total_mems : 0.0;
        const double pct_no_mem = 100.0 * static_cast<double>(reads_no_mem) / reads.size();

        // ── Fase locate+LCA por índice ───────────────────────────────────────
        std::vector<int> res_lz, res_sr;
        std::uint64_t cs_lz = 1469598103934665603ULL, cs_sr = 1469598103934665603ULL;
        const PhaseTimes pt_lz = classify_all(reads, mems_per_read, lz, classifier, tree, res_lz, cs_lz);
        const PhaseTimes pt_sr = classify_all(reads, mems_per_read, sr, classifier, tree, res_sr, cs_sr);

        const std::string row_lz = make_row(dataset, n_bytes, "lz77", 0, rep, reads, read_len,
            min_mem, mem_s, pt_lz, avg_mems, avg_mem_len, pct_no_mem, res_lz, tree, cs_lz);
        const std::string row_sr = make_row(dataset, n_bytes, "sr", s, rep, reads, read_len,
            min_mem, mem_s, pt_sr, avg_mems, avg_mem_len, pct_no_mem, res_sr, tree, cs_sr);

        bench::append_csv_row(csv_path, header, row_lz);
        bench::append_csv_row(csv_path, header, row_sr);
        std::cout << row_lz << "\n" << row_sr << "\n";
    }
    return 0;
}
