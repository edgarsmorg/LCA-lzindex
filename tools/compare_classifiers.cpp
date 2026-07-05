/**
 * compare_classifiers — pipeline LCA end-to-end: LZ77Index vs SrIndexLocator.
 *
 * Uso:
 *   compare_classifiers <ref_text> <tree_tsv> <genomes_tsv> <reads_tsv>
 *                       <sr_data_dir> [sr=16] [--out=<csv_path>] [--min-mem=31]
 *
 *   ref_text      Texto de referencia (concatenación DFS en texto plano)
 *   tree_tsv      Árbol filogenético en TSV (generado por tree_json_to_tsv.py)
 *   genomes_tsv   Rangos de genomas en TSV (generado por tree_json_to_tsv.py)
 *   reads_tsv     Reads con ground-truth (generado por gen_classification_reads.py)
 *   sr_data_dir   Directorio con archivos .sdsl pre-construidos del sr-index
 *   sr            Subsampling rate del sr-index (default: 16)
 *   --out         Ruta del CSV de salida (default: results/compare/compare_sr<s>.csv)
 *   --min-mem     Longitud mínima de MEM (default: 31)
 *
 * Categorías de equivalencia emitidas en el CSV:
 *   EQUAL              — ambos retornan el mismo nodo LCA
 *   LZ_DESCENDANT_OF_SR — LZ77 en subárbol del sr-index (limitación primary-only esperada)
 *   SR_DESCENDANT_OF_LZ — inesperado, potencial bug
 *   INCOMPARABLE        — ramas distintas, bug
 *   LZ_UNCLASSIFIED     — LZ77 no encontró ocurrencias (read invisible para el índice)
 *   SR_UNCLASSIFIED     — sr-index no encontró ocurrencias
 *   BOTH_UNCLASSIFIED   — ninguno encontró ocurrencias
 */

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "bench_common.hpp"
#include "classify_io.hpp"
#include "index.hpp"
#include "taxonomy/lca.hpp"
#include "taxonomy/classifier.hpp"
#include "mem/extractor.hpp"
#include "baseline/sr_index_locator.hpp"

#include <sdsl/suffix_arrays.hpp>

using namespace lz77tax;
namespace fs = std::filesystem;
using Clock  = std::chrono::steady_clock;

// ── Parsing de argumentos ────────────────────────────────────────────────────

static void usage(const char* prog) {
    std::cerr
        << "Uso: " << prog
        << " <ref_text> <tree_tsv> <genomes_tsv> <reads_tsv> <sr_data_dir>"
           " [sr=16] [--out=<csv>] [--min-mem=31]\n";
}

using classify_io::ReadEntry;

// ── Categoría de equivalencia ─────────────────────────────────────────────────

enum class Category {
    EQUAL,
    LZ_DESCENDANT_OF_SR,
    SR_DESCENDANT_OF_LZ,
    INCOMPARABLE,
    LZ_UNCLASSIFIED,
    SR_UNCLASSIFIED,
    BOTH_UNCLASSIFIED,
};

static const char* cat_name(Category c) {
    switch (c) {
        case Category::EQUAL:                return "EQUAL";
        case Category::LZ_DESCENDANT_OF_SR:  return "LZ_DESCENDANT_OF_SR";
        case Category::SR_DESCENDANT_OF_LZ:  return "SR_DESCENDANT_OF_LZ";
        case Category::INCOMPARABLE:         return "INCOMPARABLE";
        case Category::LZ_UNCLASSIFIED:      return "LZ_UNCLASSIFIED";
        case Category::SR_UNCLASSIFIED:      return "SR_UNCLASSIFIED";
        case Category::BOTH_UNCLASSIFIED:    return "BOTH_UNCLASSIFIED";
    }
    return "UNKNOWN";
}

static Category categorize(int lz, int sr, const PhyloTree& tree) {
    if (lz == -1 && sr == -1) return Category::BOTH_UNCLASSIFIED;
    if (lz == -1)              return Category::LZ_UNCLASSIFIED;
    if (sr == -1)              return Category::SR_UNCLASSIFIED;
    if (lz == sr)              return Category::EQUAL;
    if (tree.lca(lz, sr) == sr) return Category::LZ_DESCENDANT_OF_SR;
    if (tree.lca(lz, sr) == lz) return Category::SR_DESCENDANT_OF_LZ;
    return Category::INCOMPARABLE;
}

// ─────────────────────────────────────────────────────────────────────────────

int main(int argc, char** argv) {
    if (argc < 6) { usage(argv[0]); return 1; }

    const std::string ref_path    = argv[1];
    const std::string tree_path   = argv[2];
    const std::string genomes_path = argv[3];
    const std::string reads_path  = argv[4];
    const std::string sr_data_dir = argv[5];

    std::size_t sr      = 16;
    std::string out_csv = "";
    std::size_t min_mem = 31;

    for (int i = 6; i < argc; ++i) {
        std::string arg(argv[i]);
        const std::string val_out    = bench::option_value(arg, "out");
        const std::string val_minmem = bench::option_value(arg, "min-mem");
        if (!val_out.empty())         out_csv = val_out;
        else if (!val_minmem.empty()) min_mem = std::stoul(val_minmem);
        else { try { sr = std::stoul(arg); } catch (...) {
            std::cerr << "Argumento desconocido: " << arg << "\n";
            usage(argv[0]); return 1;
        }}
    }

    // CSV de salida por defecto
    if (out_csv.empty()) {
        fs::create_directories("results/compare");
        out_csv = "results/compare/compare_sr" + std::to_string(sr) + ".csv";
    } else {
        fs::create_directories(fs::path(out_csv).parent_path());
    }

    // ── Cargar texto de referencia ────────────────────────────────────────────
    std::cout << "Cargando texto de referencia..." << std::flush;
    std::string ref_text;
    try { ref_text = bench::read_text_file(ref_path); }
    catch (const std::exception& e) { std::cerr << "\n" << e.what() << "\n"; return 2; }
    std::cout << " " << ref_text.size() << " bytes\n";

    // ── Construir LZ77Index ────────────────────────────────────────────────────
    std::cout << "Construyendo LZ77Index..." << std::flush;
    auto t0 = Clock::now();
    LZ77Index lz_idx;
    lz_idx.build(ref_text);
    double lz_build_s = std::chrono::duration<double>(Clock::now() - t0).count();
    std::cout << " " << std::fixed << std::setprecision(1) << lz_build_s << "s\n";

    // ── Construir FM-index + MEMExtractor ─────────────────────────────────────
    std::cout << "Construyendo FM-index para MEMExtractor..." << std::flush;
    t0 = Clock::now();
    sdsl::csa_wt<> fm;
    sdsl::construct_im(fm, ref_text, 1);
    MEMExtractor ext(fm);
    double fm_build_s = std::chrono::duration<double>(Clock::now() - t0).count();
    std::cout << " " << std::setprecision(1) << fm_build_s << "s\n";

    // ── Cargar sr-index ────────────────────────────────────────────────────────
    const std::string data_name = fs::path(ref_path).filename().string();
    std::cout << "Cargando sr-index (data_name=" << data_name
              << ", sr=" << sr << ")..." << std::flush;
    t0 = Clock::now();
    SrIndexLocator sr_loc;
    try {
        sr_loc.load(data_name, sr_data_dir, sr);
    } catch (const std::exception& e) {
        std::cerr << "\nError al cargar sr-index: " << e.what() << "\n";
        return 3;
    }
    double sr_load_s = std::chrono::duration<double>(Clock::now() - t0).count();
    std::cout << " " << std::setprecision(1) << sr_load_s << "s\n";

    // ── Cargar árbol filogenético ─────────────────────────────────────────────
    PhyloTree tree;
    if (!classify_io::load_tree(tree_path, tree)) return 4;
    std::cout << "Árbol cargado: " << tree.size() << " nodos\n";

    // ── Cargar rangos de genoma + configurar Classifier ───────────────────────
    std::vector<Classifier::GenomeRange> genome_ranges;
    if (!classify_io::load_genomes(genomes_path, genome_ranges)) return 5;

    Classifier classifier;
    classifier.setup(tree, genome_ranges);
    std::cout << "Genomas cargados: " << genome_ranges.size() << " rangos\n";

    // ── Cargar reads ──────────────────────────────────────────────────────────
    std::vector<ReadEntry> reads;
    if (!classify_io::load_reads(reads_path, reads)) return 6;
    std::cout << "Reads cargados: " << reads.size() << "\n\n";

    // ── Pipeline de clasificación y comparación ───────────────────────────────
    std::ofstream csv(out_csv);
    if (!csv) { std::cerr << "No se puede escribir: " << out_csv << "\n"; return 7; }
    csv << "read_id,true_genome_idx,true_node_id,lz_node,sr_node,category,"
           "lz_depth,sr_node_depth,true_depth\n";

    std::map<Category, int> counts;

    std::cout << "Clasificando " << reads.size() << " reads...\n";
    t0 = Clock::now();
    for (const auto& r : reads) {
        const int lz = classifier.classify_read(r.sequence, lz_idx, ext, min_mem);
        const int sr = classifier.classify_read(r.sequence, sr_loc, ext, min_mem);
        const Category cat = categorize(lz, sr, tree);
        counts[cat]++;

        const int lz_depth   = (lz >= 0) ? tree.depth(lz) : -1;
        const int sr_depth   = (sr >= 0) ? tree.depth(sr) : -1;
        const int true_depth = tree.depth(r.true_node_id);

        csv << r.read_id     << ','
            << r.true_genome_idx << ','
            << r.true_node_id << ','
            << lz            << ','
            << sr            << ','
            << cat_name(cat) << ','
            << lz_depth      << ','
            << sr_depth      << ','
            << true_depth    << '\n';
    }
    double classify_s = std::chrono::duration<double>(Clock::now() - t0).count();

    // ── Resumen ───────────────────────────────────────────────────────────────
    const int total = static_cast<int>(reads.size());
    const int classified_both = total
        - counts[Category::LZ_UNCLASSIFIED]
        - counts[Category::SR_UNCLASSIFIED]
        - counts[Category::BOTH_UNCLASSIFIED];

    std::cout << "\n=== Comparación LZ77-Index vs SR-Index (sr=" << sr << ") ===\n";
    std::cout << std::fixed << std::setprecision(1);
    std::cout << "Tiempo clasificación: " << classify_s << "s"
              << "  (" << classify_s / total * 1000 << " ms/read)\n\n";

    auto pct = [&](Category c) {
        return total > 0 ? 100.0 * counts[c] / total : 0.0;
    };

    std::cout << "  EQUAL              : " << counts[Category::EQUAL]
              << "  (" << pct(Category::EQUAL) << "%)\n";
    std::cout << "  LZ_DESCENDANT_OF_SR: " << counts[Category::LZ_DESCENDANT_OF_SR]
              << "  (" << pct(Category::LZ_DESCENDANT_OF_SR) << "%)  [esperado: primary-only]\n";
    std::cout << "  SR_DESCENDANT_OF_LZ: " << counts[Category::SR_DESCENDANT_OF_LZ]
              << "  (" << pct(Category::SR_DESCENDANT_OF_LZ) << "%)  [no esperado]\n";
    std::cout << "  INCOMPARABLE       : " << counts[Category::INCOMPARABLE]
              << "  (" << pct(Category::INCOMPARABLE) << "%)  [BUG si > 0]\n";
    std::cout << "  LZ_UNCLASSIFIED    : " << counts[Category::LZ_UNCLASSIFIED]
              << "  (" << pct(Category::LZ_UNCLASSIFIED) << "%)  [LZ77 sin primarias]\n";
    std::cout << "  SR_UNCLASSIFIED    : " << counts[Category::SR_UNCLASSIFIED]
              << "  (" << pct(Category::SR_UNCLASSIFIED) << "%)\n";
    std::cout << "  BOTH_UNCLASSIFIED  : " << counts[Category::BOTH_UNCLASSIFIED]
              << "  (" << pct(Category::BOTH_UNCLASSIFIED) << "%)\n";

    const double agreement = classified_both > 0
        ? 100.0 * (counts[Category::EQUAL] + counts[Category::LZ_DESCENDANT_OF_SR])
          / classified_both
        : 0.0;
    std::cout << "\n  Tasa de acuerdo (EQUAL+DESCENDANT) / clasificados_ambos: "
              << std::setprecision(2) << agreement << "%\n";

    if (counts[Category::INCOMPARABLE] > 0 || counts[Category::SR_DESCENDANT_OF_LZ] > 0) {
        std::cout << "\n  *** ATENCIÓN: hay "
                  << counts[Category::INCOMPARABLE] + counts[Category::SR_DESCENDANT_OF_LZ]
                  << " caso(s) inesperado(s). Revisar el CSV: " << out_csv << "\n";
    }

    std::cout << "\nCSV: " << out_csv << "\n";
    return 0;
}
