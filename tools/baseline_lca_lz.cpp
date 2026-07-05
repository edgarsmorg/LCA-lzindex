/**
 * baseline_lca_lz — benchmark del pipeline locate_extremal usando LZ77-index.
 *
 * Análogo a baseline_lca_sr pero para el LZ77-index. Carga un índice pre-construido
 * (build_lz_index) y mide locate_extremal sobre un conjunto de patrones.
 *
 * Uso:
 *   baseline_lca_lz <text_file> <lz_prefix> <patterns_file> [--mode=locate|info] [--name=<dataset>]
 *
 *   text_file     Texto de referencia en texto plano (para calcular n y bpc)
 *   lz_prefix     Prefijo del índice pre-construido (e.g. indexes/migenoma)
 *                 Espera: <prefix>.meta  .grid  .rcsa_fwd  .rcsa_rev
 *   patterns_file Archivo con patrones, uno por línea
 *   --mode=locate (default) Mide locate_extremal: µs/query, hits, peak RSS
 *   --mode=info   Solo imprime tamaño del índice y sale
 *   --name=<s>    Nombre del dataset para el output (default: basename de text_file)
 *
 * Salida (stdout):
 *   - z (frases LZ77), tamaño del índice (bytes + bpc)
 *   - Patrones con/sin ocurrencias primarias
 *   - µs/query de locate_extremal
 *   - Peak RSS en MB
 */

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <sys/resource.h>

#include "bench_common.hpp"
#include "index.hpp"

using namespace lz77tax;
using Clock = std::chrono::steady_clock;
using us    = std::chrono::microseconds;
namespace fs = std::filesystem;

static void usage(const char* argv0) {
    std::cerr << "Uso: " << argv0
              << " <text_file> <lz_prefix> <patterns_file>"
              << " [--mode=locate|info] [--name=<dataset>]\n";
}

int main(int argc, char** argv) {
    if (argc < 4) { usage(argv[0]); return 1; }

    const fs::path text_path    = argv[1];
    const fs::path lz_prefix    = argv[2];
    const fs::path patterns_path = argv[3];

    std::string mode    = "locate";
    std::string dataset = bench::basename(text_path);

    for (int i = 4; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--mode=locate")    mode = "locate";
        else if (arg == "--mode=info") mode = "info";
        else if (auto v = bench::option_value(arg, "name"); !v.empty()) dataset = v;
        else { std::cerr << "Argumento desconocido: " << arg << "\n"; usage(argv[0]); return 1; }
    }

    // ── Leer texto (necesario para verificación Patricia) ─────────────────────
    std::ifstream tf(text_path, std::ios::binary);
    std::ostringstream tbuf;
    tbuf << tf.rdbuf();
    const std::string text = tbuf.str();

    // ── Cargar índice ─────────────────────────────────────────────────────────
    std::cout << "Cargando LZ77-index (" << lz_prefix.string() << ")..." << std::flush;
    auto t0 = Clock::now();
    LZ77Index idx;
    try {
        idx.load(lz_prefix, text);
    } catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << "\n";
        std::cerr << "¿Construiste el índice con build_lz_index?\n";
        return 2;
    }
    double load_s = std::chrono::duration<double>(Clock::now() - t0).count();
    std::cout << " " << std::fixed << std::setprecision(1) << load_s << "s\n";

    // ── Tamaño del índice ─────────────────────────────────────────────────────
    const std::size_t n_bytes    = fs::file_size(text_path);
    const std::size_t grid_bytes = bench::lz_grid_bytes(idx);
    const std::size_t trie_bytes = idx.trie_bytes();
    const std::size_t total      = grid_bytes + trie_bytes;
    const double bpc = n_bytes ? total * 8.0 / static_cast<double>(n_bytes) : 0.0;

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "frases LZ77 (z)   : " << idx.phrase_count() << "\n";
    std::cout << "tamaño del índice : " << total << " bytes ("
              << total / (1024.0 * 1024.0) << " MB, " << bpc << " bpc)\n";
    std::cout << "  grid            : " << grid_bytes << " bytes\n";
    std::cout << "  tries SST+rev   : " << trie_bytes << " bytes\n";

    if (mode == "info") return 0;

    // ── Leer patrones ─────────────────────────────────────────────────────────
    std::vector<std::string> patterns;
    try { patterns = bench::read_patterns(patterns_path); }
    catch (const std::exception& e) { std::cerr << e.what() << "\n"; return 3; }
    std::cout << "patrones  : " << patterns.size() << "\n\n";

    // ── Warmup ────────────────────────────────────────────────────────────────
    for (std::size_t i = 0; i < std::min(patterns.size(), std::size_t{10}); ++i)
        (void)idx.locate_extremal(patterns[i]);

    // ── Benchmark locate_extremal ─────────────────────────────────────────────
    std::size_t hits    = 0;
    std::size_t no_hits = 0;
    std::size_t d       = 0;  // sink para evitar dead-code elimination

    auto tb0 = Clock::now();
    for (const auto& p : patterns) {
        const auto [mn, mx] = idx.locate_extremal(p);
        if (mn == SIZE_MAX) { ++no_hits; }
        else                { ++hits; d += mn + mx; }
    }
    double total_us = std::chrono::duration_cast<us>(Clock::now() - tb0).count();
    (void)d;

    const double us_per_query = total_us / static_cast<double>(patterns.size());

    struct rusage ru{};
    getrusage(RUSAGE_SELF, &ru);
    const long peak_kb = ru.ru_maxrss;

    std::cout << "=== locate_extremal (LZ77-index) ===\n";
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  patrones con ocurrencias : " << hits    << "\n";
    std::cout << "  patrones sin ocurrencias : " << no_hits << "\n";
    std::cout << "  tiempo total             : " << total_us / 1e6 << " s\n";
    std::cout << "  µs/query                 : " << us_per_query << "\n";
    std::cout << "  peak RSS                 : " << peak_kb << " KB ("
              << std::setprecision(1) << peak_kb / 1024.0 << " MB)\n";

    return 0;
}
