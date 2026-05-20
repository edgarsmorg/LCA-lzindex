/**
 * baseline_lca_sr — benchmark del pipeline LCA usando sr-index como baseline.
 *
 * Uso:
 *   baseline_lca_sr <data_dir> <data_name> <patterns_file> [sr=16] [--mode=locate|info]
 *
 *   data_dir      Directorio con los archivos .sdsl pre-construidos
 *   data_name     Nombre base del dataset (e.g. "Escherichia_Coli")
 *   patterns_file Archivo con patrones, uno por línea
 *   sr            Subsampling rate usado al construir el índice (default: 16)
 *   --mode=locate (default) Mide locate_extremal: tiempo µs/query y estadísticas
 *   --mode=info   Solo imprime el tamaño del índice y sale
 *
 * Salida (stdout):
 *   - Tamaño del índice en bytes y bpc
 *   - Número de patrones con ocurrencias
 *   - Tiempo promedio de locate_extremal en µs/query
 *
 * Propósito: comparar con measure_index --patterns=<file> del LZ77-index.
 */

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "baseline/sr_index_locator.hpp"

using namespace lz77tax;
using Clock = std::chrono::steady_clock;
using us    = std::chrono::microseconds;
namespace fs = std::filesystem;

static void usage(const char* argv0) {
    std::cerr << "Uso: " << argv0
              << " <data_dir> <data_name> <patterns_file> [sr=16] [--mode=locate|info]\n";
}

int main(int argc, char** argv) {
    if (argc < 4) { usage(argv[0]); return 1; }

    const std::string data_dir      = argv[1];
    const std::string data_name     = argv[2];
    const std::string patterns_path = argv[3];

    std::size_t sr   = 16;
    std::string mode = "locate";

    for (int i = 4; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg == "--mode=locate")      mode = "locate";
        else if (arg == "--mode=info")   mode = "info";
        else {
            try { sr = std::stoul(arg); }
            catch (...) {
                std::cerr << "Argumento desconocido: " << arg << "\n";
                usage(argv[0]); return 1;
            }
        }
    }

    // ── Cargar índice ─────────────────────────────────────────────────────────
    std::cout << "Cargando sr-index (data_name=" << data_name
              << ", sr=" << sr << ")..." << std::flush;
    auto t0 = Clock::now();
    SrIndexLocator sr_loc;
    try {
        sr_loc.load(data_name, data_dir, sr);
    } catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << "\n";
        std::cerr << "¿Construiste el índice con sr=" << sr
                  << " en " << data_dir << "?\n";
        return 2;
    }
    double load_s = std::chrono::duration<double>(Clock::now() - t0).count();
    std::cout << " " << std::fixed << std::setprecision(1) << load_s << "s\n";

    // ── Tamaño del índice ─────────────────────────────────────────────────────
    // Estimamos sumando tamaños de los archivos sdsl relevantes en data_dir.
    // (serialize() reconstruye todo el índice en un stream, lo cual es lento;
    //  preferimos la suma de archivos para el benchmark.)
    std::size_t idx_bytes = 0;
    const std::string prefix = std::to_string(sr) + "_";
    for (const auto& entry : fs::directory_iterator(data_dir)) {
        const std::string fname = entry.path().filename().string();
        // Incluir solo los archivos del sr con el data_name correcto
        if (fname.find(data_name) != std::string::npos &&
            (fname.rfind(prefix, 0) == 0 ||
             fname.find("bwt_") == 0 ||
             fname.find("alphabet") != std::string::npos)) {
            idx_bytes += entry.file_size();
        }
    }

    std::cout << "tamaño índice (archivos sdsl): "
              << idx_bytes << " bytes";
    if (idx_bytes > 0) {
        // No conocemos n aquí — imprimimos solo bytes.
        std::cout << " (" << std::fixed << std::setprecision(2)
                  << idx_bytes / (1024.0 * 1024.0) << " MB)";
    }
    std::cout << "\n";

    if (mode == "info") return 0;

    // ── Leer patrones ─────────────────────────────────────────────────────────
    std::vector<std::string> patterns;
    {
        std::ifstream pf(patterns_path);
        if (!pf) {
            std::cerr << "No se puede abrir: " << patterns_path << "\n";
            return 3;
        }
        std::string line;
        while (std::getline(pf, line))
            if (!line.empty()) patterns.push_back(std::move(line));
    }
    if (patterns.empty()) {
        std::cerr << "No se encontraron patrones en " << patterns_path << "\n";
        return 4;
    }
    std::cout << "patrones : " << patterns.size() << "\n\n";

    // ── Warmup ────────────────────────────────────────────────────────────────
    for (std::size_t i = 0; i < std::min(patterns.size(), std::size_t{10}); ++i)
        (void)sr_loc.locate_extremal(patterns[i]);

    // ── Benchmark locate_extremal ─────────────────────────────────────────────
    std::size_t hits    = 0;
    std::size_t no_hits = 0;
    std::size_t d       = 0;  // sink para evitar dead-code elimination

    auto tb0 = Clock::now();
    for (const auto& p : patterns) {
        const auto [mn, mx] = sr_loc.locate_extremal(p);
        if (mn == SIZE_MAX) { ++no_hits; }
        else                { ++hits; d += mn + mx; }
    }
    double total_us = std::chrono::duration_cast<us>(Clock::now() - tb0).count();
    (void)d;

    const double us_per_query = total_us / static_cast<double>(patterns.size());

    std::cout << "=== locate_extremal (sr-index, sr=" << sr << ") ===\n";
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  patrones con ocurrencias : " << hits    << "\n";
    std::cout << "  patrones sin ocurrencias : " << no_hits << "\n";
    std::cout << "  tiempo total             : "
              << total_us / 1e6 << " s\n";
    std::cout << "  µs/query                 : " << us_per_query << "\n";

    return 0;
}
