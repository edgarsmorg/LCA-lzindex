/**
 * measure_index — mide el tamaño en memoria del LZ77-index por componente
 * y la latencia de locate_extremal.
 *
 * Uso:
 *   measure_index <text_file> [--patterns=<file>] [--name=<label>] [--csv=<path>]
 *
 * CSV (una fila, header auto-generado si el archivo no existe):
 *   dataset,n_bytes,z_phrases,z_per_n,build_time_s,
 *   total_bytes,total_bpc,csa_bpc,grid_bpc,
 *   n_patterns,locate_ext_us
 */

#include "bench_common.hpp"
#include "index.hpp"
#include "lz77/grid.hpp"

#include <filesystem>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace lz77tax;
using Clock = bench::Clock;
using us    = std::chrono::microseconds;

static void print_row(const std::string& label, size_t bytes, size_t n) {
    if (bytes == 0) return;
    double bpc  = (bytes * 8.0) / n;
    double kb   = bytes / 1024.0;
    double mb   = bytes / (1024.0 * 1024.0);
    std::string size_str;
    if (mb >= 1.0)
        size_str = std::to_string(static_cast<int>(mb + 0.5)) + " MB";
    else
        size_str = std::to_string(static_cast<int>(kb + 0.5)) + " KB";

    std::cout << "  " << std::left  << std::setw(30) << label
              << ": " << std::right << std::setw(7)  << size_str
              << "  (" << std::fixed << std::setprecision(2) << bpc << " bpc)\n";
}

static void print_index_stats(const LZ77Index& idx, size_t n) {
    const Grid2D& grid = idx.grid();
    std::cout << "=== LZ77-Index (WmMinRmq): tamaño por componente ===\n";

    const auto bd = grid.size_breakdown();

    const auto bd_min = grid.wm_min_rmq().size_breakdown();
    const auto bd_max = grid.wm_max_rmq().size_breakdown();

    print_row("wm_int (grilla compartida)", bd.wm,            n);
    print_row("sd_vector fwd",            bd.bv_fwd,          n);
    print_row("sd_vector rev",            bd.bv_rev,          n);
    print_row("rank_fwd + rank_rev",      bd.rank_fwd + bd.rank_rev, n);
    print_row("text_pos (int_vector)",    bd.text_pos,        n);
    print_row("phrase_total_len",         bd.phrase_total_len,n);
    print_row("wm_min_rmq (rmq BP)",      bd_min.rmq,         n);
    print_row("wm_max_rmq (rmq BP)",      bd_max.rmq,         n);

    const size_t trie_total  = idx.trie_bytes();
    const size_t grid_total  = bd.wm + bd.bv_fwd + bd.bv_rev
                             + bd.rank_fwd + bd.rank_rev
                             + bd.text_pos + bd.phrase_total_len
                             + bd.wm_min_rmq + bd.wm_max_rmq;
    const size_t total       = grid_total + trie_total;

    std::cout << "  " << std::string(50, '-') << "\n";
    print_row("  subtotal grilla",  grid_total,  n);
    print_row("  tries (SST+rev)",  trie_total,  n);
    print_row("  TOTAL",            total,        n);
    std::cout << "\n";
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "uso: measure_index <text_file> [--patterns=<file>] [--name=<label>] [--csv=<path>]\n";
        return 1;
    }

    std::string patterns_path;
    std::string dataset_name = bench::basename(std::filesystem::path(argv[1]));
    std::filesystem::path csv_path;
    for (int i = 2; i < argc; ++i) {
        std::string arg(argv[i]);
        if (auto v = bench::option_value(arg, "patterns"); !v.empty())
            patterns_path = v;
        else if (auto v = bench::option_value(arg, "name"); !v.empty())
            dataset_name = v;
        else if (auto v = bench::option_value(arg, "csv"); !v.empty())
            csv_path = v;
        else {
            std::cerr << "argumento desconocido: " << arg << "\n";
            return 1;
        }
    }

    std::cerr << "Leyendo texto... " << std::flush;
    std::string text;
    try { text = bench::read_text_file(argv[1]); }
    catch (const std::exception& e) { std::cerr << "\n" << e.what() << "\n"; return 1; }
    if (text.empty()) { std::cerr << "\ntexto vacío\n"; return 1; }

    const size_t n = text.size();
    std::cerr << "OK\n";
    std::cout << "archivo : " << argv[1] << "\n";
    std::cout << "n       : " << n << " bytes ("
              << std::fixed << std::setprecision(2) << n / (1024.0 * 1024.0) << " MB)\n\n";

    std::cerr << "Construyendo LZ77Index (WmMinRmq)..." << std::flush;
    auto t0 = Clock::now();
    LZ77Index idx;
    idx.build(text);
    const double build_s = std::chrono::duration<double>(Clock::now() - t0).count();
    std::cerr << " " << std::fixed << std::setprecision(1) << build_s << "s\n\n";

    std::cout << "z (frases LZ77)  = " << idx.phrase_count() << "\n";
    std::cout << "z/n              = " << std::fixed << std::setprecision(4)
              << (double)idx.phrase_count() / n << "\n";
    std::cout << "puntos en grilla = " << idx.grid_points() << "\n\n";

    print_index_stats(idx, n);

    // Benchmark locate_extremal
    std::vector<std::string> patterns;
    std::string patterns_source = "aleatorios (8–32 bp)";
    if (!patterns_path.empty()) {
        try { patterns = bench::read_patterns(patterns_path); }
        catch (const std::exception& e) { std::cerr << e.what() << "\n"; return 1; }
        patterns_source = patterns_path + " (" + std::to_string(patterns.size()) + " patrones)";
    } else {
        std::mt19937 rng(42);
        std::uniform_int_distribution<size_t> pos_dis(0, n > 33 ? n - 33 : 0);
        std::uniform_int_distribution<size_t> len_dis(8, 32);
        for (size_t i = 0; i < 1000; ++i)
            patterns.push_back(text.substr(pos_dis(rng), len_dis(rng)));
    }
    const size_t N_PAT = patterns.size();
    std::cout << "=== Bench: locate_extremal ===\n";
    std::cout << "patrones: " << patterns_source << "\n";

    // Warm-up
    for (const auto& p : patterns) (void)idx.locate_extremal(p);

    auto tb1 = Clock::now();
    size_t d2 = 0; for (const auto& p : patterns) d2 += idx.locate_extremal(p).first;
    double us_ext = std::chrono::duration_cast<us>(Clock::now() - tb1).count() / (double)N_PAT;

    (void)d2;
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  locate_extremal   : " << us_ext << " µs/query\n";

    // CSV output
    if (!csv_path.empty()) {
        const size_t z = idx.phrase_count();
        const size_t grid_bytes = bench::lz_grid_bytes(idx);
        const size_t trie_bytes = idx.trie_bytes();
        const size_t total      = grid_bytes + trie_bytes;
        const double total_bpc  = total      * 8.0 / n;
        const double trie_bpc   = trie_bytes * 8.0 / n;
        const double grid_bpc   = grid_bytes * 8.0 / n;

        std::ostringstream row;
        row << std::fixed << std::setprecision(6);
        row << bench::csv_quote(dataset_name) << ","
            << n          << ","
            << z          << ","
            << (double)z / n << ","
            << build_s    << ","
            << total      << ","
            << std::setprecision(4) << total_bpc << ","
            << trie_bpc   << ","
            << grid_bpc   << ","
            << N_PAT      << ","
            << std::setprecision(6) << us_ext;

        static constexpr const char* CSV_HEADER =
            "dataset,n_bytes,z_phrases,z_per_n,build_time_s,"
            "total_bytes,total_bpc,trie_bpc,grid_bpc,"
            "n_patterns,locate_ext_us";
        bench::append_csv_row(csv_path, CSV_HEADER, row.str());
    }

    return 0;
}
