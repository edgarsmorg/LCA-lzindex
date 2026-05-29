/**
 * measure_index — mide el tamaño en memoria del LZ77-index por componente
 * y compara la latencia de locate_min vs locate_extremal.
 *
 * Uso:
 *   measure_index <text_file> [--patterns=<file>]
 */

#include "index.hpp"
#include "lz77/grid.hpp"

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using namespace lz77tax;
using Clock = std::chrono::steady_clock;
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

    const size_t wt_bytes    = sdsl::size_in_bytes(grid.wt());
    const size_t bvfwd_bytes = sdsl::size_in_bytes(grid.bv_fwd());
    const size_t bvrev_bytes = sdsl::size_in_bytes(grid.bv_rev());

    const auto bd_min = grid.wm_min_rmq().size_breakdown();
    const auto bd_max = grid.wm_max_rmq().size_breakdown();
    const size_t wm_min_total = bd_min.bv + bd_min.rank_sel + bd_min.rmq;
    const size_t wm_max_total = bd_max.bv + bd_max.rank_sel + bd_max.rmq;

    print_row("wt_int (grilla / count)",  wt_bytes,    n);
    print_row("sd_vector fwd",            bvfwd_bytes, n);
    print_row("sd_vector rev",            bvrev_bytes, n);
    print_row("wm_min_rmq (bitvector)",   bd_min.bv,       n);
    print_row("wm_min_rmq (rank+select)", bd_min.rank_sel, n);
    print_row("wm_min_rmq (rmq BP)",      bd_min.rmq,      n);
    print_row("wm_max_rmq (bitvector)",   bd_max.bv,       n);
    print_row("wm_max_rmq (rank+select)", bd_max.rank_sel, n);
    print_row("wm_max_rmq (rmq BP)",      bd_max.rmq,      n);

    const size_t csa_f_bytes = idx.csa_fwd_bytes();
    const size_t csa_r_bytes = idx.csa_rev_bytes();
    const size_t csa_total   = csa_f_bytes + csa_r_bytes;
    const size_t grid_total  = wt_bytes + bvfwd_bytes + bvrev_bytes
                             + wm_min_total + wm_max_total;
    const size_t total       = grid_total + csa_total;

    std::cout << "  " << std::string(50, '-') << "\n";
    print_row("  subtotal grilla",  grid_total,  n);
    print_row("  CSA forward",      csa_f_bytes, n);
    print_row("  CSA reverse",      csa_r_bytes, n);
    print_row("  subtotal CSA",     csa_total,   n);
    print_row("  TOTAL",            total,        n);
    std::cout << "\n";
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "uso: measure_index <text_file> [--patterns=<file>]\n";
        return 1;
    }

    std::string patterns_path;
    for (int i = 2; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg.rfind("--patterns=", 0) == 0) {
            patterns_path = arg.substr(11);
        } else {
            std::cerr << "argumento desconocido: " << arg << "\n";
            return 1;
        }
    }

    std::cout << "Leyendo texto... " << std::flush;
    std::ifstream f(argv[1]);
    if (!f) { std::cerr << "\nno se puede abrir: " << argv[1] << "\n"; return 1; }
    const std::string text((std::istreambuf_iterator<char>(f)), {});
    if (text.empty()) { std::cerr << "\ntexto vacío\n"; return 1; }

    const size_t n = text.size();
    std::cout << "OK\n";
    std::cout << "archivo : " << argv[1] << "\n";
    std::cout << "n       : " << n << " bytes ("
              << std::fixed << std::setprecision(2) << n / (1024.0 * 1024.0) << " MB)\n\n";

    std::cout << "Construyendo LZ77Index (WmMinRmq)..." << std::flush;
    auto t0 = Clock::now();
    LZ77Index idx;
    idx.build(text);
    std::cout << " " << std::fixed << std::setprecision(1)
              << std::chrono::duration<double>(Clock::now() - t0).count() << "s\n\n";

    std::cout << "z (frases LZ77)  = " << idx.phrase_count() << "\n";
    std::cout << "z/n              = " << std::fixed << std::setprecision(4)
              << (double)idx.phrase_count() / n << "\n";
    std::cout << "puntos en grilla = " << idx.grid_points() << "\n\n";

    print_index_stats(idx, n);

    // Benchmark locate_min vs locate_extremal
    std::vector<std::string> patterns;
    std::string patterns_source = "aleatorios (8–32 bp)";
    if (!patterns_path.empty()) {
        std::ifstream pf(patterns_path);
        if (!pf) { std::cerr << "no se puede abrir patterns: " << patterns_path << "\n"; return 1; }
        std::string line;
        while (std::getline(pf, line))
            if (!line.empty()) patterns.push_back(line);
        patterns_source = patterns_path + " (" + std::to_string(patterns.size()) + " patrones)";
    } else {
        std::mt19937 rng(42);
        std::uniform_int_distribution<size_t> pos_dis(0, n - 33);
        std::uniform_int_distribution<size_t> len_dis(8, 32);
        for (size_t i = 0; i < 1000; ++i)
            patterns.push_back(text.substr(pos_dis(rng), len_dis(rng)));
    }
    const size_t N_PAT = patterns.size();
    std::cout << "=== Bench: locate_min vs locate_extremal ===\n";
    std::cout << "patrones: " << patterns_source << "\n";

    for (const auto& p : patterns) { (void)idx.locate_min(p); (void)idx.locate_extremal(p); }

    auto tb0 = Clock::now();
    size_t d1 = 0; for (const auto& p : patterns) d1 += idx.locate_min(p);
    double us_min = std::chrono::duration_cast<us>(Clock::now() - tb0).count() / (double)N_PAT;

    auto tb1 = Clock::now();
    size_t d2 = 0; for (const auto& p : patterns) d2 += idx.locate_extremal(p).first;
    double us_ext = std::chrono::duration_cast<us>(Clock::now() - tb1).count() / (double)N_PAT;

    (void)d1; (void)d2;
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  locate_min        : " << us_min << " µs/query\n";
    std::cout << "  locate_extremal   : " << us_ext << " µs/query\n";
    if (us_ext > 0.01)
        std::cout << "  speedup min/ext   : " << us_ext / us_min << "x\n";

    return 0;
}
