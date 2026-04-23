/**
 * measure_index — mide el tamaño en memoria del LZ77-index por componente
 * y compara la latencia de locate_min vs locate_extremal.
 *
 * Uso:
 *   measure_index <text_file>
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

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "uso: measure_index <text_file>\n";
        return 1;
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

    // ── Construcción ──────────────────────────────────────────────────────────
    std::cout << "Construyendo LZ77Index..." << std::flush;
    auto t0 = Clock::now();
    LZ77Index idx;
    idx.build(text);
    auto t1 = Clock::now();
    double secs = std::chrono::duration<double>(t1 - t0).count();
    std::cout << " " << std::fixed << std::setprecision(1) << secs << "s\n\n";

    std::cout << "z (frases LZ77)  = " << idx.phrase_count() << "\n";
    std::cout << "z/n              = " << std::fixed << std::setprecision(4)
              << (double)idx.phrase_count() / n << "  (compresibilidad LZ)\n";
    std::cout << "puntos en grilla = " << idx.grid_points()  << "\n\n";

    // ── Tamaños ───────────────────────────────────────────────────────────────
    const Grid2D& grid = idx.grid();
    const auto bd = grid.wt_min_rmq().size_breakdown();

    const size_t wt_bytes     = sdsl::size_in_bytes(grid.wt());
    const size_t bvfwd_bytes  = sdsl::size_in_bytes(grid.bv_fwd());
    const size_t bvrev_bytes  = sdsl::size_in_bytes(grid.bv_rev());
    // text_pos_ es privado (sdsl::int_vector compacto); usamos cota conservadora
    const size_t tpos_approx  = idx.grid_points() * sizeof(size_t);
    const size_t wt_min_total = bd.bv + bd.rank_sel + bd.rmq;
    const size_t csa_f_bytes  = idx.csa_fwd_bytes();
    const size_t csa_r_bytes  = idx.csa_rev_bytes();

    const size_t grid_total = wt_bytes + bvfwd_bytes + bvrev_bytes
                            + tpos_approx + wt_min_total;
    const size_t csa_total  = csa_f_bytes + csa_r_bytes;
    const size_t total      = grid_total + csa_total;

    std::cout << "=== LZ77-Index: tamaño por componente ===\n";
    print_row("wt_int (grilla / count)",  wt_bytes,    n);
    print_row("sd_vector fwd",            bvfwd_bytes, n);
    print_row("sd_vector rev",            bvrev_bytes, n);
    print_row("text_pos[] (z×8B aprox)", tpos_approx, n);
    std::cout << "  --- WtMinRmq (locate_min) ---\n";
    print_row("  wt_min_rmq (bitvectors)", bd.bv,      n);
    print_row("  wt_min_rmq (rank+select)", bd.rank_sel, n);
    print_row("  wt_min_rmq (rmq BP)",    bd.rmq,     n);
    print_row("  wt_min_rmq total",       wt_min_total, n);
    std::cout << "  " << std::string(50, '-') << "\n";
    print_row("  subtotal grilla",        grid_total,  n);
    std::cout << "\n";
    print_row("csa_wt forward",           csa_f_bytes, n);
    print_row("csa_wt reverse",           csa_r_bytes, n);
    std::cout << "  " << std::string(50, '-') << "\n";
    print_row("  subtotal CSA",           csa_total,   n);
    std::cout << "\n";
    print_row("TOTAL",                    total,       n);

    std::cout << "\nNOTA: En producción los CSA (~"
              << std::fixed << std::setprecision(0)
              << csa_total * 100.0 / total
              << "%) se reemplazan por ropebwt3. Proyección sin CSA:\n";
    print_row("  PROYECCIÓN (sin CSA)", grid_total, n);

    // ── Benchmark locate_min vs locate_extremal ───────────────────────────────
    std::cout << "\n=== Bench: locate_min vs locate_extremal (1000 patrones) ===\n";

    // Generar 1000 patrones aleatorios de longitud 8..32
    const size_t N_PATTERNS = 1000;
    std::mt19937 rng(42);
    std::uniform_int_distribution<size_t> pos_dis(0, n - 33);
    std::uniform_int_distribution<size_t> len_dis(8, 32);

    std::vector<std::string> patterns;
    patterns.reserve(N_PATTERNS);
    for (size_t i = 0; i < N_PATTERNS; ++i) {
        const size_t pos = pos_dis(rng);
        const size_t len = len_dis(rng);
        patterns.push_back(text.substr(pos, len));
    }

    // Warmup
    for (const auto& p : patterns) { (void)idx.locate_min(p); }

    // Bench locate_min
    auto t2 = Clock::now();
    size_t dummy1 = 0;
    for (const auto& p : patterns) dummy1 += idx.locate_min(p);
    auto t3 = Clock::now();
    double us_min = std::chrono::duration_cast<us>(t3 - t2).count() / (double)N_PATTERNS;

    // Bench locate_extremal
    auto t4 = Clock::now();
    size_t dummy2 = 0;
    for (const auto& p : patterns) dummy2 += idx.locate_extremal(p).first;
    auto t5 = Clock::now();
    double us_ext = std::chrono::duration_cast<us>(t5 - t4).count() / (double)N_PATTERNS;

    (void)dummy1; (void)dummy2;  // evitar que el compilador las elimine

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  locate_min      : " << us_min << " µs/query\n";
    std::cout << "  locate_extremal : " << us_ext << " µs/query\n";
    if (us_ext > 0.01)
        std::cout << "  speedup         : " << us_ext / us_min << "x\n";

    return 0;
}
