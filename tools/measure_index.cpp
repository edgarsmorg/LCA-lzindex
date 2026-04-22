/**
 * measure_index — mide el tamaño en memoria del LZ77-index por componente.
 *
 * Construye el LZ77Index sobre el texto dado y reporta:
 *   - z  (número de frases LZ77)
 *   - Tamaño en bytes y bpc de cada componente
 *   - Total del índice
 *
 * Uso:
 *   measure_index <text_file>
 *
 * Ejemplo:
 *   ./build/measure_index data/Escherichia_Coli
 *   ./build/measure_index data/synthetic_1mb/synthetic_1mb.txt
 */

#include "index.hpp"
#include "lz77/grid.hpp"

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace lz77tax;
using Clock = std::chrono::steady_clock;

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

    std::cout << "  " << std::left  << std::setw(26) << label
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

    // ── Tamaños (reutiliza CSA ya construidos dentro del índice) ──────────────
    const Grid2D& grid = idx.grid();
    const size_t wt_bytes    = sdsl::size_in_bytes(grid.wt());
    const size_t bvfwd_bytes = sdsl::size_in_bytes(grid.bv_fwd());
    const size_t bvrev_bytes = sdsl::size_in_bytes(grid.bv_rev());
    const size_t tpos_bytes  = idx.grid_points() * sizeof(size_t);
    const size_t csa_f_bytes = idx.csa_fwd_bytes();
    const size_t csa_r_bytes = idx.csa_rev_bytes();

    const size_t grid_total = wt_bytes + bvfwd_bytes + bvrev_bytes + tpos_bytes;
    const size_t csa_total  = csa_f_bytes + csa_r_bytes;
    const size_t total      = grid_total + csa_total;

    std::cout << "=== LZ77-Index: tamaño por componente ===\n";
    print_row("wt_int (grilla)",    wt_bytes,    n);
    print_row("sd_vector fwd",      bvfwd_bytes, n);
    print_row("sd_vector rev",      bvrev_bytes, n);
    print_row("text_pos[] (z×8B)",  tpos_bytes,  n);
    std::cout << "  " << std::string(46, '-') << "\n";
    print_row("  subtotal grilla",  grid_total,  n);
    std::cout << "\n";
    print_row("csa_wt forward",     csa_f_bytes, n);
    print_row("csa_wt reverse",     csa_r_bytes, n);
    std::cout << "  " << std::string(46, '-') << "\n";
    print_row("  subtotal CSA",     csa_total,   n);
    std::cout << "\n";
    print_row("TOTAL",              total,       n);

    std::cout << "\nNOTA: En producción los CSA (~"
              << std::fixed << std::setprecision(0)
              << csa_total * 100.0 / total
              << "%) se reemplazan por el\n"
              << "      BWT comprimido de ropebwt3. Proyección sin CSA:\n";
    print_row("  PROYECCIÓN (sin CSA)", grid_total, n);

    return 0;
}
