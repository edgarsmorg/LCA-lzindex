/**
 * measure_index — mide el tamaño en memoria del LZ77-index por componente
 * y compara la latencia de locate_min con WtMinRmq vs WmMinRmq.
 *
 * Uso:
 *   measure_index <text_file> [--patterns=<file>] [--variant=wt|wm|both]
 *
 *   --variant=wm   (defecto) solo construye WmMinRmq — rápido, para benchmarks
 *   --variant=wt   solo construye WtMinRmq
 *   --variant=both construye ambas variantes para comparación de espacio
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

static void print_index_stats(const LZ77Index& idx, RmqVariant variant,
                               const std::string& label, size_t n) {
    const Grid2D& grid = idx.grid();
    std::cout << "=== " << label << ": tamaño por componente ===\n";

    const size_t wt_bytes    = sdsl::size_in_bytes(grid.wt());
    const size_t bvfwd_bytes = sdsl::size_in_bytes(grid.bv_fwd());
    const size_t bvrev_bytes = sdsl::size_in_bytes(grid.bv_rev());
    const size_t tpos_approx = idx.grid_points() * sizeof(size_t);

    print_row("wt_int (grilla / count)",  wt_bytes,    n);
    print_row("sd_vector fwd",            bvfwd_bytes, n);
    print_row("sd_vector rev",            bvrev_bytes, n);
    print_row("text_pos[] (z×8B aprox)", tpos_approx, n);

    size_t rmq_total = 0;
    if (variant == RmqVariant::Wt) {
        const auto bd = grid.wt_min_rmq().size_breakdown();
        print_row("  wt_min_rmq (bitvectors)",  bd.bv,       n);
        print_row("  wt_min_rmq (rank+select)", bd.rank_sel, n);
        print_row("  wt_min_rmq (rmq BP)",      bd.rmq,      n);
        rmq_total = 2 * (bd.bv + bd.rank_sel + bd.rmq);  // min + max
        print_row("  wt_min_rmq ×2 (min+max)",  rmq_total,   n);
    } else {
        const auto bd = grid.wm_min_rmq().size_breakdown();
        print_row("  wm_min_rmq (bitvector)",   bd.bv,       n);
        print_row("  wm_min_rmq (rank+select)", bd.rank_sel, n);
        print_row("  wm_min_rmq (rmq BP)",      bd.rmq,      n);
        rmq_total = 2 * (bd.bv + bd.rank_sel + bd.rmq);  // min + max
        print_row("  wm_min_rmq ×2 (min+max)",  rmq_total,   n);
    }

    const size_t csa_f_bytes = idx.csa_fwd_bytes();
    const size_t csa_r_bytes = idx.csa_rev_bytes();
    const size_t csa_total   = csa_f_bytes + csa_r_bytes;
    const size_t grid_total  = wt_bytes + bvfwd_bytes + bvrev_bytes + tpos_approx + rmq_total;
    const size_t total       = grid_total + csa_total;

    std::cout << "  " << std::string(50, '-') << "\n";
    print_row("  subtotal grilla",  grid_total,  n);
    print_row("  subtotal CSA",     csa_total,   n);
    print_row("  TOTAL",            total,        n);
    std::cout << "\n";
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "uso: measure_index <text_file> [--patterns=<file>] [--variant=wt|wm|both]\n";
        return 1;
    }

    std::string patterns_path;
    bool run_both = false;
    RmqVariant single_variant = RmqVariant::Wm;

    for (int i = 2; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg.rfind("--patterns=", 0) == 0) {
            patterns_path = arg.substr(11);
        } else if (arg == "--variant=wt") {
            single_variant = RmqVariant::Wt;
        } else if (arg == "--variant=wm") {
            single_variant = RmqVariant::Wm;
        } else if (arg == "--variant=both") {
            run_both = true;
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

    if (run_both) {
        // ── Modo comparación: construye ambas variantes ────────────────────────
        std::cout << "Construyendo LZ77Index (WmMinRmq)..." << std::flush;
        auto t0 = Clock::now();
        LZ77Index idx_wm;
        idx_wm.build(text, RmqVariant::Wm);
        std::cout << " " << std::fixed << std::setprecision(1)
                  << std::chrono::duration<double>(Clock::now() - t0).count() << "s\n";

        std::cout << "Construyendo LZ77Index (WtMinRmq)..." << std::flush;
        auto t1 = Clock::now();
        LZ77Index idx_wt;
        idx_wt.build(text, RmqVariant::Wt);
        std::cout << " " << std::fixed << std::setprecision(1)
                  << std::chrono::duration<double>(Clock::now() - t1).count() << "s\n\n";

        std::cout << "z (frases LZ77)  = " << idx_wm.phrase_count() << "\n";
        std::cout << "z/n              = " << std::fixed << std::setprecision(4)
                  << (double)idx_wm.phrase_count() / n << "\n";
        std::cout << "puntos en grilla = " << idx_wm.grid_points() << "\n\n";

        print_index_stats(idx_wt, RmqVariant::Wt, "LZ77-Index (WtMinRmq, legacy)", n);
        print_index_stats(idx_wm, RmqVariant::Wm, "LZ77-Index (WmMinRmq, producción)", n);

        // Benchmark ambas
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
        std::cout << "=== Bench: WtMinRmq vs WmMinRmq ===\n";
        std::cout << "patrones: " << patterns_source << "\n";

        for (const auto& p : patterns) { (void)idx_wt.locate_min(p); (void)idx_wm.locate_min(p); }

        auto tb0 = Clock::now();
        size_t d1 = 0; for (const auto& p : patterns) d1 += idx_wt.locate_min(p);
        double us_wt = std::chrono::duration_cast<us>(Clock::now() - tb0).count() / (double)N_PAT;

        auto tb1 = Clock::now();
        size_t d2 = 0; for (const auto& p : patterns) d2 += idx_wm.locate_min(p);
        double us_wm = std::chrono::duration_cast<us>(Clock::now() - tb1).count() / (double)N_PAT;

        (void)d1; (void)d2;
        std::cout << std::fixed << std::setprecision(2);
        std::cout << "  locate_min (WtMinRmq): " << us_wt << " µs/query\n";
        std::cout << "  locate_min (WmMinRmq): " << us_wm << " µs/query\n";

    } else {
        // ── Modo single-variant (defecto: Wm) ─────────────────────────────────
        const std::string vname = (single_variant == RmqVariant::Wm) ? "WmMinRmq" : "WtMinRmq";
        std::cout << "Construyendo LZ77Index (" << vname << ")..." << std::flush;
        auto t0 = Clock::now();
        LZ77Index idx;
        idx.build(text, single_variant);
        std::cout << " " << std::fixed << std::setprecision(1)
                  << std::chrono::duration<double>(Clock::now() - t0).count() << "s\n\n";

        std::cout << "z (frases LZ77)  = " << idx.phrase_count() << "\n";
        std::cout << "z/n              = " << std::fixed << std::setprecision(4)
                  << (double)idx.phrase_count() / n << "\n";
        std::cout << "puntos en grilla = " << idx.grid_points() << "\n\n";

        print_index_stats(idx, single_variant, "LZ77-Index (" + vname + ")", n);

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
        std::cout << "=== Bench: locate_min (" << vname << ") vs locate_extremal ===\n";
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
        std::cout << "  locate_min (" << vname << "): " << us_min << " µs/query\n";
        std::cout << "  locate_extremal             : " << us_ext << " µs/query\n";
        if (us_ext > 0.01)
            std::cout << "  speedup locate_min/ext      : " << us_ext / us_min << "x\n";

        std::cout << "\nNOTA: usar --variant=both para comparar Wt vs Wm en espacio y tiempo.\n";
    }

    return 0;
}
