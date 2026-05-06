/**
 * bench_special.cpp
 *
 * Valida y mide el overhead del path special en LZ77Index sobre un corpus FASTA.
 *
 * Uso:
 *   ./bench_special <fasta_file> [num_patterns] [pattern_len]
 *
 * Salida:
 *   - z (frases LZ77), n (largo texto)
 *   - Mismatches entre idx.count y bf_primary_count
 *   - Cuántos patrones son detectados SOLO por el path special
 *   - Tiempo promedio por query count()
 */

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "index.hpp"
#include "lz77/parser.hpp"

using namespace lz77tax;
using Clock = std::chrono::steady_clock;

// ─── Leer FASTA, concatenar secuencias (sin headers ni newlines) ─────────────
static std::string read_fasta(const std::string& path) {
    std::ifstream f(path);
    if (!f) { std::cerr << "No se pudo abrir: " << path << "\n"; std::exit(1); }
    std::string text, line;
    while (std::getline(f, line)) {
        if (!line.empty() && line[0] == '>') continue;
        for (char c : line) if (c != '\r') text += c;
    }
    return text;
}

// ─── Brute-force primary count (crossings + special) ─────────────────────────
static size_t bf_primary_count(const std::string& text,
                                const std::string& pattern,
                                const LZ77Parsing& phrases) {
    const size_t m = pattern.size();
    if (m < 2) return 0;
    size_t count = 0;
    for (size_t p = 0; p + m <= text.size(); ++p) {
        if (text.compare(p, m, pattern) != 0) continue;
        for (size_t k = 0; k + 1 < phrases.size(); ++k) {
            const size_t b = phrases[k + 1].start_pos;
            if (b >= text.size()) continue;
            if (p < b && b <= p + m - 1) ++count;
        }
        for (size_t k = 0; k + 1 < phrases.size(); ++k) {
            const size_t end_k = phrases[k].start_pos + phrases[k].length;
            if (end_k >= text.size()) continue;
            if (phrases[k].length + 1 >= m && p + m - 1 == end_k) ++count;
        }
    }
    return count;
}

// ─── Brute-force special-only count (end-aligned, no crossings) ──────────────
static size_t bf_special_only_count(const std::string& text,
                                    const std::string& pattern,
                                    const LZ77Parsing& phrases) {
    const size_t m = pattern.size();
    if (m < 2) return 0;
    size_t count = 0;
    for (size_t p = 0; p + m <= text.size(); ++p) {
        if (text.compare(p, m, pattern) != 0) continue;
        bool has_crossing = false;
        for (size_t k = 0; k + 1 < phrases.size(); ++k) {
            const size_t b = phrases[k + 1].start_pos;
            if (b >= text.size()) continue;
            if (p < b && b <= p + m - 1) { has_crossing = true; break; }
        }
        bool has_special = false;
        for (size_t k = 0; k + 1 < phrases.size(); ++k) {
            const size_t end_k = phrases[k].start_pos + phrases[k].length;
            if (end_k >= text.size()) continue;
            if (phrases[k].length + 1 >= m && p + m - 1 == end_k) { has_special = true; break; }
        }
        if (has_special && !has_crossing) ++count;
    }
    return count;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Uso: bench_special <fasta> [npatterns=200] [plen=8]\n";
        return 1;
    }
    const std::string fasta_path = argv[1];
    const int npatterns = (argc >= 3) ? std::atoi(argv[2]) : 200;
    const size_t plen   = (argc >= 4) ? (size_t)std::atoi(argv[3]) : 8;

    std::cout << "Leyendo " << fasta_path << " ...\n";
    const std::string text = read_fasta(fasta_path);
    std::cout << "  n = " << text.size() << " chars\n";

    std::cout << "Construyendo índice ...\n";
    auto t0 = Clock::now();
    LZ77Index idx;
    idx.build(text);
    double build_ms = std::chrono::duration<double,std::milli>(Clock::now()-t0).count();

    std::cout << "  z = " << idx.phrase_count() << " frases\n";
    std::cout << "  grid_points = " << idx.grid_points() << "\n";
    std::cout << "  build_time = " << (int)build_ms << " ms\n\n";

    // ── Generar patrones aleatorios ─────────────────────────────────────────
    std::mt19937 rng(12345);
    std::vector<std::string> patterns;
    patterns.reserve(npatterns);
    if (text.size() < plen) { std::cerr << "Texto demasiado corto\n"; return 1; }
    for (int i = 0; i < npatterns; ++i) {
        size_t pos = rng() % (text.size() - plen + 1);
        patterns.emplace_back(text.substr(pos, plen));
    }

    // ── Contar y validar ────────────────────────────────────────────────────
    int mismatches = 0;
    int special_only_detected = 0;   // patrones donde idx.count > 0 SOLO por special
    size_t total_count = 0;

    std::cout << "Validando " << npatterns << " patrones (plen=" << plen << ") ...\n";
    t0 = Clock::now();
    for (const auto& p : patterns) {
        const size_t got = idx.count(p);
        total_count += got;
        const size_t exp = bf_primary_count(text, p, idx.phrases());
        if (got != exp) {
            ++mismatches;
            std::cerr << "MISMATCH pattern='" << p
                      << "' got=" << got << " exp=" << exp << "\n";
        }
        // ¿es primario SOLO por special?
        const size_t sonly = bf_special_only_count(text, p, idx.phrases());
        if (sonly > 0) ++special_only_detected;
    }
    double query_ms = std::chrono::duration<double,std::milli>(Clock::now()-t0).count();

    std::cout << "  mismatches = " << mismatches << " / " << npatterns << "\n";
    std::cout << "  patrones con special-only = " << special_only_detected
              << " / " << npatterns << "\n";
    std::cout << "  total_count (sum primaries) = " << total_count << "\n";
    std::cout << "  query_time = " << query_ms << " ms total"
              << "  (" << (query_ms/npatterns*1000) << " µs/query)\n";

    return mismatches > 0 ? 1 : 0;
}
