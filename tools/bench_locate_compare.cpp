#include "bench_common.hpp"
#include "baseline/sr_index_locator.hpp"

#include <climits>
#include <filesystem>
#include <iostream>

using lz77tax::LZ77Index;
using lz77tax::SrIndexLocator;
namespace fs = std::filesystem;

static void usage(const char* argv0) {
    std::cerr << "usage: " << argv0
              << " <text_file> <patterns_file> <lz_prefix> <sr_dir>"
              << " [--s=16] [--name=<data_name>] [--reps=3] [--csv=<path>]\n";
}

template <class Index>
static std::string run_one(const std::string& dataset, const std::string& index_name,
                           std::size_t s, std::size_t rep,
                           const std::vector<std::string>& patterns,
                           const Index& idx) {
    std::size_t hits = 0;
    std::size_t no_hits = 0;
    std::size_t checksum = 1469598103934665603ULL;

    const auto t0 = bench::Clock::now();
    for (const auto& p : patterns) {
        const auto [mn, mx] = idx.locate_extremal(p);
        if (mn == SIZE_MAX) {
            ++no_hits;
            checksum ^= p.size() + 0x9e3779b97f4a7c15ULL;
        } else {
            ++hits;
            checksum ^= mn + 0x9e3779b97f4a7c15ULL + (checksum << 6) + (checksum >> 2);
            checksum ^= mx + 0x9e3779b97f4a7c15ULL + (checksum << 6) + (checksum >> 2);
        }
    }
    const double total_s = std::chrono::duration<double>(bench::Clock::now() - t0).count();
    const double us_per_query = total_s * 1e6 / static_cast<double>(patterns.size());

    return bench::csv_quote(dataset) + "," + index_name + "," + std::to_string(s) + "," +
        std::to_string(rep) + "," + std::to_string(patterns.size()) + "," +
        std::to_string(hits) + "," + std::to_string(no_hits) + "," +
        bench::seconds(total_s) + "," + bench::seconds(us_per_query) + "," +
        std::to_string(checksum);
}

int main(int argc, char** argv) {
    if (argc < 5) { usage(argv[0]); return 1; }

    const fs::path text_path = argv[1];
    const fs::path patterns_path = argv[2];
    const fs::path lz_prefix = argv[3];
    const fs::path sr_dir = argv[4];
    std::size_t s = 16;
    std::size_t reps = 3;
    std::string dataset = bench::basename(text_path);
    fs::path csv_path;

    for (int i = 5; i < argc; ++i) {
        const std::string arg = argv[i];
        if (auto v = bench::option_value(arg, "s"); !v.empty()) s = std::stoull(v);
        else if (auto v = bench::option_value(arg, "name"); !v.empty()) dataset = v;
        else if (auto v = bench::option_value(arg, "reps"); !v.empty()) reps = std::stoull(v);
        else if (auto v = bench::option_value(arg, "csv"); !v.empty()) csv_path = v;
        else { std::cerr << "unknown argument: " << arg << "\n"; usage(argv[0]); return 1; }
    }

    const auto patterns = bench::read_patterns(patterns_path);
    LZ77Index lz;
    lz.load(lz_prefix);
    SrIndexLocator sr;
    sr.load(dataset, sr_dir.string(), s);

    for (std::size_t i = 0; i < std::min<std::size_t>(patterns.size(), 10); ++i) {
        (void)lz.locate_extremal(patterns[i]);
        (void)sr.locate_extremal(patterns[i]);
    }

    const std::string header =
        "dataset,index,s,rep,n_patterns,hits,no_hits,total_seconds,us_per_query,checksum";
    for (std::size_t rep = 0; rep < reps; ++rep) {
        const auto lz_row = run_one(dataset, "lz77", 0, rep, patterns, lz);
        const auto sr_row = run_one(dataset, "sr", s, rep, patterns, sr);
        bench::append_csv_row(csv_path, header, lz_row);
        bench::append_csv_row(csv_path, header, sr_row);
        std::cout << lz_row << "\n" << sr_row << "\n";
    }
    return 0;
}
