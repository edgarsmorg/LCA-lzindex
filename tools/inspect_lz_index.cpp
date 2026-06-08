#include "bench_common.hpp"

#include <filesystem>
#include <iostream>

using lz77tax::LZ77Index;
namespace fs = std::filesystem;

static void usage(const char* argv0) {
    std::cerr << "usage: " << argv0 << " <text_file> <index_prefix> [--name=<data_name>] [--csv=<path>]\n";
}

int main(int argc, char** argv) {
    if (argc < 3) { usage(argv[0]); return 1; }

    const fs::path text_path = argv[1];
    const fs::path prefix = argv[2];
    std::string dataset = bench::basename(text_path);
    fs::path csv_path;
    for (int i = 3; i < argc; ++i) {
        const std::string arg = argv[i];
        if (auto v = bench::option_value(arg, "name"); !v.empty()) dataset = v;
        else if (auto v = bench::option_value(arg, "csv"); !v.empty()) csv_path = v;
        else { std::cerr << "unknown argument: " << arg << "\n"; usage(argv[0]); return 1; }
    }

    const auto n_bytes = fs::file_size(text_path);
    LZ77Index idx;
    idx.load(prefix);
    const auto z = idx.phrase_count();
    const auto grid_bytes = bench::lz_grid_bytes(idx);
    const auto csa_bytes = idx.csa_fwd_bytes() + idx.csa_rev_bytes();
    const auto total = grid_bytes + csa_bytes;
    const auto disk = bench::sum_existing_files(bench::lz_files(prefix));
    const double bpc = n_bytes ? total * 8.0 / static_cast<double>(n_bytes) : 0.0;

    const std::string row =
        bench::csv_quote(dataset) + ",lz77,0," +
        std::to_string(n_bytes) + "," + std::to_string(z) + "," +
        std::to_string(total) + "," + std::to_string(disk) + "," +
        std::to_string(grid_bytes) + "," + std::to_string(csa_bytes) + "," +
        bench::seconds(bpc);

    bench::append_csv_row(csv_path,
        "dataset,index,s,n_bytes,z,serialized_bytes,disk_bytes,grid_bytes,csa_bytes,bpc", row);
    std::cout << row << "\n";
    return 0;
}
