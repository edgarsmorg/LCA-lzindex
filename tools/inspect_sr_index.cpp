#include "bench_common.hpp"
#include "baseline/sr_index_locator.hpp"

#include <filesystem>
#include <iostream>

using lz77tax::SrIndexLocator;
namespace fs = std::filesystem;

static void usage(const char* argv0) {
    std::cerr << "usage: " << argv0
              << " <text_file> <sr_dir> [--s=16] [--name=<data_name>] [--csv=<path>]\n";
}

int main(int argc, char** argv) {
    if (argc < 3) { usage(argv[0]); return 1; }

    const fs::path text_path = argv[1];
    const fs::path sr_dir = argv[2];
    std::size_t s = 16;
    std::string dataset = bench::basename(text_path);
    fs::path csv_path;

    for (int i = 3; i < argc; ++i) {
        const std::string arg = argv[i];
        if (auto v = bench::option_value(arg, "s"); !v.empty()) s = std::stoull(v);
        else if (auto v = bench::option_value(arg, "name"); !v.empty()) dataset = v;
        else if (auto v = bench::option_value(arg, "csv"); !v.empty()) csv_path = v;
        else { std::cerr << "unknown argument: " << arg << "\n"; usage(argv[0]); return 1; }
    }

    const auto n_bytes = fs::file_size(text_path);
    // Valida que el índice es cargable pero no usamos size_in_bytes(): el
    // sr-index almacena su data en GenericStorage y serialize() solo emite
    // metadata, subestimando el tamaño real por ~5×.
    SrIndexLocator idx;
    idx.load(dataset, sr_dir.string(), s);
    const auto serialized = bench::sr_index_bytes(sr_dir, dataset);
    const auto disk = bench::sum_files_with_extension(sr_dir, ".sdsl");
    const double bpc = n_bytes ? serialized * 8.0 / static_cast<double>(n_bytes) : 0.0;

    const std::string row =
        bench::csv_quote(dataset) + ",sr," + std::to_string(s) + "," +
        std::to_string(n_bytes) + ",," + std::to_string(serialized) + "," +
        std::to_string(disk) + ",,," + bench::seconds(bpc);

    bench::append_csv_row(csv_path,
        "dataset,index,s,n_bytes,z,serialized_bytes,disk_bytes,grid_bytes,csa_bytes,bpc", row);
    std::cout << row << "\n";
    return 0;
}
