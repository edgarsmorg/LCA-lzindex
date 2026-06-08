#include "bench_common.hpp"
#include "baseline/sr_index_locator.hpp"

#include <filesystem>
#include <iostream>

using lz77tax::SrIndexLocator;
namespace fs = std::filesystem;

static void usage(const char* argv0) {
    std::cerr << "usage: " << argv0
              << " <text_file> <out_dir> [--s=16] [--name=<data_name>] [--csv=<path>]\n";
}

int main(int argc, char** argv) {
    if (argc < 3) { usage(argv[0]); return 1; }

    const fs::path text_path = argv[1];
    const fs::path out_dir = argv[2];
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
    fs::create_directories(out_dir);
    fs::path build_text_path = text_path;
    if (dataset != bench::basename(text_path)) {
        build_text_path = out_dir / dataset;
        fs::copy_file(text_path, build_text_path, fs::copy_options::overwrite_existing);
    }

    auto t0 = bench::Clock::now();
    std::string status = "ok";
    try {
        SrIndexLocator idx;
        idx.build(build_text_path.string(), out_dir.string(), s);
    } catch (const std::exception& e) {
        status = std::string("error:") + e.what();
    }
    const double build_s = std::chrono::duration<double>(bench::Clock::now() - t0).count();

    const std::string row =
        bench::csv_quote(dataset) + ",sr," + std::to_string(s) + "," +
        std::to_string(n_bytes) + ",,," + bench::seconds(build_s) + "," +
        bench::csv_quote(out_dir.string()) + "," + bench::csv_quote(status);

    bench::append_csv_row(csv_path,
        "dataset,index,s,n_bytes,z,z_per_n,build_seconds,out_path,status", row);
    std::cout << row << "\n";
    return status == "ok" ? 0 : 2;
}
