#include "bench_common.hpp"

#include <filesystem>
#include <iostream>

using lz77tax::LZ77Index;
namespace fs = std::filesystem;

static void usage(const char* argv0) {
    std::cerr << "usage: " << argv0 << " <text_file> <out_prefix> [--name=<data_name>] [--csv=<path>]\n";
}

int main(int argc, char** argv) {
    if (argc < 3) { usage(argv[0]); return 1; }

    const fs::path text_path = argv[1];
    const fs::path out_prefix = argv[2];
    std::string dataset = bench::basename(text_path);
    fs::path csv_path;
    for (int i = 3; i < argc; ++i) {
        const std::string arg = argv[i];
        if (auto v = bench::option_value(arg, "name"); !v.empty()) dataset = v;
        else if (auto v = bench::option_value(arg, "csv"); !v.empty()) csv_path = v;
        else { std::cerr << "unknown argument: " << arg << "\n"; usage(argv[0]); return 1; }
    }

    const auto text = bench::read_text_file(text_path);
    fs::create_directories(out_prefix.parent_path().empty() ? fs::path(".") : out_prefix.parent_path());

    auto t0 = bench::Clock::now();
    std::string status = "ok";
    std::size_t z = 0;
    try {
        LZ77Index idx;
        idx.build(text);
        idx.build_oracle();          // accesor → el índice guardado es self-index
        z = idx.phrase_count();
        idx.save(out_prefix);
    } catch (const std::exception& e) {
        status = std::string("error:") + e.what();
    }
    const double build_s = std::chrono::duration<double>(bench::Clock::now() - t0).count();
    const double z_per_n = text.empty() ? 0.0 : z / static_cast<double>(text.size());

    const std::string row =
        bench::csv_quote(dataset) + ",lz77,0," +
        std::to_string(text.size()) + "," + std::to_string(z) + "," +
        bench::seconds(z_per_n) + "," + bench::seconds(build_s) + "," +
        bench::csv_quote(out_prefix.string()) + "," + bench::csv_quote(status);

    bench::append_csv_row(csv_path,
        "dataset,index,s,n_bytes,z,z_per_n,build_seconds,out_path,status", row);
    std::cout << row << "\n";
    return status == "ok" ? 0 : 2;
}
