#include <filesystem>
#include <iostream>
#include <string>

#include <sdsl/io.hpp>

#include "sr-index/config.h"
#include "sr-index/construct.h"
#include "sr-index/sr_index.h"

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <data.txt> <index.sdsl> [sr]\n";
        return 1;
    }

    const std::string data_path = argv[1];
    const std::string index_path = argv[2];
    const std::size_t sr = (argc >= 4) ? static_cast<std::size_t>(std::stoul(argv[3])) : 16;

    sri::Config config(
        data_path,
        std::filesystem::path(index_path).parent_path(),
        sri::SDSL_LIBDIVSUFSORT,
        false,
        8);

    sri::SrIndexValidArea<> index(sr);
    sri::construct(index, data_path, config);
    sdsl::store_to_file(index, index_path);
    std::cout << "SR index construido y guardado en: " << index_path << "\n";

    return 0;
}
