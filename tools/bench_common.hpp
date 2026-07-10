#pragma once

#include "index.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace bench {

namespace fs = std::filesystem;
using Clock = std::chrono::steady_clock;

inline bool starts_with(const std::string& s, const std::string& prefix) {
    return s.rfind(prefix, 0) == 0;
}

inline std::string option_value(const std::string& arg, const std::string& name) {
    const std::string prefix = "--" + name + "=";
    if (!starts_with(arg, prefix)) return {};
    return arg.substr(prefix.size());
}

inline std::string basename(const fs::path& path) {
    return path.filename().string();
}

inline std::string read_text_file(const fs::path& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot open text file: " + path.string());
    return std::string(std::istreambuf_iterator<char>(in), {});
}

inline std::vector<std::string> read_patterns(const fs::path& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open patterns file: " + path.string());
    std::vector<std::string> patterns;
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (!line.empty()) patterns.push_back(std::move(line));
    }
    if (patterns.empty()) throw std::runtime_error("patterns file is empty: " + path.string());
    return patterns;
}

inline bool file_empty_or_missing(const fs::path& path) {
    std::error_code ec;
    return !fs::exists(path, ec) || fs::file_size(path, ec) == 0;
}

inline void append_csv_row(const fs::path& path,
                           const std::string& header,
                           const std::string& row) {
    if (path.empty()) return;
    fs::create_directories(path.parent_path().empty() ? fs::path(".") : path.parent_path());
    const bool write_header = file_empty_or_missing(path);
    std::ofstream out(path, std::ios::app);
    if (!out) throw std::runtime_error("cannot open csv: " + path.string());
    if (write_header) out << header << '\n';
    out << row << '\n';
}

inline std::string csv_quote(const std::string& value) {
    if (value.find_first_of(",\"\n\r") == std::string::npos) return value;
    std::string out = "\"";
    for (char c : value) {
        if (c == '"') out += "\"\"";
        else out += c;
    }
    out += '"';
    return out;
}

inline std::string seconds(double value) {
    std::ostringstream os;
    os << std::fixed << std::setprecision(6) << value;
    return os.str();
}

inline std::uintmax_t sum_files_with_extension(const fs::path& dir,
                                               const std::string& ext) {
    std::uintmax_t total = 0;
    std::error_code ec;
    if (!fs::exists(dir, ec)) return 0;
    for (const auto& entry : fs::recursive_directory_iterator(dir, ec)) {
        if (ec) break;
        if (!entry.is_regular_file(ec)) continue;
        if (!ext.empty() && entry.path().extension() != ext) continue;
        total += entry.file_size(ec);
    }
    return total;
}

inline std::uintmax_t sum_existing_files(const std::vector<fs::path>& paths) {
    std::uintmax_t total = 0;
    std::error_code ec;
    for (const auto& path : paths) {
        if (fs::is_regular_file(path, ec)) total += fs::file_size(path, ec);
    }
    return total;
}

// Tamaño total del índice consultable (grilla + tries de ambas direcciones).
inline std::size_t lz_grid_bytes(const lz77tax::LZ77Index& idx) {
    return idx.index_bytes();
}

inline std::vector<fs::path> lz_files(const fs::path& prefix) {
    return {
        fs::path(prefix).replace_extension(".meta"),
        fs::path(prefix).replace_extension(".grid"),
        fs::path(prefix).replace_extension(".trie"),
        // sub-índice reverso
        fs::path(prefix.string() + "_rev").replace_extension(".meta"),
        fs::path(prefix.string() + "_rev").replace_extension(".grid"),
        fs::path(prefix.string() + "_rev").replace_extension(".trie"),
    };
}

// NOTA: para el tamaño real del sr-index consultable usar
// SrIndexLocator::size_in_bytes() (serializa solo las estructuras cargadas en
// memoria). Sumar los .sdsl del directorio sobreestima ~5× porque incluye
// intermedios de construcción (p.ej. versiones planas con gemelo comprimido).

}  // namespace bench
