#include "index.hpp"
#include "lz77/parser.hpp"
#include "lz77/grid.hpp"

namespace lz77tax {

bool LZ77Index::build(const std::string& reference_file) {
    // TODO: Implement full pipeline:
    // 1. Load reference from FASTA file
    // 2. Parse LZ77 decomposition
    // 3. Build 2D grid
    // 4. Construct Wavelet Tree with RMQ
    // 5. Build bitvectors for BWT↔frase mapping
    return true;
}

bool LZ77Index::save(const std::string& output_file) const {
    // TODO: Serialize index to disk
    return true;
}

bool LZ77Index::load(const std::string& input_file) {
    // TODO: Load index from disk
    return true;
}

}  // namespace lz77tax
