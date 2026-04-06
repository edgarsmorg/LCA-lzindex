#include "grid.hpp"

namespace lz77tax {

void Grid2D::build(const LZ77Parsing& phrases, size_t text_length) {
    // TODO: Implement grid construction from LZ77 phrases
    // 1. Compute SA of text and SA of reverse text
    // 2. For each phrase, compute coordinates (X, Y)
    // 3. Build sequence R[i] = Y-coordinate where X-coordinate = i
    points_.clear();
}

}  // namespace lz77tax
