#pragma once

#include <cstddef>
#include <vector>

namespace lz77tax {

/**
 * Representa un Maximal Exact Match (MEM) entre una query y la referencia.
 */
struct MEM {
    size_t query_start;    ///< Posición en el query donde inicia el match
    size_t ref_start;      ///< Posición en la referencia donde inicia el match
    size_t length;         ///< Largo del match
};

/// Conjunto de MEMs encontrados en una query
using MEMList = std::vector<MEM>;

}  // namespace lz77tax
