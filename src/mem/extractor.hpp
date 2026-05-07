#pragma once

#include "mem.hpp"
#include <sdsl/suffix_arrays.hpp>
#include <string>

namespace lz77tax {

/**
 * Extractor de MEMs (Maximal Exact Matches) usando backward_search sobre el FM-index.
 *
 * Para cada endpoint derecho j de la query, determina el menor i tal que Q[i..j]
 * ocurre en T (left-maximal por construcción), luego filtra por right-maximality.
 *
 * Complejidad: O(m² · t_SA) donde m=|query|, t_SA = O(log n). Aceptable para reads ≤ 1000 bp.
 */
class MEMExtractor {
public:
    explicit MEMExtractor(const sdsl::csa_wt<>& csa_fwd) : csa_(csa_fwd) {}

    /**
     * Extrae todas las MEMs de length ≥ min_length entre query y la referencia indexada.
     *
     * @param query       Secuencia query (lectura de ADN)
     * @param min_length  Largo mínimo de MEM a reportar (default: 31)
     * @return            Lista de MEMs: cada una con query_start, ref_start (posición en T),
     *                    y length. ref_start es una ocurrencia representativa, no necesariamente
     *                    la mínima.
     */
    MEMList extract(const std::string& query, size_t min_length = 31) const;

private:
    const sdsl::csa_wt<>& csa_;
};

}  // namespace lz77tax
