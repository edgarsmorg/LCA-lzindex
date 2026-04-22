#pragma once

#include "phrase.hpp"
#include <string>

namespace lz77tax {

/**
 * Parser LZ77 greedy left-to-right desde Suffix Array.
 *
 * Construye SA con divsufsort (O(n log n)) y LCP con Kasai (O(n)).
 * El greedy parsing itera sobre las z frases resultantes.
 *
 * Límite de tamaño: usa int32_t internamente → texto máximo 2^31-1 ≈ 2 GB.
 * Para los corpus objetivo (E. coli 108 MB, cere 461 MB) es suficiente.
 */
class LZ77Parser {
public:
    /**
     * Parsea el texto y retorna la secuencia de frases LZ77.
     *
     * @param text El texto a parsear
     * @return Vector de frases en orden
     */
    LZ77Parsing parse(const std::string& text);

    /// Cantidad de frases en el último parsing
    size_t phrase_count() const { return phrase_count_; }

private:
    size_t phrase_count_ = 0;
};

}  // namespace lz77tax
