#pragma once

#include "phrase.hpp"
#include <string>

namespace lz77tax {

/**
 * Parser LZ77 desde un Suffix Array.
 *
 * Implementa parsing LZ77 greedy left-to-right.
 * Requiere: SA, ISA, LCP arrays (construidos por sdsl-lite/libdivsufsort).
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
