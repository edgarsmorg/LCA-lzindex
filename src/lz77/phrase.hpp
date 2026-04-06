#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace lz77tax {

/**
 * Representa una frase en el parsing LZ77.
 *
 * Una frase es un factor de la descomposición LZ77 del texto.
 * Almacena solo información de ocurrencias primarias (no propagamos secundarias).
 */
struct Phrase {
    size_t start_pos;      ///< Posición en T donde comienza esta frase
    size_t length;         ///< Largo de la frase (0 para literal puro)
    uint8_t next_char;     ///< Carácter explícito al final (separador)

    // NOTA: NO almacenamos source (no necesario para nuestro caso de uso)
};

/// Secuencia de frases que forma el parsing LZ77 completo
using LZ77Parsing = std::vector<Phrase>;

}  // namespace lz77tax
