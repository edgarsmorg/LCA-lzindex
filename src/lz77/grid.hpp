#pragma once

#include "phrase.hpp"
#include <vector>
#include <cstddef>

namespace lz77tax {

/**
 * Grilla 2D bidimensional que codifica cruces de límites de frase.
 *
 * Cada punto (x, y) en la grilla representa un cruce de límite de frase:
 *   - x: coordenada en el SA reverso (rango final de la frase)
 *   - y: coordenada en el SA directo (rango inicial de la frase)
 *
 * Se almacena como secuencia R[1..z] donde R[i] = coordenada Y del punto con coord X = i.
 */
class Grid2D {
public:
    /**
     * Construye la grilla a partir de las frases LZ77.
     *
     * @param phrases Secuencia de frases del parsing LZ77
     * @param text_length Largo del texto original
     */
    void build(const LZ77Parsing& phrases, size_t text_length);

    /// Cantidad de puntos en la grilla (= cantidad de frases z)
    size_t point_count() const { return points_.size(); }

    /// Secuencia R[i] = coordenada Y del punto con coord X = i
    const std::vector<size_t>& r_sequence() const { return points_; }

private:
    std::vector<size_t> points_;  ///< Secuencia R de valores Y
};

}  // namespace lz77tax
