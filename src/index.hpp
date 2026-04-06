#pragma once

#include <string>
#include <vector>
#include <memory>
#include "lz77/phrase.hpp"
#include "mem/mem.hpp"

namespace lz77tax {

/**
 * Índice LZ77 comprimido para clasificación taxonómica de lecturas de ADN.
 *
 * Implementa:
 * 1. Parsing LZ77 del texto de referencia
 * 2. Grilla 2D de puntos primarios con Wavelet Tree
 * 3. RMQ sucinto integrado en el Wavelet Tree
 * 4. Búsqueda de MEMs y mapeo a genomas
 * 5. Cómputo de LCA en el árbol filogenético
 */
class LZ77Index {
public:
    LZ77Index() = default;
    ~LZ77Index() = default;

    /**
     * Construye el índice a partir de un archivo de referencia.
     *
     * @param reference_file Ruta al archivo FASTA con la concatenación de genomas
     * @return true si la construcción fue exitosa
     */
    bool build(const std::string& reference_file);

    /**
     * Serializa el índice a disco.
     *
     * @param output_file Ruta donde guardar el índice
     * @return true si la serialización fue exitosa
     */
    bool save(const std::string& output_file) const;

    /**
     * Carga el índice desde disco.
     *
     * @param input_file Ruta del índice guardado
     * @return true si la carga fue exitosa
     */
    bool load(const std::string& input_file);

    /// Cantidad de frases en el parsing LZ77
    size_t phrase_count() const { return phrases_.size(); }

    /// Largo del texto de referencia original
    size_t text_length() const { return text_length_; }

private:
    LZ77Parsing phrases_;
    size_t text_length_ = 0;
    // TODO: agregar Wavelet Tree, RMQ, bitvectors
};

}  // namespace lz77tax
