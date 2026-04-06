#pragma once

#include "lca.hpp"
#include <vector>
#include <string>
#include <limits>

namespace lz77tax {

/**
 * Clasificador taxonómico: dado un conjunto de posiciones en el texto
 * de referencia, determina el LCA de los genomas donde ocurren.
 *
 * El texto de referencia es la concatenación en orden DFS del árbol filogenético.
 * Cada genoma ocupa un rango [start, end) en el texto.
 */
class Classifier {
public:
    struct GenomeRange {
        size_t start;    ///< Posición inicial en el texto concatenado
        size_t end;      ///< Posición final (exclusiva) en el texto
        int    node_id;  ///< Nodo hoja en el árbol filogenético
    };

    /**
     * Configura el clasificador con el árbol y los rangos de genoma.
     *
     * @param tree     Árbol filogenético construido
     * @param genomes  Rangos de cada genoma en orden DFS (= orden del texto)
     */
    void setup(const PhyloTree& tree, const std::vector<GenomeRange>& genomes) {
        tree_    = &tree;
        genomes_ = genomes;
    }

    /**
     * Clasifica un conjunto de posiciones en el texto.
     *
     * Para cada posición determina a qué genoma pertenece y calcula
     * el LCA de todos los genomas involucrados.
     *
     * @param positions  Posiciones en el texto de referencia
     * @return           ID del nodo LCA, o -1 si no hay posiciones
     */
    int classify(const std::vector<size_t>& positions) const {
        if (positions.empty()) return -1;

        int result = -1;
        for (size_t pos : positions) {
            int genome_node = genome_of(pos);
            if (genome_node < 0) continue;
            result = (result < 0) ? genome_node : tree_->lca(result, genome_node);
        }
        return result;
    }

    /**
     * Versión optimizada: solo necesita la posición mínima y máxima.
     *
     * Válido porque el texto sigue el orden DFS: los genomas con
     * posiciones extremales determinan el LCA de todos los intermedios.
     */
    int classify_extremal(size_t pos_min, size_t pos_max) const {
        int node_min = genome_of(pos_min);
        int node_max = genome_of(pos_max);
        if (node_min < 0 || node_max < 0) return -1;
        return tree_->lca(node_min, node_max);
    }

    /** Retorna el node_id del genoma al que pertenece una posición. */
    int genome_of(size_t pos) const {
        for (const auto& g : genomes_) {
            if (pos >= g.start && pos < g.end) return g.node_id;
        }
        return -1;  // posición en separador
    }

private:
    const PhyloTree*        tree_    = nullptr;
    std::vector<GenomeRange> genomes_;
};

}  // namespace lz77tax
