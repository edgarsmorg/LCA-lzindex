#pragma once

#include "lz77/parser.hpp"
#include "lz77/grid.hpp"

#include <bitset>
#include <cstddef>
#include <memory>
#include <string>
#include <filesystem>

namespace lz77tax {

/**
 * Índice LZ77 para ocurrencias primarias.
 *
 * count(P) retorna el número de pares (ocurrencia p, boundary k) tales que
 * P ocurre en p y el boundary k cae dentro de P (p < boundary_k <= p+m-1).
 * Esto es estrictamente menor o igual al total de ocurrencias de P.
 *
 * Para clasificación taxonómica, las posiciones extremales (min/max) de estas
 * ocurrencias primarias son suficientes para determinar el LCA correcto.
 *
 * El índice usa tries Patricia DFUDS (Kreft–Navarro) para búsqueda en ambas
 * direcciones. El texto original debe estar disponible en query-time para
 * verificar los saltos Patricia y evitar falsos positivos.
 */
class LZ77Index {
public:
    LZ77Index();
    ~LZ77Index();

    // ── Interfaz small-scale / tests (≤ ~500 MB) ─────────────────────────────
    /**
     * Construye el índice completo desde el texto crudo.
     * Internamente añade '\0' como centinela, parsea, construye grid y tries.
     * El texto se almacena internamente para verificación Patricia en queries.
     *
     * @param text  Texto sin centinela (e.g. ADN puro, sin '\0')
     */
    void build(const std::string& text);


    // ── Persistencia ──────────────────────────────────────────────────────────
    /**
     * Guarda el índice a disco en varios archivos con prefijo dado.
     * Escribe: <prefix>.meta, <prefix>.grid, <prefix>.trie
     */
    void save(const std::filesystem::path& prefix) const;

    /**
     * Carga el índice desde los archivos generados por save().
     * Requiere el texto original para verificación Patricia en queries.
     *
     * @param prefix    Prefijo de archivos (mismo que en save())
     * @param text      Texto original sin centinela (se usa para verificación)
     */
    void load(const std::filesystem::path& prefix, const std::string& text);

    // ── Query ─────────────────────────────────────────────────────────────────
    /**
     * Cuenta ocurrencias primarias de pattern.
     */
    size_t count(const std::string& pattern) const;

    /**
     * Localiza la posición mínima y máxima entre TODAS las ocurrencias del patrón.
     *
     * La mínima (más a la izquierda) es siempre primaria y se obtiene del índice
     * directo. La máxima (más a la derecha) NO es en general el máximo de las
     * primarias del texto directo (puede ser secundaria); se obtiene como la
     * ocurrencia más a la izquierda de P^R en el índice sobre el texto reverso
     * (que sí es primaria allí) y se remapea a coordenadas del texto directo.
     *
     * Retorna {SIZE_MAX, 0} si el patrón no ocurre.
     */
    std::pair<size_t, size_t> locate_extremal(const std::string& pattern) const;

    // ── Accesores ─────────────────────────────────────────────────────────────
    /// Tamaño del texto indexado (incluyendo centinela '\0')
    size_t text_size()    const { return n_; }

    /// Número de frases LZ77 (cacheado; phrases_ se libera tras build)
    size_t phrase_count() const { return z_; }

    /// Número de puntos en la grilla (= z-1 para z frases)
    size_t grid_points()  const { return grid_.point_count(); }

    const Grid2D& grid() const { return grid_; }

    /// Bytes de las estructuras de tries (SST + rev_trie + DAC skips)
    size_t trie_bytes() const;

    /// Bytes totales del índice: grilla + tries de esta dirección MÁS el
    /// sub-índice reverso (si existe). Es el tamaño real para calcular bpc.
    size_t index_bytes() const;

private:
    struct TrieImpl;               // PIMPL: dfuds* sst + rev_trie + FTRep* skips

    /// Construye el índice sobre `text` (una sola dirección, sin sub-índice reverso).
    void build_core(const std::string& text);

    /// Posición de la ocurrencia primaria más a la izquierda de `pattern`, o
    /// SIZE_MAX si no ocurre. Base de locate_extremal en ambas direcciones.
    size_t locate_leftmost(const std::string& pattern) const;

    Grid2D                     grid_;
    size_t                     n_ = 0;
    size_t                     z_ = 0;
    std::bitset<256>           alphabet_{};
    std::string                text_s_;           ///< texto con centinela (para verificación Patricia)
    mutable std::unique_ptr<TrieImpl>  trie_;  // mutable: dfuds no tiene metodos const
    std::unique_ptr<LZ77Index> rev_index_;     ///< índice sobre el texto reverso (nullptr en el propio sub-índice reverso)
};

}  // namespace lz77tax
