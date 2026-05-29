#pragma once

#include "lz77/parser.hpp"
#include "lz77/grid.hpp"

#include <bitset>
#include <cstddef>
#include <memory>
#include <string>

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
 */
class LZ77Index {
public:
    LZ77Index();
    ~LZ77Index();

    // ── Interfaz small-scale / tests (≤ ~500 MB) ─────────────────────────────
    /**
     * Construye el índice completo desde el texto crudo.
     * Internamente añade '\0' como centinela, parsea, construye grid y CSAs.
     *
     * @param text  Texto sin centinela (e.g. ADN puro, sin '\0')
     */
    void build(const std::string& text,
               RmqVariant variant = RmqVariant::Wm);


    // ── Query ─────────────────────────────────────────────────────────────────
    /**
     * Cuenta ocurrencias primarias de pattern.
     *
     * Complejidad: O(m² · t_ψ + occ · log z)
     * donde m=|pattern|, t_ψ = intervalo de muestreo del CSA, occ = resultado.
     *
     * Retorna 0 si pattern está vacío, tiene longitud 1, o el índice está vacío.
     */
    size_t count(const std::string& pattern) const;

    /**
     * Localiza la posición mínima y máxima entre todas las ocurrencias primarias
     * de pattern (aquellas que cruzan al menos un boundary de frase).
     *
     * Retorna {pos_min, pos_max} como posiciones de inicio del patrón en el texto.
     * Retorna {SIZE_MAX, 0} si no hay ocurrencias primarias (pattern no existe,
     * longitud < 2, o todas las ocurrencias caen dentro de una sola frase).
     *
     * Uso para clasificación taxonómica:
     *   auto [lo, hi] = idx.locate_extremal(pattern);
     *   int lca_node  = classifier.classify_extremal(lo, hi);
     *
     * Complejidad: O(m² · t_ψ + m · log² z) — usa query_min_2d + query_max_2d,
     * sin enumerar puntos. Más eficiente que la versión anterior O(occ_prim · log z).
     */
    std::pair<size_t, size_t> locate_extremal(const std::string& pattern) const;

    /**
     * Posición mínima entre ocurrencias primarias de pattern.
     *
     * Complejidad: O(m² · t_ψ + m · log² z) — no enumera puntos.
     * Retorna SIZE_MAX si no hay ocurrencias primarias.
     */
    size_t locate_min(const std::string& pattern) const;

    /**
     * Posición máxima entre ocurrencias primarias de pattern.
     *
     * Simétrico a locate_min() usando WtMinRmq sobre valores de posición invertidos.
     * Complejidad: O(m² · t_ψ + m · log² z) — no enumera puntos.
     * Retorna SIZE_MAX si no hay ocurrencias primarias.
     */
    size_t locate_max(const std::string& pattern) const;

    // ── Accesores ─────────────────────────────────────────────────────────────
    /// Tamaño del texto indexado (incluyendo centinela '\0')
    size_t text_size()    const { return n_; }

    /// Número de frases LZ77 (cacheado; phrases_ se libera tras build)
    size_t phrase_count() const { return z_; }

    /// Número de puntos en la grilla (= z-1 para z frases)
    size_t grid_points()  const { return grid_.point_count(); }

    const Grid2D&      grid()    const { return grid_; }
    size_t csa_fwd_bytes()       const;
    size_t csa_rev_bytes()       const;

private:
    struct RCSAImpl;              // PIMPL: definido en index.cpp (r_csa.h solo allí)

    Grid2D                     grid_;
    size_t                     n_ = 0;
    size_t                     z_ = 0;       ///< número de frases LZ77 (cacheado post-build)
    std::bitset<256>           alphabet_{};  ///< chars presentes en el texto indexado
    std::unique_ptr<RCSAImpl>  rcsa_;
};

}  // namespace lz77tax
