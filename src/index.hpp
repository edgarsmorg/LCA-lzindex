#pragma once

#include "lz77/parser.hpp"
#include "lz77/grid.hpp"

#include <sdsl/int_vector.hpp>

#include <bitset>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>
#include <filesystem>

namespace lz77tax {

class TextOracle;
class MirrorOracle;

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
 * El índice usa dos sparse suffix arrays (PSA) con búsqueda binaria para el
 * matching en ambas direcciones (forward y reverso). El texto original debe
 * estar disponible en query-time: la búsqueda compara los patrones contra él,
 * lo que la hace exacta (sin falsos positivos).
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
     * Guarda el índice a disco con prefijo dado. Escribe solo <prefix>.meta y
     * <prefix>.grid: las PSA no se serializan (psa_fwd = grid.text_pos;
     * psa_rev se reconstruye desde la grilla en load).
     */
    void save(const std::filesystem::path& prefix) const;

    /**
     * Carga el índice desde los archivos generados por save().
     * Requiere el texto original para el matching por búsqueda binaria en queries.
     *
     * @param prefix    Prefijo de archivos (mismo que en save())
     * @param text      Texto original sin centinela (se usa en la búsqueda)
     */
    /// Carga el índice. Un self-index (con `.orac`) se carga sin `text`; para
    /// índices antiguos sin accesor serializado, `text` es obligatorio.
    void load(const std::filesystem::path& prefix, const std::string& text = std::string());

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

    /// Accesor al texto derivado del parsing (evita almacenar T). Experimental:
    /// build_oracle() lo construye; use_oracle(true) hace que las consultas lo
    /// usen en lugar del texto plano.
    template <class TextT>
    size_t locate_leftmost_impl(const TextT& text_s, const std::string& pattern,
                                const uint8_t* fsmp, const uint8_t* bsmp) const;

    template <class TextT>
    size_t count_impl(const TextT& text_s, const std::string& pattern,
                      const uint8_t* fsmp, const uint8_t* bsmp) const;

    void   build_oracle();
    void   build_samples();
    void   use_oracle(bool on);
    size_t oracle_bytes() const;

    /// Bytes del sustrato de búsqueda de esta dirección (las dos PSA).
    /// Conserva el nombre trie_bytes por compatibilidad con las tools.
    size_t trie_bytes() const;

    /// Bytes totales del índice: grilla + tries de esta dirección MÁS el
    /// sub-índice reverso (si existe). Es el tamaño real para calcular bpc.
    size_t index_bytes() const;

private:
    /// Construye el índice sobre `text` (una sola dirección, sin sub-índice reverso).
    void build_core(const std::string& text);

    /// Reconstruye psa_rev_ desde la grilla (no se serializa). Se llama al final
    /// de build_core y en load, tras cargar la grilla.
    void rebuild_psa_rev();

    /// Inicios de frase derivados de la grilla: [0] + sorted(grid.text_pos).
    /// Los usa el accesor al cargar (start_ no se serializa).
    std::vector<size_t> derive_starts() const;

    /// Carga recursiva. En modo self-index (self_index=true) el índice se arma sin
    /// texto: el directo carga su accesor de `.orac`, el reverso recibe un espejo
    /// del accesor del padre (parent_oracle). En modo compatibilidad carga el texto.
    void load_rec(const std::filesystem::path& prefix, const std::string& text,
                  bool self_index, const TextOracle* parent_oracle);

    /// Posición de la ocurrencia primaria más a la izquierda de `pattern`, o
    /// SIZE_MAX si no ocurre. Base de locate_extremal en ambas direcciones.
    size_t locate_leftmost(const std::string& pattern) const;

    Grid2D                     grid_;
    size_t                     n_ = 0;
    size_t                     z_ = 0;
    std::bitset<256>           alphabet_{};
    std::unique_ptr<TextOracle>   oracle_;    ///< accesor sobre T (solo en el índice directo)
    std::unique_ptr<MirrorOracle> mirror_;    ///< vista espejada de oracle_ (solo en el sub-índice reverso)
    std::vector<uint8_t>       fwd_smp_;      ///< muestra por punto de grilla (orden X), kSmp bytes c/u
    std::vector<uint8_t>       bwd_smp_;      ///< muestra por punto de grilla (orden Y), kSmp bytes c/u
    bool                       use_oracle_ = false;
    std::string                text_s_;           ///< texto con centinela (base del matching por búsqueda binaria)
    sdsl::int_vector<>         psa_rev_;          ///< posiciones end_k en orden de sufijo hacia atrás (Y-rank); reconstruida en load, NO serializada
    std::unique_ptr<LZ77Index> rev_index_;     ///< índice sobre el texto reverso (nullptr en el propio sub-índice reverso)
};

}  // namespace lz77tax
