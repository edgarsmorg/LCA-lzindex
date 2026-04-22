#pragma once

#include "phrase.hpp"

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include <sdsl/wt_int.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/suffix_arrays.hpp>

namespace lz77tax {

/**
 * Grilla 2D de ocurrencias primarias del LZ-index.
 *
 * Cada uno de los z-1 puntos (X_k, Y_k) representa el boundary entre la frase k
 * y la frase k+1:
 *   X_k = ISA[ start_{k+1} ]       rank en SA forward  de T[start_{k+1}..n-1]
 *   Y_k = ISA_rev[ n-2 - end_k ]   rank en rev-SA de T^R[n-2-end_k..n-1]
 *   donde end_k = start_{k+1} - 1 = posición del next_char de la frase k
 *
 * Internamente se almacena como:
 *   wt_    : wt_int<> sobre R[0..z-2] donde R[rank_fwd(X_k)] = rank_rev(Y_k)
 *   bv_fwd_: sd_vector con 1s en las posiciones X_k del SA global
 *   bv_rev_: sd_vector con 1s en las posiciones Y_k del rev-SA global
 *
 * Los bitvectors convierten rangos BWT globales [sp,ep] a rangos relativos
 * [rank(sp), rank(ep+1)-1] sobre el wt_, sin materializar SA ni ISA completos
 * en la consulta.
 *
 * ESCALA: la interfaz de producción recibe csa_wt (FM-index ya construido por
 * ropebwt3) y resuelve cada ISA[i] en O(t_psi) sin materializar el array entero.
 * La interfaz de test acepta el texto crudo y construye el SA en RAM (solo para
 * textos pequeños, <= ~500 MB).
 */
class Grid2D {
public:
    // ── Interfaz de producción (GB+) ─────────────────────────────────────────
    /**
     * Construye la grilla usando FM-indexes ya construidos.
     * Usa csa.isa[i] que evita materializar el ISA completo — seguro a escala GB+.
     *
     * @param phrases   Parsing LZ77 del texto (z frases)
     * @param csa_fwd   FM-index de T (SA de T forward)
     * @param csa_rev   FM-index de T^R (SA del texto reverso)
     */
    void build(const LZ77Parsing& phrases,
               const sdsl::csa_wt<>& csa_fwd,
               const sdsl::csa_wt<>& csa_rev);

    // ── Interfaz de test / small-scale (≤ ~500 MB) ───────────────────────────
    /**
     * Construye la grilla desde el texto crudo.
     * Construye SA + ISA + SA_rev + ISA_rev en RAM (~5x texto) — solo para tests.
     *
     * @param phrases     Parsing LZ77 del texto
     * @param text        Texto original (debe terminar con el carácter más pequeño,
     *                    e.g. '\0' o '$', para que el SA sea correcto)
     */
    void build(const LZ77Parsing& phrases, const std::string& text);

    // ── Consultas ─────────────────────────────────────────────────────────────
    /**
     * Consulta rectangular en la grilla usando rangos BWT globales.
     *
     * Para el split P[0..i-1] | P[i..m-1]:
     *   [sp_right, ep_right] = backward_search(csa_fwd,  P[i..m-1])
     *   [sp_left,  ep_left]  = backward_search(csa_rev,  reverse(P[0..i-1]))
     *
     * Retorna {count, vector<{wt_idx, y_rel}>}.
     * sdsl::range_search_2d (wt_int.hpp línea 717: emplace_back(i-1, path)) devuelve
     * pares {position, value}:
     *   position = wt_idx = rank_fwd_(X_k)  (índice original en el array R)
     *   value    = y_rel  = rank_rev_(Y_k)   (valor almacenado en R[wt_idx])
     */
    using RangeResult = std::pair<
        size_t,
        std::vector<std::pair<sdsl::wt_int<>::value_type,
                              sdsl::wt_int<>::size_type>>>;

    RangeResult query(size_t sp_right, size_t ep_right,
                      size_t sp_left,  size_t ep_left) const;

    /**
     * Consulta extremal: igual que query() pero devuelve solo el mínimo y máximo
     * de start_{k+1} (posición del boundary en el texto) entre los puntos hallados.
     *
     * Complejidad: O(occ · log z) — enumera puntos y hace lookup en text_pos_.
     * Para clasificación taxonómica: la posición real de la ocurrencia del patrón
     * es boundary - split_i, que el llamador (LZ77Index) computa externamente.
     */
    struct ExtremalResult {
        size_t count;        ///< Número de puntos encontrados en el rectángulo
        size_t boundary_min; ///< Mínimo start_{k+1} entre los puntos; SIZE_MAX si count=0
        size_t boundary_max; ///< Máximo start_{k+1} entre los puntos; 0 si count=0
    };

    ExtremalResult query_extremal(size_t sp_right, size_t ep_right,
                                  size_t sp_left,  size_t ep_left) const;

    // ── Accesores ─────────────────────────────────────────────────────────────
    /// Número de puntos en la grilla (= z-1 para z frases)
    size_t point_count() const { return wt_.size(); }

    const sdsl::wt_int<>&    wt()      const { return wt_; }
    const sdsl::sd_vector<>& bv_fwd()  const { return bv_fwd_; }
    const sdsl::sd_vector<>& bv_rev()  const { return bv_rev_; }
    /// Posición en el texto del boundary k+1 para el punto con índice WT wt_idx.
    size_t text_pos(size_t wt_idx) const { return text_pos_[wt_idx]; }

    // rank_fwd_/rank_rev_ guardan punteros raw a bv_fwd_/bv_rev_: no es seguro
    // copiar ni mover Grid2D con los constructores implícitos. Prohibimos copia;
    // si en el futuro se necesita mover, implementar con sdsl::util::init_support.
    Grid2D() = default;
    Grid2D(const Grid2D&) = delete;
    Grid2D& operator=(const Grid2D&) = delete;

private:
    sdsl::wt_int<>                  wt_;
    sdsl::sd_vector<>               bv_fwd_;
    sdsl::sd_vector<>               bv_rev_;
    sdsl::sd_vector<>::rank_1_type  rank_fwd_;
    sdsl::sd_vector<>::rank_1_type  rank_rev_;
    std::vector<size_t>             text_pos_; ///< text_pos_[j] = start_{k+1} del punto con índice WT j

    /// Núcleo compartido: recibe las z-1 coordenadas y los boundaries start_{k+1},
    /// construye bv_fwd_, bv_rev_, wt_ y text_pos_.
    void build_from_coords(const std::vector<std::pair<size_t, size_t>>& coords,
                           const std::vector<size_t>& boundaries,
                           size_t n);
};

}  // namespace lz77tax
