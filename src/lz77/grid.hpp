#pragma once

#include "phrase.hpp"
#include "../wavelet/wm_shared.hpp"
#include "../wavelet/wm_rmq_min.hpp"

#include <cstddef>
#include <iosfwd>
#include <utility>
#include <vector>

#include <sdsl/int_vector.hpp>

namespace lz77tax {

/**
 * Grilla 2D de ocurrencias primarias del LZ-index.
 *
 * Cada uno de los z-1 puntos (X_k, Y_k) representa el boundary entre la frase k
 * y la frase k+1:
 *   X_k = rank del sufijo T[start_{k+1}..] entre los inicios de frase
 *   Y_k = rank del prefijo invertido que termina en end_k = start_{k+1} - 1
 *
 * Las coordenadas llegan ya relativas (0-indexed sobre los z-1 puntos): las
 * calcula el llamador con los arreglos de sufijos dispersos. La grilla nunca ve
 * rangos BWT globales ni materializa SA/ISA.
 *
 * Internamente se almacena como:
 *   wm_              : wm_int<> sobre R[0..z-2], donde R[X_k] = Y_k
 *   text_pos_        : start_{k+1} de cada punto, en orden X
 *   phrase_total_len_: largo total de la frase k de cada punto, en orden X
 *   wm_min_rmq_      : RMQ de minimo por nivel sobre text_pos_
 */
class Grid2D {
public:
    /**
     * Construye la grilla desde coordenadas ya relativas.
     *
     * @param R           R[i] = coordenada Y (0-indexed) del punto con X = i
     * @param boundaries  boundaries[i] = start_{k+1} del punto con X = i
     * @param phrase_lens phrase_lens[i] = largo total de la frase k de ese punto
     */
    void build(const std::vector<size_t>& R,
               const std::vector<size_t>& boundaries,
               const std::vector<size_t>& phrase_lens);

    // ── Consultas (rangos relativos, 0-indexed inclusivos) ────────────────────

    /**
     * Consulta rectangular [x_lo,x_hi] x [y_lo,y_hi].
     *
     * Retorna {count, vector<{x, y}>}: sdsl::range_search_2d devuelve pares
     * {position, value} = {X_k, Y_k}.
     */
    using RangeResult = std::pair<
        size_t,
        std::vector<std::pair<sdsl::wm_int<>::value_type,
                              sdsl::wm_int<>::size_type>>>;

    RangeResult query(size_t x_lo, size_t x_hi,
                      size_t y_lo, size_t y_hi) const;

    /**
     * Consulta extremal minima (O(log^2 z), sin enumerar puntos).
     *
     * Retorna {count, boundary_min} donde boundary_min = start_{k+1} minimo
     * entre los puntos del rectangulo. Si count==0, boundary_min = SIZE_MAX.
     */
    struct MinResult {
        size_t count;        ///< 0 si rectangulo vacio
        size_t boundary_min; ///< minimo start_{k+1}; SIZE_MAX si count==0
        size_t wt_idx;       ///< indice del punto extremal; SIZE_MAX si count==0
    };

    MinResult query_min(size_t x_lo, size_t x_hi,
                        size_t y_lo, size_t y_hi) const;

    /// Como query_min, pero solo considera puntos con phrase_total_len >= min_phrase_len.
    MinResult query_min_filtered(size_t x_lo, size_t x_hi,
                                 size_t y_lo, size_t y_hi,
                                 size_t min_phrase_len) const;

    /**
     * Consulta especial para patrones end-aligned.
     * Retorna ocurrencias donde el patron termina exactamente al final de la frase k,
     * siempre y cuando quepa en la frase (phrase_total_len_ >= plen).
     */
    struct SpecialResult {
        size_t count;
        size_t occ_min_pos;
        size_t occ_max_pos;
    };

    SpecialResult query_special(size_t y_lo, size_t y_hi, size_t plen) const;

    // Accesores
    /// Numero de puntos en la grilla (= z-1 para z frases)
    size_t point_count() const { return wm_.size(); }

    const SharedWm& wm()          const { return wm_; }
    const WmMinRmq& wm_min_rmq()  const { return wm_min_rmq_; }
    /// Posicion en el texto del boundary k+1 para el punto con coordenada X wt_idx.
    size_t text_pos(size_t wt_idx) const { return text_pos_[wt_idx]; }
    /// Total de caracteres de la frase k para el punto con coordenada X wt_idx (= length + 1).
    size_t phrase_total_len(size_t wt_idx) const { return phrase_total_len_[wt_idx]; }

    struct SizeBreakdown {
        size_t wm;
        size_t text_pos;
        size_t phrase_total_len;
        size_t wm_min_rmq;
    };
    SizeBreakdown size_breakdown() const;

    // Serializacion
    size_t serialize(std::ostream& out) const;
    void   load(std::istream& in);

    /// Reconstruye phrase_total_len_ desde text_pos_ (no se serializa). Público
    /// para poder reconstruir tras build() en el path de construcción directa.
    void   reconstruct_phrase_total_len();

    // wm_min_rmq_ no es copiable (contiene los RMQ por nivel); Grid2D hereda esa
    // restriccion. Mover si es necesario.
    Grid2D() = default;
    Grid2D(const Grid2D&) = delete;
    Grid2D& operator=(const Grid2D&) = delete;

private:
    SharedWm                        wm_;
    sdsl::int_vector<>              text_pos_;        ///< text_pos_[j] = start_{k+1}, compacto
    sdsl::int_vector<>              phrase_total_len_;///< phrase_total_len_[j] = total_len de la frase k
    WmMinRmq                        wm_min_rmq_;
};

}  // namespace lz77tax
