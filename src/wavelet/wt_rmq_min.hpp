#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include <sdsl/bit_vector_il.hpp>
#include <sdsl/rank_support_v.hpp>
#include <sdsl/select_support_mcl.hpp>
#include <sdsl/rmq_succinct_sct.hpp>
#include <sdsl/bits.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>

namespace lz77tax {

/**
 * Wavelet Tree binario con RMQ topológico por nodo (Ferrada-Navarro style).
 *
 * Resuelve argmin 2D en una secuencia de puntos (x, y, val) en O(log² σ):
 *   dado un rectángulo [x_lo, x_hi] × [y_lo, y_hi],
 *   retorna el índice global del punto con valor val mínimo en ese rectángulo.
 *
 * Diseño:
 *   - Cada nodo interno almacena bv (routing), rank1, sel0, sel1 y un RMQ
 *     topológico sobre los valores del nodo (construido desde int_vector temporal
 *     que se descarta tras build — el RMQ solo guarda el Cartesian tree en BP).
 *   - Los valores (text_pos) viven UNA SOLA VEZ en Grid2D. El WT no los retiene.
 *   - Para recuperar el índice global desde un nodo hoja: unwind_to_root via
 *     select_0/select_1 en cadena hacia la raíz, O(log σ) por unwind.
 *   - Consulta: O(log σ) nodos fully-covered × O(log σ) unwind = O(log² σ).
 *   - Espacio: O(z log σ) bits — bitvectors + supports + RMQs, sin duplicar values.
 */
class WtMinRmq {
public:
    struct ArgminResult {
        size_t count;         ///< 0 si rectángulo vacío
        size_t argmin_global; ///< SIZE_MAX si count==0; else índice en [0, n)
    };

    struct SizeBreakdown { size_t bv; size_t rank_sel; size_t rmq; };

    // ── Construcción ─────────────────────────────────────────────────────────
    /**
     * Construye el WT desde la secuencia de coordenadas Y en orden X.
     *
     * @param y_values  y de cada punto en orden X, valores en [0, sigma-1]
     * @param text_pos  valor asociado (text position); usado para construir los
     *                  RMQ por nodo — NO se retiene después de build()
     * @param sigma     cardinalidad del alfabeto de y (= número de puntos, z-1)
     */
    void build(const std::vector<size_t>& y_values,
               const std::vector<size_t>& text_pos,
               size_t sigma);

    // ── Consulta ──────────────────────────────────────────────────────────────
    /**
     * Argmin 2D en rectángulo cerrado [x_lo, x_hi] × [y_lo, y_hi].
     *
     * @param text_pos  array global de valores — Grid2D lo pasa en cada llamada.
     *                  El WT no lo almacena internamente.
     */
    ArgminResult range_argmin_2d(size_t x_lo, size_t x_hi,
                                 size_t y_lo, size_t y_hi,
                                 const std::vector<size_t>& text_pos) const;

    // ── Métricas ──────────────────────────────────────────────────────────────
    size_t size()          const { return n_; }
    size_t alphabet_size() const { return sigma_; }
    size_t size_in_bytes() const;
    SizeBreakdown size_breakdown() const;

    WtMinRmq() = default;
    WtMinRmq(const WtMinRmq&)            = delete;
    WtMinRmq& operator=(const WtMinRmq&) = delete;
    WtMinRmq(WtMinRmq&&)                 = default;
    WtMinRmq& operator=(WtMinRmq&&)      = default;

private:
    struct Node {
        sdsl::bit_vector            bv;
        sdsl::rank_support_v<1>     rank1;
        sdsl::select_support_mcl<0> sel0;
        sdsl::select_support_mcl<1> sel1;
        sdsl::rmq_succinct_sct<true> rmq;   // topológico: sin valores internos
        size_t lo = 0, hi = 0;
        Node*  parent         = nullptr;
        bool   is_right_child = false;
        std::unique_ptr<Node> left, right;
    };

    std::unique_ptr<Node> root_;
    size_t n_     = 0;
    size_t sigma_ = 0;

    static std::unique_ptr<Node> build_rec(std::vector<size_t> ys,
                                           std::vector<size_t> tps,
                                           size_t lo, size_t hi,
                                           Node* parent, bool is_right);

    static size_t unwind_to_root(const Node* node, size_t local_j);

    void query_rec(const Node* node,
                   size_t x_lo, size_t x_hi,
                   size_t y_lo, size_t y_hi,
                   const std::vector<size_t>& tp_global,
                   size_t& count_acc,
                   size_t& best_value,
                   size_t& best_global) const;

    static void accumulate_size(const Node* node,
                                size_t& bv_bytes,
                                size_t& rank_sel_bytes,
                                size_t& rmq_bytes);
};

}  // namespace lz77tax
