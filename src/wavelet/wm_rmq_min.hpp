#pragma once

#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <vector>

#include "wm_shared.hpp"

#include <sdsl/rmq_succinct_sct.hpp>
#include <sdsl/int_vector.hpp>

namespace lz77tax {

// API idéntica a WtRmq<t_min> (wt_rmq_min.hpp) — drop-in replacement.
/**
 * Wavelet Matrix (flat) con RMQ por nivel.
 *
 * Resuelve argmin (t_min=true) o argmax (t_min=false) 2D en O(log² σ).
 * Usa O(log σ) RMQs (uno por nivel) en lugar de O(z) RMQs (uno por nodo).
 * Para z=1.75M reduce de ~63 bpc (WtRmq) a ~2 bpc.
 *
 * Diseño:
 *   - sdsl::wm_int<> como estructura de navegación.
 *   - max_level instancias de sdsl::rmq_succinct_sct<t_min>: una por nivel.
 *   - Consulta: descenso estándar de wavelet matrix; cuando el rango de y
 *     cubre el nodo virtual → RMQ en O(1) + unwind en O(log σ) via select.
 *   - API idéntica para t_min=true/false (drop-in replacement).
 */
template <bool t_min = true>
class WmRmq {
public:
    struct ArgminResult {
        size_t count;         ///< 0 si rectángulo vacío
        size_t argmin_global; ///< SIZE_MAX si count==0; else índice en [0, n)
    };

    struct SizeBreakdown { size_t bv; size_t rank_sel; size_t rmq; };

    // ── Construcción ─────────────────────────────────────────────────────────
    /**
     * @param wm        wavelet matrix compartida, construida sobre los valores y
     * @param text_pos  pesos para la minimización/maximización
     */
    void build(const SharedWm& wm, const sdsl::int_vector<>& text_pos);

    // ── Consulta ──────────────────────────────────────────────────────────────
    ArgminResult range_argmin_2d(const SharedWm& wm,
                                 size_t x_lo, size_t x_hi,
                                 size_t y_lo, size_t y_hi,
                                 const sdsl::int_vector<>& text_pos) const;

    // ── Métricas ──────────────────────────────────────────────────────────────
    size_t size()          const { return n_; }
    size_t alphabet_size() const { return sigma_; }
    size_t size_in_bytes() const;
    SizeBreakdown size_breakdown() const;

    // ── Serialización ─────────────────────────────────────────────────────────
    size_t serialize(std::ostream& out) const;
    void   load(std::istream& in);

    WmRmq() = default;
    WmRmq(const WmRmq&)            = delete;
    WmRmq& operator=(const WmRmq&) = delete;
    WmRmq(WmRmq&&)                 = default;
    WmRmq& operator=(WmRmq&&)      = default;

private:
    std::vector<sdsl::rmq_succinct_sct<t_min>>     rmqs_;  // rmqs_[k]: nivel k

    size_t n_     = 0;
    size_t sigma_ = 0;

    void build_level_rmqs(const SharedWm& wm, const sdsl::int_vector<>& text_pos);

    size_t unwind(const SharedWm& wm, size_t j, uint32_t from_depth) const;

    size_t query_rec(const SharedWm& wm,
                     uint32_t k,
                     size_t cur_lo, size_t cur_hi,
                     size_t sym_lo, size_t sym_hi,
                     size_t y_lo,   size_t y_hi,
                     const sdsl::int_vector<>& text_pos,
                     size_t& count_acc) const;
};

using WmMinRmq = WmRmq<true>;
using WmMaxRmq = WmRmq<false>;

}  // namespace lz77tax
