#pragma once

#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <vector>

#include <sdsl/wm_int.hpp>
#include <sdsl/rmq_succinct_sct.hpp>
#include <sdsl/int_vector.hpp>

namespace lz77tax {

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
     * @param y_values  y de cada punto en orden X, valores en [0, sigma-1]
     * @param text_pos  pesos para la minimización/maximización
     * @param sigma     cardinalidad del alfabeto de y (= número de puntos)
     */
    void build(const std::vector<size_t>& y_values,
               const sdsl::int_vector<>& text_pos,
               size_t sigma);

    // ── Consulta ──────────────────────────────────────────────────────────────
    ArgminResult range_argmin_2d(size_t x_lo, size_t x_hi,
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
    struct WmAccess : public sdsl::wm_int<> {
        size_t zero_cnt_at(size_t k)    const { return m_zero_cnt[k]; }
        size_t rank_level_at(size_t k)  const { return m_rank_level[k]; }
        size_t tree_rank1(size_t i)     const { return m_tree_rank(i); }
        size_t tree_sel0(size_t i)      const { return m_tree_select0(i); }
        size_t tree_sel1(size_t i)      const { return m_tree_select1(i); }
    };

    WmAccess                                       wm_;
    std::vector<sdsl::rmq_succinct_sct<t_min>>     rmqs_;  // rmqs_[k]: nivel k

    size_t n_     = 0;
    size_t sigma_ = 0;

    size_t unwind(size_t j, uint32_t from_depth) const;

    size_t query_rec(uint32_t k,
                     size_t cur_lo, size_t cur_hi,
                     size_t sym_lo, size_t sym_hi,
                     size_t y_lo,   size_t y_hi,
                     const sdsl::int_vector<>& text_pos,
                     size_t& count_acc) const;
};

using WmMinRmq = WmRmq<true>;
using WmMaxRmq = WmRmq<false>;

}  // namespace lz77tax
