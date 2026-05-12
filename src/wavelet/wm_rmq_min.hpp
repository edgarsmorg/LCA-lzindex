#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include <sdsl/wm_int.hpp>
#include <sdsl/rmq_succinct_sct.hpp>
#include <sdsl/int_vector.hpp>

namespace lz77tax {

/**
 * Wavelet Matrix (flat) con RMQ por nivel.
 *
 * Resuelve argmin 2D en O(log² σ) igual que WtMinRmq, pero con espacio
 * radicalmente menor: usa O(log σ) RMQs (uno por nivel) en lugar de O(z)
 * RMQs (uno por nodo). Para z=1.75M reduce de ~63 bpc a ~2 bpc.
 *
 * Diseño:
 *   - sdsl::wm_int<> como estructura de navegación (bitvector plano único,
 *     ONE rank + ONE select_0 + ONE select_1 TOTAL, vs z estructuras en WtMinRmq).
 *   - max_level RMQs succintos: uno por nivel, cada uno sobre n elementos.
 *     Construido sobre los pesos en orden nivel-k (permutación tracking en build).
 *   - Consulta: descenso estándar de wavelet matrix; cuando el rango de y
 *     cubre el nodo virtual → RMQ en O(1) + unwind en O(log σ) via select.
 *   - API idéntica a WtMinRmq para drop-in replacement.
 */
class WmMinRmq {
public:
    struct ArgminResult {
        size_t count;         ///< 0 si rectángulo vacío
        size_t argmin_global; ///< SIZE_MAX si count==0; else índice en [0, n)
    };

    struct SizeBreakdown { size_t bv; size_t rank_sel; size_t rmq; };

    // ── Construcción ─────────────────────────────────────────────────────────
    /**
     * Misma signatura que WtMinRmq::build.
     * @param y_values  y de cada punto en orden X, valores en [0, sigma-1]
     * @param text_pos  pesos para la minimización (text positions)
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

    WmMinRmq() = default;
    WmMinRmq(const WmMinRmq&)            = delete;
    WmMinRmq& operator=(const WmMinRmq&) = delete;
    WmMinRmq(WmMinRmq&&)                 = default;
    WmMinRmq& operator=(WmMinRmq&&)      = default;

private:
    // Subclase thin para exponer miembros protegidos de wm_int sin duplicarlos.
    struct WmAccess : public sdsl::wm_int<> {
        size_t zero_cnt_at(size_t k)    const { return m_zero_cnt[k]; }
        size_t rank_level_at(size_t k)  const { return m_rank_level[k]; }
        size_t tree_rank1(size_t i)     const { return m_tree_rank(i); }
        size_t tree_sel0(size_t i)      const { return m_tree_select0(i); }
        size_t tree_sel1(size_t i)      const { return m_tree_select1(i); }
    };

    WmAccess                                  wm_;
    std::vector<sdsl::rmq_succinct_sct<true>> rmqs_;  // rmqs_[k]: nivel k

    size_t n_     = 0;
    size_t sigma_ = 0;

    // Unwind: posición local j al nivel from_depth → índice original (nivel 0).
    size_t unwind(size_t j, uint32_t from_depth) const;

    // Recursión principal de la consulta 2D.
    // Retorna argmin_global (SIZE_MAX si vacío); actualiza count_acc.
    size_t query_rec(uint32_t k,
                     size_t cur_lo, size_t cur_hi,
                     size_t sym_lo, size_t sym_hi,
                     size_t y_lo,   size_t y_hi,
                     const sdsl::int_vector<>& text_pos,
                     size_t& count_acc) const;
};

}  // namespace lz77tax
