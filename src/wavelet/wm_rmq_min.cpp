#include "wm_rmq_min.hpp"

#include <algorithm>
#include <cassert>
#include <climits>
#include <numeric>

#include <sdsl/bits.hpp>
#include <sdsl/construct.hpp>
#include <sdsl/util.hpp>

namespace lz77tax {

// ─────────────────────────────────────────────────────────────────────────────
// Construcción
// ─────────────────────────────────────────────────────────────────────────────

void WmMinRmq::build(const std::vector<size_t>& y_values,
                     const sdsl::int_vector<>& text_pos,
                     size_t sigma) {
    assert(y_values.size() == text_pos.size());
    n_     = y_values.size();
    sigma_ = sigma;
    rmqs_.clear();
    wm_   = WmAccess{};
    if (n_ == 0 || sigma_ == 0) return;

    // ── 1. Construir wavelet matrix desde y_values ────────────────────────────
    {
        sdsl::int_vector<> R(n_, 0, sdsl::bits::hi(sigma_ - 1) + 1);
        for (size_t i = 0; i < n_; ++i) R[i] = y_values[i];
        sdsl::construct_im(static_cast<sdsl::wm_int<>&>(wm_), R);
    }

    const uint32_t L = wm_.max_level;
    // max_level RMQs para niveles 0..L-1 (dentro del bitvector)
    // +1 RMQ extra para el nivel hoja L (fuera del bitvector, símbolos únicos).
    rmqs_.resize(L + 1);

    // ── 2. Construir un RMQ por nivel, trazando la permutación ────────────────
    // perm[i] = índice original del elemento en la posición i del nivel actual.
    // Nivel 0: orden original → perm = identidad.
    std::vector<size_t> perm(n_);
    std::iota(perm.begin(), perm.end(), 0);

    for (uint32_t k = 0; k <= L; ++k) {
        // Pesos al nivel k en el orden de la permutación actual.
        {
            size_t max_w = 0;
            for (size_t i = 0; i < n_; ++i)
                max_w = std::max(max_w, static_cast<size_t>(text_pos[perm[i]]));

            const uint8_t w = max_w ? static_cast<uint8_t>(sdsl::bits::hi(max_w) + 1) : 1;
            sdsl::int_vector<> weights(n_, 0, w);
            for (size_t i = 0; i < n_; ++i)
                weights[i] = text_pos[perm[i]];

            rmqs_[k] = sdsl::rmq_succinct_sct<true>(&weights);
            // weights se destruye aquí; rmq solo retiene el Cartesian tree en BP.
        }

        if (k == L) break;

        // Actualizar perm para nivel k+1: partición estable por el bit del nivel k.
        // wm_.tree[k*n_ + i] == 0 → va a la izquierda; == 1 → a la derecha.
        std::vector<size_t> left_perm, right_perm;
        left_perm.reserve(n_);
        right_perm.reserve(n_);
        const size_t kn = k * n_;
        for (size_t i = 0; i < n_; ++i) {
            if (!wm_.tree[kn + i])
                left_perm.push_back(perm[i]);
            else
                right_perm.push_back(perm[i]);
        }
        perm.clear();
        perm.insert(perm.end(), left_perm.begin(),  left_perm.end());
        perm.insert(perm.end(), right_perm.begin(), right_perm.end());
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Unwind: posición local j en depth from_depth → índice original (depth 0).
// ─────────────────────────────────────────────────────────────────────────────
//
// La wavelet matrix almacena max_level niveles en m_tree[0..max_level*n-1].
// Al nivel k, posición local j (0-indexed) corresponde a posición global k*n+j.
// Para ir de depth d a depth d-1 (d >= 1):
//   - Nivel parent = d-1, offset en m_tree = parent*n_.
//   - Si j < zero_cnt_at(parent): el elemento vino de un bit=0 en ese nivel.
//     → su posición en el nivel parent = select_0_en_nivel_parent(j + 1).
//   - Si j >= zero_cnt_at(parent): vino de un bit=1; rank_ones = j - zero_cnt.
//     → su posición en el nivel parent = select_1_en_nivel_parent(rank_ones + 1).

size_t WmMinRmq::unwind(size_t j, uint32_t from_depth) const {
    for (uint32_t d = from_depth; d >= 1; --d) {
        const size_t parent = d - 1;
        const size_t pn     = parent * n_;          // offset de nivel parent en m_tree
        if (j < wm_.zero_cnt_at(parent)) {
            // j es el índice (0-based) entre los ceros del nivel parent.
            const size_t zeros_before = pn - wm_.rank_level_at(parent);
            // tree_sel0 es 1-indexed en sdsl.
            const size_t abs = wm_.tree_sel0(zeros_before + j + 1);
            j = abs - pn;
        } else {
            const size_t rank_ones = j - wm_.zero_cnt_at(parent);
            const size_t abs = wm_.tree_sel1(wm_.rank_level_at(parent) + rank_ones + 1);
            j = abs - pn;
        }
    }
    return j;  // índice en el nivel 0 = posición original
}

// ─────────────────────────────────────────────────────────────────────────────
// Consulta recursiva
// ─────────────────────────────────────────────────────────────────────────────
//
// Parámetros:
//   k          nivel actual en la wavelet matrix (0 = raíz, max_level = hoja)
//   cur_lo/hi  rango X al nivel k (indices 0-based dentro del nivel)
//   sym_lo/hi  rango de y-valores que cubre el nodo virtual actual
//   y_lo/hi    rango de y del rectángulo de consulta
//
// Retorna el índice global del argmin (SIZE_MAX si el rectángulo está vacío).
// Actualiza count_acc con el número de elementos encontrados.

size_t WmMinRmq::query_rec(uint32_t k,
                            size_t cur_lo, size_t cur_hi,
                            size_t sym_lo, size_t sym_hi,
                            size_t y_lo,   size_t y_hi,
                            const sdsl::int_vector<>& text_pos,
                            size_t& count_acc) const {
    if (cur_lo > cur_hi)                    return SIZE_MAX;
    if (sym_hi < y_lo || sym_lo > y_hi)    return SIZE_MAX;  // nodo fuera del rango y

    if (y_lo <= sym_lo && sym_hi <= y_hi) {
        // Nodo virtual completamente cubierto: RMQ O(1) + unwind O(log σ).
        // A nivel k, todos los elementos de [cur_lo, cur_hi] tienen y ∈ [sym_lo, sym_hi] ⊆ [y_lo, y_hi].
        const size_t j_local  = rmqs_[k](cur_lo, cur_hi);
        const size_t j_global = unwind(j_local, k);
        count_acc += cur_hi - cur_lo + 1;
        return j_global;
    }

    // Nodo hoja (depth == max_level, sym_lo == sym_hi): si llega aquí es porque
    // la cobertura parcial lo escapó. sym_lo != sym_hi no ocurre en depth max_level,
    // así que esto solo pasa en casos de rango y vacío (ya manejado arriba).
    // Defensa adicional:
    assert(k < wm_.max_level && "hoja sin cobertura total — rango y inválido");

    // Proyectar rango X al siguiente nivel via rank en el bitvector plano.
    const size_t kn     = k * n_;
    const size_t rl_lo  = wm_.tree_rank1(kn + cur_lo)       - wm_.rank_level_at(k);
    const size_t rl_hi1 = wm_.tree_rank1(kn + cur_hi + 1)   - wm_.rank_level_at(k);
    const size_t r0_lo  = cur_lo  - rl_lo;
    const size_t r0_hi1 = (cur_hi + 1) - rl_hi1;
    const size_t zc     = wm_.zero_cnt_at(k);

    // Midpoint: en el nivel k se comprueba el bit de posición (max_level-1-k).
    // El hijo izquierdo cubre [sym_lo, mid], el derecho [mid+1, sym_hi].
    const size_t bit_w = size_t(1) << (wm_.max_level - 1 - k);
    const size_t mid   = sym_lo + bit_w - 1;

    size_t cand_L = SIZE_MAX, cand_R = SIZE_MAX;

    if (y_lo <= mid && r0_lo < r0_hi1)
        cand_L = query_rec(k + 1, r0_lo,      r0_hi1 - 1,
                           sym_lo, mid,        y_lo, y_hi,
                           text_pos, count_acc);

    if (y_hi > mid && rl_lo < rl_hi1)
        cand_R = query_rec(k + 1, zc + rl_lo, zc + rl_hi1 - 1,
                           mid + 1, sym_hi,    y_lo, y_hi,
                           text_pos, count_acc);

    if (cand_L == SIZE_MAX) return cand_R;
    if (cand_R == SIZE_MAX) return cand_L;
    return (text_pos[cand_L] <= text_pos[cand_R]) ? cand_L : cand_R;
}

WmMinRmq::ArgminResult
WmMinRmq::range_argmin_2d(size_t x_lo, size_t x_hi,
                           size_t y_lo, size_t y_hi,
                           const sdsl::int_vector<>& text_pos) const {
    if (n_ == 0 || wm_.empty() || x_lo > x_hi) return {0, SIZE_MAX};
    // El rango de y del nodo raíz abarca todos los valores de max_level bits.
    const size_t sym_hi_root = (size_t(1) << wm_.max_level) - 1;
    size_t count = 0;
    const size_t best_global = query_rec(0, x_lo, x_hi,
                                         0, sym_hi_root,
                                         y_lo, y_hi,
                                         text_pos, count);
    if (count == 0) return {0, SIZE_MAX};
    return {count, best_global};
}

// ─────────────────────────────────────────────────────────────────────────────
// Métricas
// ─────────────────────────────────────────────────────────────────────────────

WmMinRmq::SizeBreakdown WmMinRmq::size_breakdown() const {
    const size_t bv_bytes  = sdsl::size_in_bytes(wm_.tree);
    const size_t wm_total  = sdsl::size_in_bytes(wm_);
    const size_t rs_bytes  = wm_total > bv_bytes ? wm_total - bv_bytes : 0;
    size_t rmq_bytes = 0;
    for (const auto& r : rmqs_) rmq_bytes += sdsl::size_in_bytes(r);
    return {bv_bytes, rs_bytes, rmq_bytes};
}

size_t WmMinRmq::size_in_bytes() const {
    const auto bd = size_breakdown();
    return bd.bv + bd.rank_sel + bd.rmq;
}

}  // namespace lz77tax
