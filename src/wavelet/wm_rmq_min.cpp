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

template <bool t_min>
void WmRmq<t_min>::build(const std::vector<size_t>& y_values,
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
    rmqs_.resize(L + 1);

    // ── 2. Construir un RMQ por nivel, trazando la permutación ────────────────
    std::vector<size_t> perm(n_);
    std::iota(perm.begin(), perm.end(), 0);

    for (uint32_t k = 0; k <= L; ++k) {
        {
            size_t max_w = 0;
            for (size_t i = 0; i < n_; ++i)
                max_w = std::max(max_w, static_cast<size_t>(text_pos[perm[i]]));

            const uint8_t w = max_w ? static_cast<uint8_t>(sdsl::bits::hi(max_w) + 1) : 1;
            sdsl::int_vector<> weights(n_, 0, w);
            for (size_t i = 0; i < n_; ++i)
                weights[i] = text_pos[perm[i]];

            rmqs_[k] = sdsl::rmq_succinct_sct<t_min>(&weights);
        }

        if (k == L) break;

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
// Unwind: posición local j al nivel from_depth → índice original (nivel 0).
// ─────────────────────────────────────────────────────────────────────────────

template <bool t_min>
size_t WmRmq<t_min>::unwind(size_t j, uint32_t from_depth) const {
    for (uint32_t d = from_depth; d >= 1; --d) {
        const size_t parent = d - 1;
        const size_t pn     = parent * n_;
        if (j < wm_.zero_cnt_at(parent)) {
            const size_t zeros_before = pn - wm_.rank_level_at(parent);
            const size_t abs = wm_.tree_sel0(zeros_before + j + 1);
            j = abs - pn;
        } else {
            const size_t rank_ones = j - wm_.zero_cnt_at(parent);
            const size_t abs = wm_.tree_sel1(wm_.rank_level_at(parent) + rank_ones + 1);
            j = abs - pn;
        }
    }
    return j;
}

// ─────────────────────────────────────────────────────────────────────────────
// Consulta recursiva
// ─────────────────────────────────────────────────────────────────────────────

template <bool t_min>
size_t WmRmq<t_min>::query_rec(uint32_t k,
                                size_t cur_lo, size_t cur_hi,
                                size_t sym_lo, size_t sym_hi,
                                size_t y_lo,   size_t y_hi,
                                const sdsl::int_vector<>& text_pos,
                                size_t& count_acc) const {
    if (cur_lo > cur_hi)                    return SIZE_MAX;
    if (sym_hi < y_lo || sym_lo > y_hi)    return SIZE_MAX;

    if (y_lo <= sym_lo && sym_hi <= y_hi) {
        const size_t j_local  = rmqs_[k](cur_lo, cur_hi);
        const size_t j_global = unwind(j_local, k);
        count_acc += cur_hi - cur_lo + 1;
        return j_global;
    }

    assert(k < wm_.max_level && "hoja sin cobertura total — rango y inválido");

    const size_t kn     = k * n_;
    const size_t rl_lo  = wm_.tree_rank1(kn + cur_lo)       - wm_.rank_level_at(k);
    const size_t rl_hi1 = wm_.tree_rank1(kn + cur_hi + 1)   - wm_.rank_level_at(k);
    const size_t r0_lo  = cur_lo  - rl_lo;
    const size_t r0_hi1 = (cur_hi + 1) - rl_hi1;
    const size_t zc     = wm_.zero_cnt_at(k);

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
    if constexpr (t_min)
        return (text_pos[cand_L] <= text_pos[cand_R]) ? cand_L : cand_R;
    else
        return (text_pos[cand_L] >= text_pos[cand_R]) ? cand_L : cand_R;
}

template <bool t_min>
typename WmRmq<t_min>::ArgminResult
WmRmq<t_min>::range_argmin_2d(size_t x_lo, size_t x_hi,
                               size_t y_lo, size_t y_hi,
                               const sdsl::int_vector<>& text_pos) const {
    if (n_ == 0 || wm_.empty() || x_lo > x_hi) return {0, SIZE_MAX};
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

template <bool t_min>
typename WmRmq<t_min>::SizeBreakdown WmRmq<t_min>::size_breakdown() const {
    const size_t bv_bytes  = sdsl::size_in_bytes(wm_.tree);
    const size_t wm_total  = sdsl::size_in_bytes(wm_);
    const size_t rs_bytes  = wm_total > bv_bytes ? wm_total - bv_bytes : 0;
    size_t rmq_bytes = 0;
    for (const auto& r : rmqs_) rmq_bytes += sdsl::size_in_bytes(r);
    return {bv_bytes, rs_bytes, rmq_bytes};
}

template <bool t_min>
size_t WmRmq<t_min>::size_in_bytes() const {
    const auto bd = size_breakdown();
    return bd.bv + bd.rank_sel + bd.rmq;
}

// ─────────────────────────────────────────────────────────────────────────────
// Instanciaciones explícitas
// ─────────────────────────────────────────────────────────────────────────────

template class WmRmq<true>;
template class WmRmq<false>;

}  // namespace lz77tax
