#include "wm_rmq_min.hpp"

#include <algorithm>
#include <cassert>
#include <climits>
#include <numeric>

#include <sdsl/bits.hpp>
#include <sdsl/util.hpp>

namespace lz77tax {

// ─────────────────────────────────────────────────────────────────────────────
// Construcción
// ─────────────────────────────────────────────────────────────────────────────

template <bool t_min>
void WmRmq<t_min>::build_level_rmqs(const SharedWm& wm,
                                    const sdsl::int_vector<>& text_pos,
                                    bool distinct_symbols) {
    const uint32_t L = wm.max_level;
    // Con simbolos distintos, el nivel L cubre un solo punto por nodo virtual y
    // su RMQ es innecesario (el argmin es el propio elemento): se omite. Si hay
    // simbolos repetidos el nodo hoja puede tener varios puntos y hace falta.
    rmqs_.resize(distinct_symbols ? L : L + 1);

    // Trazar la permutación de índices que induce la wavelet matrix nivel a nivel.
    std::vector<size_t> perm(n_);
    std::iota(perm.begin(), perm.end(), 0);

    for (uint32_t k = 0; k <= L; ++k) {
        // Construir RMQ del nivel k sobre los pesos reordenados por perm.
        size_t max_w = 0;
        for (size_t i = 0; i < n_; ++i)
            max_w = std::max(max_w, static_cast<size_t>(text_pos[perm[i]]));

        const uint8_t w = max_w ? static_cast<uint8_t>(sdsl::bits::hi(max_w) + 1) : 1;
        sdsl::int_vector<> weights(n_, 0, w);
        for (size_t i = 0; i < n_; ++i)
            weights[i] = text_pos[perm[i]];

        if (k < rmqs_.size()) rmqs_[k] = sdsl::rmq_succinct_sct<t_min>(&weights);

        if (k == L) break;

        // Avanzar perm al siguiente nivel: primero los 0-bits, luego los 1-bits.
        std::vector<size_t> left_perm, right_perm;
        left_perm.reserve(n_);
        right_perm.reserve(n_);
        const size_t kn = k * n_;
        for (size_t i = 0; i < n_; ++i) {
            if (!wm.tree[kn + i])
                left_perm.push_back(perm[i]);
            else
                right_perm.push_back(perm[i]);
        }
        perm.clear();
        perm.insert(perm.end(), left_perm.begin(),  left_perm.end());
        perm.insert(perm.end(), right_perm.begin(), right_perm.end());
    }
}

template <bool t_min>
void WmRmq<t_min>::build(const SharedWm& wm, const sdsl::int_vector<>& text_pos,
                         bool distinct_symbols) {
    assert(wm.size() == text_pos.size());
    n_     = wm.size();
    sigma_ = n_;
    rmqs_.clear();
    if (n_ == 0 || wm.empty()) return;

    build_level_rmqs(wm, text_pos, distinct_symbols);
}

// ─────────────────────────────────────────────────────────────────────────────
// Unwind: posición local j al nivel from_depth → índice original (nivel 0).
// ─────────────────────────────────────────────────────────────────────────────

template <bool t_min>
size_t WmRmq<t_min>::unwind(const SharedWm& wm, size_t j, uint32_t from_depth) const {
    for (uint32_t d = from_depth; d >= 1; --d) {
        const size_t parent = d - 1;
        const size_t pn     = parent * n_;
        if (j < wm.zero_cnt_at(parent)) {
            const size_t zeros_before = pn - wm.rank_level_at(parent);
            const size_t abs = wm.tree_sel0(zeros_before + j + 1);
            j = abs - pn;
        } else {
            const size_t rank_ones = j - wm.zero_cnt_at(parent);
            const size_t abs = wm.tree_sel1(wm.rank_level_at(parent) + rank_ones + 1);
            j = abs - pn;
        }
    }
    return j;
}

// ─────────────────────────────────────────────────────────────────────────────
// Consulta recursiva
// ─────────────────────────────────────────────────────────────────────────────

template <bool t_min>
size_t WmRmq<t_min>::query_rec(const SharedWm& wm,
                                uint32_t k,
                                size_t cur_lo, size_t cur_hi,
                                size_t sym_lo, size_t sym_hi,
                                size_t y_lo,   size_t y_hi,
                                const sdsl::int_vector<>& text_pos,
                                size_t& count_acc) const {
    if (cur_lo > cur_hi)                    return SIZE_MAX;
    if (sym_hi < y_lo || sym_lo > y_hi)    return SIZE_MAX;

    if (y_lo <= sym_lo && sym_hi <= y_hi) {
        // En el ultimo nivel el nodo cubre 1 simbolo = 1 punto (Y es permutacion):
        // el argmin es el propio elemento, sin consultar estructura alguna.
        // Sin RMQ para el nivel hoja (simbolos distintos): el nodo cubre un
        // unico punto y el argmin es el propio elemento.
        const size_t j_local = (k < rmqs_.size()) ? rmqs_[k](cur_lo, cur_hi)
                                                  : cur_lo;
        assert((k < rmqs_.size() || cur_lo == cur_hi) &&
               "nivel hoja sin RMQ con mas de un punto");
        const size_t j_global = unwind(wm, j_local, k);
        count_acc += cur_hi - cur_lo + 1;
        return j_global;
    }

    assert(k < wm.max_level && "hoja sin cobertura total — rango y inválido");

    const size_t kn     = k * n_;
    const size_t rl_lo  = wm.tree_rank1(kn + cur_lo)       - wm.rank_level_at(k);
    const size_t rl_hi1 = wm.tree_rank1(kn + cur_hi + 1)   - wm.rank_level_at(k);
    const size_t r0_lo  = cur_lo  - rl_lo;
    const size_t r0_hi1 = (cur_hi + 1) - rl_hi1;
    const size_t zc     = wm.zero_cnt_at(k);

    const size_t bit_w = size_t(1) << (wm.max_level - 1 - k);
    const size_t mid   = sym_lo + bit_w - 1;

    size_t cand_L = SIZE_MAX, cand_R = SIZE_MAX;

    if (y_lo <= mid && r0_lo < r0_hi1)
        cand_L = query_rec(wm, k + 1, r0_lo,      r0_hi1 - 1,
                           sym_lo, mid,        y_lo, y_hi,
                           text_pos, count_acc);

    if (y_hi > mid && rl_lo < rl_hi1)
        cand_R = query_rec(wm, k + 1, zc + rl_lo, zc + rl_hi1 - 1,
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
WmRmq<t_min>::range_argmin_2d(const SharedWm& wm,
                               size_t x_lo, size_t x_hi,
                               size_t y_lo, size_t y_hi,
                               const sdsl::int_vector<>& text_pos) const {
    if (n_ == 0 || wm.empty() || x_lo > x_hi) return {0, SIZE_MAX};
    const size_t sym_hi_root = (size_t(1) << wm.max_level) - 1;
    size_t count = 0;
    const size_t best_global = query_rec(wm, 0, x_lo, x_hi,
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
    size_t rmq_bytes = 0;
    for (const auto& r : rmqs_) rmq_bytes += sdsl::size_in_bytes(r);
    return {0, 0, rmq_bytes};
}

template <bool t_min>
size_t WmRmq<t_min>::size_in_bytes() const {
    const auto bd = size_breakdown();
    return bd.bv + bd.rank_sel + bd.rmq;
}

// ─────────────────────────────────────────────────────────────────────────────
// Serialización
// ─────────────────────────────────────────────────────────────────────────────

template <bool t_min>
size_t WmRmq<t_min>::serialize(std::ostream& out) const {
    size_t written = 0;
    written += sdsl::write_member(n_,     out);
    written += sdsl::write_member(sigma_, out);
    const uint64_t nrmq = rmqs_.size();
    written += sdsl::write_member(nrmq, out);
    for (const auto& r : rmqs_)
        written += sdsl::serialize(r, out);
    return written;
}

template <bool t_min>
void WmRmq<t_min>::load(std::istream& in) {
    sdsl::read_member(n_,     in);
    sdsl::read_member(sigma_, in);
    uint64_t nrmq = 0;
    sdsl::read_member(nrmq, in);
    rmqs_.resize(nrmq);
    for (auto& r : rmqs_)
        sdsl::load(r, in);
}

// ─────────────────────────────────────────────────────────────────────────────
// Instanciaciones explícitas
// ─────────────────────────────────────────────────────────────────────────────

// Solo mínimo: con el sub-índice sobre T^R el máximo ES un mínimo.
template class WmRmq<true>;

}  // namespace lz77tax
