#include "wt_rmq_min.hpp"

#include <algorithm>
#include <cassert>
#include <climits>

namespace lz77tax {

// ─────────────────────────────────────────────────────────────────────────────
// Construcción
// ─────────────────────────────────────────────────────────────────────────────

void WtMinRmq::build(const std::vector<size_t>& y_values,
                     const sdsl::int_vector<>& text_pos,
                     size_t sigma) {
    assert(y_values.size() == text_pos.size());
    n_     = y_values.size();
    sigma_ = sigma;
    root_.reset();
    if (n_ == 0 || sigma_ == 0) return;
    // Convertir a vector<size_t> para build_rec (temporal, liberado al salir)
    std::vector<size_t> tps_init(text_pos.begin(), text_pos.end());
    root_ = build_rec(y_values, std::move(tps_init), 0, sigma_ - 1,
                      /*parent=*/nullptr, /*is_right=*/false);
}

std::unique_ptr<WtMinRmq::Node>
WtMinRmq::build_rec(std::vector<size_t> ys,
                    std::vector<size_t> tps,
                    size_t lo, size_t hi,
                    Node* parent, bool is_right) {
    auto node = std::make_unique<Node>();
    node->lo             = lo;
    node->hi             = hi;
    node->parent         = parent;
    node->is_right_child = is_right;

    // RMQ topológico: construido desde int_vector temporal que se libera al
    // salir del bloque. rmq_succinct_sct solo retiene el Cartesian tree en BP
    // (2n+o(n) bits), sin los valores originales.
    if (!tps.empty()) {
        size_t max_v = *std::max_element(tps.begin(), tps.end());
        const uint8_t w = max_v ? static_cast<uint8_t>(sdsl::bits::hi(max_v) + 1) : 1;
        sdsl::int_vector<> iv(tps.size(), 0, w);
        for (size_t i = 0; i < tps.size(); ++i) iv[i] = tps[i];
        node->rmq = sdsl::rmq_succinct_sct<true>(&iv);
        // iv se destruye aquí — rmq no lo retiene
    }

    if (lo == hi || ys.size() <= 1) return node;  // hoja

    const size_t mid = lo + (hi - lo) / 2;
    node->bv = sdsl::bit_vector(ys.size(), 0);

    std::vector<size_t> ys_l, tps_l, ys_r, tps_r;
    ys_l.reserve(ys.size());  tps_l.reserve(ys.size());
    ys_r.reserve(ys.size());  tps_r.reserve(ys.size());

    for (size_t i = 0; i < ys.size(); ++i) {
        if (ys[i] <= mid) {
            node->bv[i] = 0;
            ys_l.push_back(ys[i]);  tps_l.push_back(tps[i]);
        } else {
            node->bv[i] = 1;
            ys_r.push_back(ys[i]);  tps_r.push_back(tps[i]);
        }
    }

    sdsl::util::init_support(node->rank1, &node->bv);
    sdsl::util::init_support(node->sel0,  &node->bv);
    sdsl::util::init_support(node->sel1,  &node->bv);

    if (!ys_l.empty())
        node->left  = build_rec(std::move(ys_l), std::move(tps_l),
                                lo, mid, node.get(), false);
    if (!ys_r.empty())
        node->right = build_rec(std::move(ys_r), std::move(tps_r),
                                mid + 1, hi, node.get(), true);
    return node;
}

// ─────────────────────────────────────────────────────────────────────────────
// Unwind (nodo, índice local) → índice global en raíz. O(log σ).
// ─────────────────────────────────────────────────────────────────────────────

size_t WtMinRmq::unwind_to_root(const Node* node, size_t local_j) {
    size_t j = local_j;
    while (node->parent) {
        const Node* p = node->parent;
        // select es 1-indexed en sdsl; retorna posición 0-indexed en bv del padre.
        j = node->is_right_child ? p->sel1(j + 1) : p->sel0(j + 1);
        node = p;
    }
    return j;
}

// ─────────────────────────────────────────────────────────────────────────────
// Consulta recursiva
// ─────────────────────────────────────────────────────────────────────────────

void WtMinRmq::query_rec(const Node* node,
                         size_t x_lo, size_t x_hi,
                         size_t y_lo, size_t y_hi,
                         const sdsl::int_vector<>& tp_global,
                         size_t& count_acc,
                         size_t& best_value,
                         size_t& best_global) const {
    if (!node || x_lo > x_hi) return;
    // Nodo completamente fuera del rango de y
    if (node->hi < y_lo || node->lo > y_hi) return;

    if (y_lo <= node->lo && node->hi <= y_hi) {
        // Nodo fully-covered: RMQ topológico en O(1), luego unwind en O(log σ)
        const size_t j_local  = node->rmq(x_lo, x_hi);
        const size_t j_global = unwind_to_root(node, j_local);
        const size_t v        = tp_global[j_global];
        count_acc += x_hi - x_lo + 1;
        if (v < best_value) {
            best_value  = v;
            best_global = j_global;
        }
        return;
    }

    if (!node->left && !node->right) return;  // hoja fuera de rango

    // Proyectar rango X a hijos via rank
    const size_t r1_lo = node->rank1(x_lo);
    const size_t r1_hi = node->rank1(x_hi + 1);
    const size_t r0_lo = x_lo  - r1_lo;
    const size_t r0_hi = (x_hi + 1) - r1_hi;

    if (node->left && r0_lo < r0_hi)
        query_rec(node->left.get(),  r0_lo, r0_hi - 1,
                  y_lo, y_hi, tp_global, count_acc, best_value, best_global);
    if (node->right && r1_lo < r1_hi)
        query_rec(node->right.get(), r1_lo, r1_hi - 1,
                  y_lo, y_hi, tp_global, count_acc, best_value, best_global);
}

WtMinRmq::ArgminResult
WtMinRmq::range_argmin_2d(size_t x_lo, size_t x_hi,
                          size_t y_lo, size_t y_hi,
                          const sdsl::int_vector<>& text_pos) const {
    if (n_ == 0 || !root_ || x_lo > x_hi) return {0, SIZE_MAX};
    size_t count = 0, best_v = SIZE_MAX, best_i = SIZE_MAX;
    query_rec(root_.get(), x_lo, x_hi, y_lo, y_hi,
              text_pos, count, best_v, best_i);
    return {count, best_i};
}

// ─────────────────────────────────────────────────────────────────────────────
// Métricas de tamaño
// ─────────────────────────────────────────────────────────────────────────────

void WtMinRmq::accumulate_size(const Node* node,
                                size_t& bv_bytes,
                                size_t& rank_sel_bytes,
                                size_t& rmq_bytes) {
    if (!node) return;
    bv_bytes       += sdsl::size_in_bytes(node->bv);
    rank_sel_bytes += sdsl::size_in_bytes(node->rank1)
                   +  sdsl::size_in_bytes(node->sel0)
                   +  sdsl::size_in_bytes(node->sel1);
    rmq_bytes      += sdsl::size_in_bytes(node->rmq);
    accumulate_size(node->left.get(),  bv_bytes, rank_sel_bytes, rmq_bytes);
    accumulate_size(node->right.get(), bv_bytes, rank_sel_bytes, rmq_bytes);
}

size_t WtMinRmq::size_in_bytes() const {
    size_t bv = 0, rs = 0, rmq = 0;
    accumulate_size(root_.get(), bv, rs, rmq);
    return bv + rs + rmq;
}

WtMinRmq::SizeBreakdown WtMinRmq::size_breakdown() const {
    SizeBreakdown bd{};
    accumulate_size(root_.get(), bd.bv, bd.rank_sel, bd.rmq);
    return bd;
}

}  // namespace lz77tax
