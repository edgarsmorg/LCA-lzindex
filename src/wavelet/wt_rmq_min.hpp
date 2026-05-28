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
 * Resuelve argmin (t_min=true) o argmax (t_min=false) 2D en O(log² σ).
 * Espacio O(z log σ) bits — un nodo sdsl por elemento.
 * Variante legacy: solo para benchmarks comparativos.
 */
template <bool t_min = true>
class WtRmq {
public:
    struct ArgminResult {
        size_t count;
        size_t argmin_global; ///< SIZE_MAX si count==0
    };

    struct SizeBreakdown { size_t bv; size_t rank_sel; size_t rmq; };

    void build(const std::vector<size_t>& y_values,
               const sdsl::int_vector<>& text_pos,
               size_t sigma);

    ArgminResult range_argmin_2d(size_t x_lo, size_t x_hi,
                                 size_t y_lo, size_t y_hi,
                                 const sdsl::int_vector<>& text_pos) const;

    size_t size()          const { return n_; }
    size_t alphabet_size() const { return sigma_; }
    size_t size_in_bytes() const;
    SizeBreakdown size_breakdown() const;

    WtRmq() = default;
    WtRmq(const WtRmq&)            = delete;
    WtRmq& operator=(const WtRmq&) = delete;
    WtRmq(WtRmq&&)                 = default;
    WtRmq& operator=(WtRmq&&)      = default;

private:
    struct Node {
        sdsl::bit_vector                  bv;
        sdsl::rank_support_v<1>           rank1;
        sdsl::select_support_mcl<0>       sel0;
        sdsl::select_support_mcl<1>       sel1;
        sdsl::rmq_succinct_sct<t_min>     rmq;
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
                   const sdsl::int_vector<>& tp_global,
                   size_t& count_acc,
                   size_t& best_value,
                   size_t& best_global) const;

    static void accumulate_size(const Node* node,
                                size_t& bv_bytes,
                                size_t& rank_sel_bytes,
                                size_t& rmq_bytes);
};

using WtMinRmq = WtRmq<true>;
using WtMaxRmq = WtRmq<false>;

}  // namespace lz77tax
