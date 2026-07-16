#include "grid.hpp"

#include <algorithm>
#include <cassert>

#include <sdsl/construct.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/bits.hpp>
#include <sdsl/io.hpp>

namespace lz77tax {

// -----------------------------------------------------------------------------
// build -- R y arrays ya en orden X (0-indexed)
// -----------------------------------------------------------------------------

void Grid2D::build(const std::vector<size_t>& R,
                   const std::vector<size_t>& boundaries,
                   const std::vector<size_t>& phrase_lens) {
    const size_t z1 = R.size();
    assert(boundaries.size() == z1 && phrase_lens.size() == z1);
    if (z1 == 0) return;

    const uint8_t width = sdsl::bits::hi(z1) + 1;
    sdsl::int_vector<> Rv(z1, 0, width);
    for (size_t j = 0; j < z1; ++j) Rv[j] = R[j];

    const size_t n_val = *std::max_element(boundaries.begin(), boundaries.end());
    const uint8_t tp_w = n_val > 0 ? static_cast<uint8_t>(sdsl::bits::hi(n_val) + 1) : 1;
    text_pos_ = sdsl::int_vector<>(z1, 0, tp_w);
    for (size_t j = 0; j < z1; ++j) text_pos_[j] = boundaries[j];

    const size_t max_pl = *std::max_element(phrase_lens.begin(), phrase_lens.end());
    const uint8_t pl_w = max_pl > 0 ? static_cast<uint8_t>(sdsl::bits::hi(max_pl) + 1) : 1;
    phrase_total_len_ = sdsl::int_vector<>(z1, 0, pl_w);
    for (size_t j = 0; j < z1; ++j) phrase_total_len_[j] = phrase_lens[j];

    sdsl::construct_im(static_cast<sdsl::wm_int<>&>(wm_), Rv);

    // Las coordenadas Y de la grilla son una permutacion (una por frase), de modo
    // que el nivel hoja de la wavelet matrix no necesita RMQ.
    wm_min_rmq_.build(wm_, text_pos_, /*distinct_symbols=*/true);
}

// -----------------------------------------------------------------------------
// Consultas -- rangos 0-indexed ya relativos
// -----------------------------------------------------------------------------

Grid2D::RangeResult Grid2D::query(size_t x_lo, size_t x_hi,
                                  size_t y_lo, size_t y_hi) const {
    if (wm_.size() == 0 || x_lo > x_hi || y_lo > y_hi) return {0, {}};
    if (x_hi >= wm_.size() || y_hi >= wm_.size()) return {0, {}};
    return wm_.range_search_2d(x_lo, x_hi, y_lo, y_hi);
}

Grid2D::MinResult Grid2D::query_min(size_t x_lo, size_t x_hi,
                                    size_t y_lo, size_t y_hi) const {
    if (wm_min_rmq_.size() == 0 || x_lo > x_hi || y_lo > y_hi) return {0, SIZE_MAX, SIZE_MAX};
    if (x_hi >= wm_.size() || y_hi >= wm_.size()) return {0, SIZE_MAX, SIZE_MAX};
    const auto r = wm_min_rmq_.range_argmin_2d(wm_, x_lo, x_hi, y_lo, y_hi, text_pos_);
    if (r.count == 0) return {0, SIZE_MAX, SIZE_MAX};
    return {r.count, static_cast<size_t>(text_pos_[r.argmin_global]), r.argmin_global};
}

Grid2D::MinResult Grid2D::query_min_filtered(size_t x_lo, size_t x_hi,
                                             size_t y_lo, size_t y_hi,
                                             size_t min_phrase_len) const {
    const auto res = query(x_lo, x_hi, y_lo, y_hi);
    if (res.first == 0) return {0, SIZE_MAX, SIZE_MAX};

    size_t count = 0;
    size_t best_boundary = SIZE_MAX;
    size_t best_idx = SIZE_MAX;
    for (const auto& [wt_idx, y_rel] : res.second) {
        (void)y_rel;
        if (phrase_total_len_[wt_idx] < min_phrase_len) continue;
        ++count;
        const size_t boundary = text_pos_[wt_idx];
        if (boundary < best_boundary) {
            best_boundary = boundary;
            best_idx = wt_idx;
        }
    }
    if (count == 0) return {0, SIZE_MAX, SIZE_MAX};
    return {count, best_boundary, best_idx};
}

Grid2D::SpecialResult Grid2D::query_special(size_t y_lo, size_t y_hi,
                                            size_t plen) const {
    if (wm_.size() == 0 || y_lo > y_hi) return {0, SIZE_MAX, 0};
    if (y_hi >= wm_.size()) return {0, SIZE_MAX, 0};

    const auto pts = wm_.range_search_2d(0, wm_.size() - 1, y_lo, y_hi).second;

    size_t count = 0;
    size_t occ_min = SIZE_MAX;
    size_t occ_max = 0;

    for (const auto& [wt_idx, y_rel] : pts) {
        (void)y_rel;
        if (phrase_total_len_[wt_idx] >= plen) {
            ++count;
            const size_t occ_pos = text_pos_[wt_idx] - plen;
            if (occ_pos < occ_min) occ_min = occ_pos;
            if (occ_pos > occ_max) occ_max = occ_pos;
        }
    }
    return {count, occ_min, occ_max};
}

// -----------------------------------------------------------------------------
// Serializacion
// -----------------------------------------------------------------------------

size_t Grid2D::serialize(std::ostream& out) const {
    // phrase_total_len_ NO se serializa: es derivable de text_pos_ (largo de la
    // frase k = inicio_{k+1} - inicio_k = diferencia entre inicios de frase
    // consecutivos, y los inicios son {0} U text_pos). Se reconstruye en load().
    size_t written = 0;
    written += sdsl::serialize(static_cast<const sdsl::wm_int<>&>(wm_), out);
    written += sdsl::serialize(text_pos_,        out);
    written += wm_min_rmq_.serialize(out);
    return written;
}

void Grid2D::load(std::istream& in) {
    sdsl::load(static_cast<sdsl::wm_int<>&>(wm_), in);
    sdsl::load(text_pos_,         in);
    reconstruct_phrase_total_len();
    wm_min_rmq_.load(in);
}

// phrase_total_len_[i] = text_pos_[i] - (mayor inicio de frase < text_pos_[i]).
// Los inicios de frase son {0} U {text_pos_}; en orden creciente, el predecesor
// de text_pos_[i] es el valor inmediatamente menor (o 0 para el mínimo).
void Grid2D::reconstruct_phrase_total_len() {
    const size_t z1 = text_pos_.size();
    if (z1 == 0) { phrase_total_len_ = sdsl::int_vector<>(); return; }

    // Orden por valor de text_pos, guardando el índice original.
    std::vector<std::pair<size_t, size_t>> sv(z1);   // (valor, índice)
    for (size_t i = 0; i < z1; ++i) sv[i] = {static_cast<size_t>(text_pos_[i]), i};
    std::sort(sv.begin(), sv.end());

    std::vector<size_t> plen(z1);
    size_t prev = 0, max_pl = 0;                     // inicio de frase 0 = pos 0
    for (const auto& [val, idx] : sv) {
        const size_t len = val - prev;               // largo total de la frase
        plen[idx] = len;
        max_pl = std::max(max_pl, len);
        prev = val;
    }
    const uint8_t w = max_pl > 0 ? static_cast<uint8_t>(sdsl::bits::hi(max_pl) + 1) : 1;
    phrase_total_len_ = sdsl::int_vector<>(z1, 0, w);
    for (size_t i = 0; i < z1; ++i) phrase_total_len_[i] = plen[i];
}

Grid2D::SizeBreakdown Grid2D::size_breakdown() const {
    SizeBreakdown bd{};
    bd.wm              = sdsl::size_in_bytes(wm_);
    bd.text_pos        = sdsl::size_in_bytes(text_pos_);
    bd.phrase_total_len= sdsl::size_in_bytes(phrase_total_len_);
    bd.wm_min_rmq      = wm_min_rmq_.size_in_bytes();
    return bd;
}

}  // namespace lz77tax
