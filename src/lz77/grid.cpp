#include "grid.hpp"

#include <algorithm>
#include <cassert>
#include <stdexcept>

#include <sdsl/construct.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/bits.hpp>
#include <sdsl/io.hpp>
#include <sdsl/util.hpp>

namespace lz77tax {

// -----------------------------------------------------------------------------
// build_core -- nucleo sin bv_fwd_/bv_rev_ (R y arrays ya en orden WT)
// -----------------------------------------------------------------------------

void Grid2D::build_core(const std::vector<size_t>& R,
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

    wm_min_rmq_.build(wm_, text_pos_);
}

// -----------------------------------------------------------------------------
// build_from_trie_coords -- produccion: coordenadas ya relativas (trie-based)
// -----------------------------------------------------------------------------

void Grid2D::build_from_trie_coords(const std::vector<size_t>& R,
                                     const std::vector<size_t>& boundaries,
                                     const std::vector<size_t>& phrase_lens) {
    // bv_fwd_/bv_rev_ se dejan vacias: las queries usan rangos ya relativos.
    build_core(R, boundaries, phrase_lens);
}

// -----------------------------------------------------------------------------
// query_direct / query_min_direct / query_max_direct / query_special_direct
// Toman rangos 0-indexed ya relativos (sin clamp_ranges).
// -----------------------------------------------------------------------------

Grid2D::RangeResult Grid2D::query_direct(size_t x_lo, size_t x_hi,
                                          size_t y_lo, size_t y_hi) const {
    if (wm_.size() == 0 || x_lo > x_hi || y_lo > y_hi) return {0, {}};
    if (x_hi >= wm_.size() || y_hi >= wm_.size()) return {0, {}};
    return wm_.range_search_2d(x_lo, x_hi, y_lo, y_hi);
}

Grid2D::MinResult Grid2D::query_min_direct(size_t x_lo, size_t x_hi,
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
    const auto res = query_direct(x_lo, x_hi, y_lo, y_hi);
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

Grid2D::SpecialResult Grid2D::query_special_direct(size_t y_lo, size_t y_hi,
                                                    size_t plen) const {
    if (wm_.size() == 0 || y_lo > y_hi) return {0, SIZE_MAX, 0};
    if (y_hi >= wm_.size()) return {0, SIZE_MAX, 0};

    const size_t lb = 0;
    const size_t rb = wm_.size() - 1;

    const auto pts = wm_.range_search_2d(lb, rb, y_lo, y_hi).second;

    size_t count = 0;
    size_t occ_min = SIZE_MAX;
    size_t occ_max = 0;

    for (const auto& [wt_idx, y_rel] : pts) {
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
// Nucleo compartido (ISA-based, preservado para tests)
// -----------------------------------------------------------------------------

void Grid2D::build_from_coords(
        const std::vector<std::pair<size_t, size_t>>& coords,
        const std::vector<size_t>& boundaries,
        const std::vector<size_t>& phrase_lens,
        size_t n) {

    const size_t z1 = coords.size();  // z-1 puntos
    assert(boundaries.size() == z1 && "coords y boundaries deben tener el mismo tamanio");
    assert(phrase_lens.size() == z1 && "coords y phrase_lens deben tener el mismo tamanio");

    if (z1 == 0) {
        // Texto de 1 frase: grilla vacia, estructuras en estado default.
        return;
    }

    // 1. Bitvector forward: un 1 en cada X_k
    {
        sdsl::bit_vector bfwd(n, 0);
        for (const auto& [x, y] : coords) {
            assert(x < n && "X_k fuera del rango del SA");
            bfwd[x] = 1;
        }
        bv_fwd_ = sdsl::sd_vector<>(bfwd);
    }
    rank_fwd_ = sdsl::sd_vector<>::rank_1_type(&bv_fwd_);

    // 2. Bitvector reverse: un 1 en cada Y_k
    {
        sdsl::bit_vector brev(n, 0);
        for (const auto& [x, y] : coords) {
            assert(y < n && "Y_k fuera del rango del rev-SA");
            brev[y] = 1;
        }
        bv_rev_ = sdsl::sd_vector<>(brev);
    }
    rank_rev_ = sdsl::sd_vector<>::rank_1_type(&bv_rev_);

    // 3. Secuencia R[0..z-2]: R[rank_fwd(X_k)] = rank_rev(Y_k)
    // boundaries_wt y phrase_lens_wt reordenados por rank_fwd (orden WT).
    std::vector<size_t> R_vec(z1), boundaries_wt(z1), phrase_lens_wt(z1);
    for (size_t k = 0; k < z1; ++k) {
        const size_t x = coords[k].first;
        const size_t y = coords[k].second;
        const size_t j     = rank_fwd_(x);
        const size_t y_rel = rank_rev_(y);
        assert(j < z1 && "indice wt_ fuera de rango");
        assert(y_rel < z1 && "y_rel fuera de rango");
        R_vec[j]          = y_rel;
        boundaries_wt[j]  = boundaries[k];
        phrase_lens_wt[j] = phrase_lens[k];
    }

    // 4. Construir wm_int, text_pos_, phrase_total_len_, WmMinRmq
    build_core(R_vec, boundaries_wt, phrase_lens_wt);
}

// -----------------------------------------------------------------------------
// build() desde FM-index (produccion, GB+)
// -----------------------------------------------------------------------------

void Grid2D::build(const LZ77Parsing& phrases,
                   const sdsl::int_vector<>& isa_fwd,
                   const sdsl::int_vector<>& isa_rev,
                   size_t n) {
    const size_t z = phrases.size();

    if (z <= 1) return;

    // Para k = 0..z-2:
    //   end_k    = phrases[k].start_pos + phrases[k].length
    //   start_k1 = phrases[k+1].start_pos  (== end_k + 1)
    //   X_k = ISA_fwd[ start_k1 ]
    //   Y_k = ISA_rev[ n-2 - end_k ]
    //
    // La formula Y_k = n-2-end_k viene de que text_rev[j] = text[n-2-j],
    // por lo que para buscar text[end_k] en el SA reverso usamos j = n-2-end_k.
    std::vector<std::pair<size_t, size_t>> coords(z - 1);
    std::vector<size_t> boundaries(z - 1);
    std::vector<size_t> phrase_lens(z - 1);
    for (size_t k = 0; k < z - 1; ++k) {
        const size_t end_k    = phrases[k].start_pos + phrases[k].length;
        const size_t start_k1 = phrases[k + 1].start_pos;

        assert(start_k1 == end_k + 1 &&
               "gap entre frases consecutivas: start_k1 != end_k + 1");
        assert(start_k1 < n    && "start_k1 fuera del texto");
        assert(end_k < n - 1   && "end_k apunta al centinela");
        assert(n >= 2 + end_k  && "n-2-end_k subflow");

        coords[k]      = { isa_fwd[start_k1], isa_rev[n - 2 - end_k] };
        boundaries[k]  = start_k1;
        phrase_lens[k] = phrases[k].length + 1;
    }

    build_from_coords(coords, boundaries, phrase_lens, n);
}

// -----------------------------------------------------------------------------
// build() desde texto crudo (small-scale / tests)
// -----------------------------------------------------------------------------

void Grid2D::build(const LZ77Parsing& phrases, const std::string& text) {
    const size_t n = text.size();
    const size_t z = phrases.size();

    if (z <= 1) return;

    // SA y ISA del texto
    sdsl::int_vector<> sa(n);
    sdsl::algorithm::calculate_sa(
            reinterpret_cast<const unsigned char*>(text.data()), n, sa);

    sdsl::int_vector<> isa(n);
    for (size_t i = 0; i < n; ++i) isa[sa[i]] = i;

    // SA y ISA del texto reverso.
    // text termina con el centinela '\0' en posicion n-1. Saltamos ese centinela
    // y revertimos solo text[0..n-2], luego aniadimos '\0' al final.
    // Resultado: text_rev_s[j] = text[n-2-j] para j=0..n-2, text_rev_s[n-1]='\0'.
    const std::string text_rev_sa(text.rbegin() + 1, text.rend());
    const std::string text_rev_sa_s = text_rev_sa + '\0';

    sdsl::int_vector<> sa_rev(n);
    sdsl::algorithm::calculate_sa(
            reinterpret_cast<const unsigned char*>(text_rev_sa_s.data()), n, sa_rev);

    sdsl::int_vector<> isa_rev(n);
    for (size_t i = 0; i < n; ++i) isa_rev[sa_rev[i]] = i;

    // Delegar la construccion de coordenadas al overload de produccion.
    build(phrases, isa, isa_rev, n);
}

// -----------------------------------------------------------------------------
// clamp_ranges() -- helper privado compartido por todos los metodos de consulta
// -----------------------------------------------------------------------------

// Convierte rangos BWT globales [sp_r,ep_r] x [sp_l,ep_l] a rangos relativos
// en wm_. rank_1_type(i) = #{j < i : bv[j]=1} (exclusivo en i, 0-indexed).
// Para inputs validos ep_* < n = bv_*.size(); las guardas de desbordamiento son
// defensivas ante llamadas con rangos invalidos.
// Devuelve false si el rectangulo resultante esta vacio.
bool Grid2D::clamp_ranges(size_t sp_r, size_t ep_r,
                           size_t sp_l, size_t ep_l,
                           size_t& lb, size_t& rb,
                           size_t& vlb, size_t& vrb) const {
    lb  = rank_fwd_(sp_r);
    if (ep_r + 1 > bv_fwd_.size()) return false;
    rb  = rank_fwd_(ep_r + 1);
    if (rb == 0 || lb >= rb) return false;

    vlb = rank_rev_(sp_l);
    if (ep_l + 1 > bv_rev_.size()) return false;
    vrb = rank_rev_(ep_l + 1);
    if (vrb == 0 || vlb >= vrb) return false;

    return true;
}

// -----------------------------------------------------------------------------
// query()
// -----------------------------------------------------------------------------

Grid2D::RangeResult Grid2D::query(size_t sp_right, size_t ep_right,
                                   size_t sp_left,  size_t ep_left) const {
    if (wm_.size() == 0) return {0, {}};

    size_t lb, rb, vlb, vrb;
    if (!clamp_ranges(sp_right, ep_right, sp_left, ep_left, lb, rb, vlb, vrb))
        return {0, {}};

    return wm_.range_search_2d(lb, rb - 1, vlb, vrb - 1);
}

// -----------------------------------------------------------------------------
// query_min_2d()
// -----------------------------------------------------------------------------

Grid2D::MinResult Grid2D::query_min_2d(size_t sp_right, size_t ep_right,
                                        size_t sp_left,  size_t ep_left) const {
    size_t lb, rb, vlb, vrb;
    if (!clamp_ranges(sp_right, ep_right, sp_left, ep_left, lb, rb, vlb, vrb))
        return {0, SIZE_MAX, SIZE_MAX};

    if (wm_min_rmq_.size() == 0) return {0, SIZE_MAX, SIZE_MAX};
    const auto r = wm_min_rmq_.range_argmin_2d(wm_, lb, rb - 1, vlb, vrb - 1, text_pos_);
    if (r.count == 0) return {0, SIZE_MAX, SIZE_MAX};
    return {r.count, static_cast<size_t>(text_pos_[r.argmin_global]), r.argmin_global};
}

// -----------------------------------------------------------------------------
// query_special()
// -----------------------------------------------------------------------------

Grid2D::SpecialResult Grid2D::query_special(size_t sp_rev, size_t ep_rev, size_t plen) const {
    if (wm_.size() == 0) return {0, SIZE_MAX, 0};

    // Para ocurrencias especiales, Y_k esta en [sp_rev, ep_rev] y X_k es libre.
    const size_t vlb = rank_rev_(sp_rev);
    if (ep_rev + 1 > bv_rev_.size()) return {0, SIZE_MAX, 0};
    const size_t vrb = rank_rev_(ep_rev + 1);
    if (vrb == 0 || vlb >= vrb) return {0, SIZE_MAX, 0};

    const size_t lb = 0;
    const size_t rb = wm_.size() - 1;

    const auto pts = wm_.range_search_2d(lb, rb, vlb, vrb - 1).second;

    size_t count = 0;
    size_t occ_min = SIZE_MAX;
    size_t occ_max = 0;

    for (const auto& [wt_idx, y_rel] : pts) {
        if (phrase_total_len_[wt_idx] >= plen) {
            count++;
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
    size_t written = 0;
    written += sdsl::serialize(static_cast<const sdsl::wm_int<>&>(wm_), out);
    written += sdsl::serialize(bv_fwd_,          out);
    written += sdsl::serialize(bv_rev_,          out);
    written += sdsl::serialize(text_pos_,        out);
    written += sdsl::serialize(phrase_total_len_,out);
    written += wm_min_rmq_.serialize(out);
    return written;
}

void Grid2D::load(std::istream& in) {
    sdsl::load(static_cast<sdsl::wm_int<>&>(wm_), in);
    sdsl::load(bv_fwd_,           in);
    sdsl::load(bv_rev_,           in);
    sdsl::load(text_pos_,         in);
    sdsl::load(phrase_total_len_, in);
    wm_min_rmq_.load(in);
    // rank_fwd_ y rank_rev_ guardan puntero a bv_*_: reconstruir tras cargar.
    sdsl::util::init_support(rank_fwd_, &bv_fwd_);
    sdsl::util::init_support(rank_rev_, &bv_rev_);
}

Grid2D::SizeBreakdown Grid2D::size_breakdown() const {
    SizeBreakdown bd{};
    bd.wm              = sdsl::size_in_bytes(wm_);
    bd.bv_fwd          = sdsl::size_in_bytes(bv_fwd_);
    bd.bv_rev          = sdsl::size_in_bytes(bv_rev_);
    bd.rank_fwd        = sdsl::size_in_bytes(rank_fwd_);
    bd.rank_rev        = sdsl::size_in_bytes(rank_rev_);
    bd.text_pos        = sdsl::size_in_bytes(text_pos_);
    bd.phrase_total_len= sdsl::size_in_bytes(phrase_total_len_);
    bd.wm_min_rmq      = wm_min_rmq_.size_in_bytes();
    return bd;
}

}  // namespace lz77tax
