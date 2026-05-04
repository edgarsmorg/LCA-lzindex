#include "grid.hpp"

#include <algorithm>
#include <cassert>
#include <stdexcept>

#include <sdsl/construct.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/bits.hpp>

namespace lz77tax {

// ─────────────────────────────────────────────────────────────────────────────
// Núcleo compartido
// ─────────────────────────────────────────────────────────────────────────────

void Grid2D::build_from_coords(
        const std::vector<std::pair<size_t, size_t>>& coords,
        const std::vector<size_t>& boundaries,
        const std::vector<size_t>& phrase_lens,
        size_t n) {

    const size_t z1 = coords.size();  // z-1 puntos
    assert(boundaries.size() == z1 && "coords y boundaries deben tener el mismo tamaño");
    assert(phrase_lens.size() == z1 && "coords y phrase_lens deben tener el mismo tamaño");

    if (z1 == 0) {
        // Texto de 1 frase: grilla vacía, estructuras en estado default.
        return;
    }

    // ── 1. Bitvector forward: un 1 en cada X_k ───────────────────────────────
    {
        sdsl::bit_vector bfwd(n, 0);
        for (const auto& [x, y] : coords) {
            assert(x < n && "X_k fuera del rango del SA");
            bfwd[x] = 1;
        }
        bv_fwd_ = sdsl::sd_vector<>(bfwd);
    }
    rank_fwd_ = sdsl::sd_vector<>::rank_1_type(&bv_fwd_);

    // ── 2. Bitvector reverse: un 1 en cada Y_k ───────────────────────────────
    {
        sdsl::bit_vector brev(n, 0);
        for (const auto& [x, y] : coords) {
            assert(y < n && "Y_k fuera del rango del rev-SA");
            brev[y] = 1;
        }
        bv_rev_ = sdsl::sd_vector<>(brev);
    }
    rank_rev_ = sdsl::sd_vector<>::rank_1_type(&bv_rev_);

    // ── 3. Secuencia R[0..z-2] y array text_pos_ ─────────────────────────────
    // R[rank_fwd(X_k)] = rank_rev(Y_k)
    // Ambos ranks son 0-indexed y forman permutaciones de {0..z-2} ya que
    // los X_k son distintos (ISA es permutación) e idem para Y_k.
    // text_pos_[j] = boundaries[k] = start_{k+1} para el punto con índice WT j.
    const uint8_t width = sdsl::bits::hi(z1) + 1;  // floor(log2(z1)) + 1: bits para [0..z1]
    sdsl::int_vector<> R(z1, 0, width);

    // text_pos_ compact: valores en [0, n-1], width = ceil(log2(n)) bits
    const uint8_t tp_w = n > 1 ? static_cast<uint8_t>(sdsl::bits::hi(n - 1) + 1) : 1;
    text_pos_ = sdsl::int_vector<>(z1, 0, tp_w);
    phrase_total_len_ = sdsl::int_vector<>(z1, 0, tp_w);

    for (size_t k = 0; k < z1; ++k) {
        const size_t x = coords[k].first;
        const size_t y = coords[k].second;
        // rank_fwd_(x) = #{X' < x entre todos los X_k} = índice en wt_
        const size_t j     = rank_fwd_(x);
        // rank_rev_(y) = #{Y' < y entre todos los Y_k} = valor en wt_
        const size_t y_rel = rank_rev_(y);
        assert(j < z1 && "índice wt_ fuera de rango");
        assert(y_rel < z1 && "y_rel fuera de rango");
        R[j]                 = y_rel;
        text_pos_[j]         = boundaries[k];
        phrase_total_len_[j] = phrase_lens[k];
    }

    // shadow plain para WtMinRmq (pasa en cada consulta)
    text_pos_plain_.resize(z1);
    for (size_t j = 0; j < z1; ++j)
        text_pos_plain_[j] = static_cast<size_t>(text_pos_[j]);

    // ── 4. Construir wt_int sobre R ───────────────────────────────────────────
    sdsl::construct_im(wt_, R);

    // ── 5. Construir WtMinRmq sobre (ys=R, text_pos_plain_, sigma=z1) ─────────
    {
        std::vector<size_t> ys(z1);
        for (size_t j = 0; j < z1; ++j)
            ys[j] = static_cast<size_t>(R[j]);
        wt_min_rmq_.build(ys, text_pos_plain_, z1);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// build() desde texto crudo (small-scale / tests)
// ─────────────────────────────────────────────────────────────────────────────

void Grid2D::build(const LZ77Parsing& phrases, const std::string& text) {
    const size_t n = text.size();
    const size_t z = phrases.size();

    if (z <= 1) return;  // 0 o 1 frase → grilla vacía

    // ── SA y ISA del texto ────────────────────────────────────────────────────
    sdsl::int_vector<> sa(n);
    sdsl::algorithm::calculate_sa(
            reinterpret_cast<const unsigned char*>(text.data()), n, sa);

    sdsl::int_vector<> isa(n);
    for (size_t i = 0; i < n; ++i) isa[sa[i]] = i;

    // ── SA y ISA del texto reverso ────────────────────────────────────────────
    // IMPORTANTE: no revertimos el centinela '\0' al inicio — calculate_sa
    // requiere que el carácter más pequeño esté al FINAL.
    // text = "abracadabra\0"  (n=12, '\0' en posición n-1)
    // Revertimos solo text[0..n-2] y añadimos '\0' al final:
    //   text_rev_sa = "arbadacarba\0"
    // Así text_rev_sa[j] = text[n-2-j] para j=0..n-2, y '\0' en j=n-1.
    // text termina con el centinela '\0' en posición n-1 (= text.rbegin()).
    // Saltamos ese centinela y revertimos solo text[0..n-2], luego añadimos '\0' al final.
    // Resultado: text_rev_sa_s[j] = text[n-2-j] para j=0..n-2, text_rev_sa_s[n-1]='\0'.
    const std::string text_rev_sa(text.rbegin() + 1, text.rend());  // salta centinela final
    const std::string text_rev_sa_s = text_rev_sa + '\0';

    sdsl::int_vector<> sa_rev(n);
    sdsl::algorithm::calculate_sa(
            reinterpret_cast<const unsigned char*>(text_rev_sa_s.data()), n, sa_rev);

    sdsl::int_vector<> isa_rev(n);
    for (size_t i = 0; i < n; ++i) isa_rev[sa_rev[i]] = i;

    // ── Coordenadas de los z-1 puntos ─────────────────────────────────────────
    // Para k = 0..z-2:
    //   end_k    = phrases[k].start_pos + phrases[k].length
    //   start_k1 = phrases[k+1].start_pos  (== end_k + 1)
    //   X_k = ISA[ start_k1 ]
    //   Y_k = ISA_rev[ n-2 - end_k ]
    //
    // La fórmula Y_k = n-2-end_k viene de que text_rev_sa[j] = text[n-2-j],
    // por lo que para buscar el carácter text[end_k] en el SA reverso usamos
    // el índice j = n-2-end_k.
    std::vector<std::pair<size_t, size_t>> coords(z - 1);
    std::vector<size_t> boundaries(z - 1);
    std::vector<size_t> phrase_lens(z - 1);
    for (size_t k = 0; k < z - 1; ++k) {
        const size_t end_k    = phrases[k].start_pos + phrases[k].length;
        const size_t start_k1 = phrases[k + 1].start_pos;

        assert(start_k1 == end_k + 1 &&
               "gap entre frases consecutivas: start_k1 != end_k + 1");
        assert(start_k1 < n    && "start_k1 fuera del texto");
        assert(end_k < n - 1   && "end_k apunta al centinela — la última frase no debe incluirse");
        assert(n >= 2 + end_k  && "n-2-end_k subflow");

        coords[k]      = { isa[start_k1], isa_rev[n - 2 - end_k] };
        boundaries[k]  = start_k1;
        phrase_lens[k] = phrases[k].length + 1; // total length including trailing char
    }

    build_from_coords(coords, boundaries, phrase_lens, n);
}

// ─────────────────────────────────────────────────────────────────────────────
// build() desde FM-index (producción, GB+)
// ─────────────────────────────────────────────────────────────────────────────

void Grid2D::build(const LZ77Parsing& phrases,
                   const sdsl::csa_wt<>& csa_fwd,
                   const sdsl::csa_wt<>& csa_rev) {
    const size_t n = csa_fwd.size();
    const size_t z = phrases.size();

    if (z <= 1) return;

    std::vector<std::pair<size_t, size_t>> coords(z - 1);
    std::vector<size_t> boundaries(z - 1);
    std::vector<size_t> phrase_lens(z - 1);
    for (size_t k = 0; k < z - 1; ++k) {
        const size_t end_k    = phrases[k].start_pos + phrases[k].length;
        const size_t start_k1 = phrases[k + 1].start_pos;

        // csa_fwd.isa[i] resuelve ISA[i] en O(t_psi) sin materializar el array.
        // Y_k = csa_rev.isa[n-2-end_k]: csa_rev debe estar construido sobre
        // reverse(text_raw) + '\0' (no reverse(text_with_sentinel)).
        coords[k]      = { csa_fwd.isa[start_k1], csa_rev.isa[n - 2 - end_k] };
        boundaries[k]  = start_k1;
        phrase_lens[k] = phrases[k].length + 1;
    }

    build_from_coords(coords, boundaries, phrase_lens, n);
}

// ─────────────────────────────────────────────────────────────────────────────
// query()
// ─────────────────────────────────────────────────────────────────────────────

Grid2D::RangeResult Grid2D::query(size_t sp_right, size_t ep_right,
                                   size_t sp_left,  size_t ep_left) const {
    if (wt_.size() == 0) return {0, {}};

    // Convertir rangos BWT globales [sp, ep] a rangos relativos en wt_.
    // rank_1_type(i) = #{j < i : bv[j]=1}  (exclusivo en i, 0-indexed)
    //   → frases con X ∈ [sp_right, ep_right] tienen índice WT en [lb, rb]
    const size_t lb  = rank_fwd_(sp_right);
    // Para inputs válidos ep_right < n = bv_fwd_.size(), así que esta guarda
    // no se activa en uso normal; es defensiva ante llamadas con rangos inválidos.
    if (ep_right + 1 > bv_fwd_.size()) return {0, {}};
    const size_t rb  = rank_fwd_(ep_right + 1);
    if (rb == 0 || lb >= rb) return {0, {}};
    const size_t rb_idx = rb - 1;

    const size_t vlb = rank_rev_(sp_left);
    if (ep_left + 1 > bv_rev_.size()) return {0, {}};
    const size_t vrb = rank_rev_(ep_left + 1);
    if (vrb == 0 || vlb >= vrb) return {0, {}};
    const size_t vrb_val = vrb - 1;

    if (lb > rb_idx || vlb > vrb_val) return {0, {}};

    return wt_.range_search_2d(lb, rb_idx, vlb, vrb_val);
}

// ─────────────────────────────────────────────────────────────────────────────
// query_extremal()
// ─────────────────────────────────────────────────────────────────────────────

Grid2D::ExtremalResult Grid2D::query_extremal(size_t sp_right, size_t ep_right,
                                               size_t sp_left,  size_t ep_left) const {
    if (wt_.size() == 0 || text_pos_.empty()) return {0, SIZE_MAX, 0};

    const auto [cnt, pts] = query(sp_right, ep_right, sp_left, ep_left);
    if (cnt == 0) return {0, SIZE_MAX, 0};

    size_t bmin = SIZE_MAX, bmax = 0;
    for (const auto& [wt_idx, y_rel] : pts) {
        // sdsl devuelve {position, value} (ver wt_int.hpp emplace_back(i-1, path)):
        //   wt_idx = índice original en R = rank_fwd_(X_k)
        //   y_rel  = valor en R[wt_idx]  = rank_rev_(Y_k)
        const size_t bp = text_pos_[wt_idx];
        if (bp < bmin) bmin = bp;
        if (bp > bmax) bmax = bp;
    }
    return {cnt, bmin, bmax};
}

// ─────────────────────────────────────────────────────────────────────────────
// query_min_2d()
// ─────────────────────────────────────────────────────────────────────────────

Grid2D::MinResult Grid2D::query_min_2d(size_t sp_right, size_t ep_right,
                                        size_t sp_left,  size_t ep_left) const {
    if (wt_min_rmq_.size() == 0) return {0, SIZE_MAX};

    const size_t lb = rank_fwd_(sp_right);
    if (ep_right + 1 > bv_fwd_.size()) return {0, SIZE_MAX};
    const size_t rb = rank_fwd_(ep_right + 1);
    if (rb == 0 || lb >= rb) return {0, SIZE_MAX};

    const size_t vlb = rank_rev_(sp_left);
    if (ep_left + 1 > bv_rev_.size()) return {0, SIZE_MAX};
    const size_t vrb = rank_rev_(ep_left + 1);
    if (vrb == 0 || vlb >= vrb) return {0, SIZE_MAX};

    const auto r = wt_min_rmq_.range_argmin_2d(lb, rb - 1, vlb, vrb - 1,
                                                text_pos_plain_);
    if (r.count == 0) return {0, SIZE_MAX};
    return {r.count, text_pos_plain_[r.argmin_global]};
}

// ─────────────────────────────────────────────────────────────────────────────
// query_special()
// ─────────────────────────────────────────────────────────────────────────────

Grid2D::SpecialResult Grid2D::query_special(size_t sp_rev, size_t ep_rev, size_t plen) const {
    if (wt_.size() == 0) return {0, SIZE_MAX, 0};

    // Para occurrences special, iteramos sobre los Y_k en [sp_rev, ep_rev].
    // X_k no está restringido (puede ser cualquier parte del texto).
    const size_t lb  = 0;
    const size_t rb  = wt_.size() - 1;

    const size_t vlb = rank_rev_(sp_rev);
    if (ep_rev + 1 > bv_rev_.size()) return {0, SIZE_MAX, 0};
    const size_t vrb = rank_rev_(ep_rev + 1);
    if (vrb == 0 || vlb >= vrb) return {0, SIZE_MAX, 0};
    const size_t vrb_val = vrb - 1;

    if (vlb > vrb_val) return {0, SIZE_MAX, 0};

    const auto pts = wt_.range_search_2d(lb, rb, vlb, vrb_val).second;

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

}  // namespace lz77tax
