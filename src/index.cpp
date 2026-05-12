#include "index.hpp"

#include <algorithm>
#include <cassert>
#include <climits>
#include <string>

#include <sdsl/suffix_arrays.hpp>

namespace lz77tax {

// ─────────────────────────────────────────────────────────────────────────────
// build() — versión small-scale / tests
// ─────────────────────────────────────────────────────────────────────────────

void LZ77Index::build(const std::string& text) {
    // 1. Centinela: la grid y el parser lo necesitan al final.
    const std::string text_s = text + '\0';

    // 2. Parsear sobre texto con centinela.
    phrases_ = LZ77Parser().parse(text_s);

    // 3. Construir grid (text_s incluye '\0' → SA tamaño n = text.size()+1).
    grid_.build(phrases_, text_s);

    // 4. CSA forward sobre texto SIN centinela.
    //    sdsl::construct_im rechaza '\0' en el input y lo añade internamente,
    //    por lo que el SA resultante es idéntico al usado por la grid.
    sdsl::construct_im(csa_fwd_, text, 1);  // 1 = byte alphabet

    // 5. CSA reverso sobre reverse(text) SIN centinela.
    //    Consistente con grid.build(): ambos invierten text completo (sin centinela).
    //    sdsl::construct_im añade '\0' internamente → SA de reverse(text)+'\0'.
    const std::string text_rev(text.rbegin(), text.rend());
    sdsl::construct_im(csa_rev_, text_rev, 1);
}

// ─────────────────────────────────────────────────────────────────────────────
// build() — versión producción (GB+)
// ─────────────────────────────────────────────────────────────────────────────

void LZ77Index::build(const LZ77Parsing& phrases,
                      const sdsl::csa_wt<>& csa_fwd,
                      const sdsl::csa_wt<>& csa_rev) {
    phrases_ = phrases;
    csa_fwd_ = csa_fwd;
    csa_rev_ = csa_rev;
    grid_.build(phrases_, csa_fwd_, csa_rev_);
}

// ─────────────────────────────────────────────────────────────────────────────
// count()
// ─────────────────────────────────────────────────────────────────────────────

size_t LZ77Index::count(const std::string& pattern) const {
    const size_t m = pattern.size();

    // Un patrón de longitud 0 o 1 no puede ser primario (crossing ni end-aligned).
    if (m < 2 || grid_.point_count() == 0) return 0;

    const size_t n = csa_fwd_.size();
    size_t total = 0;

    // Precomputar reverse(P) una vez — usado en el caso especial y reutilizable.
    const std::string rev_p(pattern.rbegin(), pattern.rend());

    // ── Special: P termina en el trailing_char de alguna frase k ─────────────
    // Equivalente al split degenerado i=m: P[0..m-1] | vacío.
    {
        size_t sp_rev = 0, ep_rev = n - 1;
        sdsl::backward_search(csa_rev_, 0, n - 1,
                              rev_p.begin(), rev_p.end(),
                              sp_rev, ep_rev);
        if (sp_rev <= ep_rev) {
            const auto sp = grid_.query_special(sp_rev, ep_rev, m);
            total += sp.count;
        }
    }

    // ── Crossings: P[0..i-1] | P[i..m-1] con i = 1..m-1 ────────────────────
    for (size_t i = 1; i < m; ++i) {
        // ── Lado derecho: P[i..m-1] en SA forward ────────────────────────────
        size_t sp_r = 0, ep_r = n - 1;
        sdsl::backward_search(csa_fwd_, 0, n - 1,
                              pattern.begin() + i, pattern.end(),
                              sp_r, ep_r);
        if (sp_r > ep_r) continue;  // P[i..m-1] no aparece en T

        // ── Lado izquierdo: reverse(P[0..i-1]) en SA reverso ─────────────────
        // left_rev = reverse(P[0..i-1]) = P[i-1], P[i-2], ..., P[0]
        const std::string left_rev(pattern.rbegin() + (m - i), pattern.rend());
        size_t sp_l = 0, ep_l = n - 1;
        sdsl::backward_search(csa_rev_, 0, n - 1,
                              left_rev.begin(), left_rev.end(),
                              sp_l, ep_l);
        if (sp_l > ep_l) continue;  // P[0..i-1] no aparece en T

        // ── Consulta 2D en la grilla ──────────────────────────────────────────
        const auto [cnt, pts] = grid_.query(sp_r, ep_r, sp_l, ep_l);
        total += cnt;
    }

    return total;
}

// ─────────────────────────────────────────────────────────────────────────────
// locate_extremal()
// ─────────────────────────────────────────────────────────────────────────────

std::pair<size_t, size_t> LZ77Index::locate_extremal(const std::string& pattern) const {
    const size_t m = pattern.size();
    if (m < 2 || grid_.point_count() == 0) return {SIZE_MAX, 0};

    const size_t n = csa_fwd_.size();
    size_t pos_min = SIZE_MAX, pos_max = 0;
    bool found = false;

    const std::string rev_p(pattern.rbegin(), pattern.rend());

    // ── Special: P termina en el trailing_char de alguna frase k ─────────────
    {
        size_t sp_rev = 0, ep_rev = n - 1;
        sdsl::backward_search(csa_rev_, 0, n - 1,
                              rev_p.begin(), rev_p.end(),
                              sp_rev, ep_rev);
        if (sp_rev <= ep_rev) {
            const auto sp = grid_.query_special(sp_rev, ep_rev, m);
            if (sp.count > 0) {
                if (sp.occ_min_pos < pos_min) pos_min = sp.occ_min_pos;
                if (sp.occ_max_pos > pos_max) pos_max = sp.occ_max_pos;
                found = true;
            }
        }
    }

    // ── Crossings: P[0..i-1] | P[i..m-1] con i = 1..m-1 ────────────────────
    for (size_t i = 1; i < m; ++i) {
        // ── Lado derecho: P[i..m-1] en SA forward ────────────────────────────
        size_t sp_r = 0, ep_r = n - 1;
        sdsl::backward_search(csa_fwd_, 0, n - 1,
                              pattern.begin() + i, pattern.end(),
                              sp_r, ep_r);
        if (sp_r > ep_r) continue;

        // ── Lado izquierdo: reverse(P[0..i-1]) en SA reverso ─────────────────
        const std::string left_rev(pattern.rbegin() + (m - i), pattern.rend());
        size_t sp_l = 0, ep_l = n - 1;
        sdsl::backward_search(csa_rev_, 0, n - 1,
                              left_rev.begin(), left_rev.end(),
                              sp_l, ep_l);
        if (sp_l > ep_l) continue;

        // ── query_min_2d + query_max_2d: O(log² z) sin enumerar puntos ────────
        const auto rmin = grid_.query_min_2d(sp_r, ep_r, sp_l, ep_l);
        if (rmin.count == 0) continue;

        const auto rmax = grid_.query_max_2d(sp_r, ep_r, sp_l, ep_l);

        if (rmin.boundary_min >= i) {
            const size_t p_min = rmin.boundary_min - i;
            if (p_min < pos_min) pos_min = p_min;
        }
        if (rmax.count > 0 && rmax.boundary_max >= i) {
            const size_t p_max = rmax.boundary_max - i;
            if (p_max > pos_max) pos_max = p_max;
        }
        found = true;
    }

    return found ? std::make_pair(pos_min, pos_max)
                 : std::make_pair(SIZE_MAX, size_t(0));
}

// ─────────────────────────────────────────────────────────────────────────────
// locate_min()
// ─────────────────────────────────────────────────────────────────────────────

size_t LZ77Index::locate_min(const std::string& pattern) const {
    const size_t m = pattern.size();
    if (m < 2 || grid_.point_count() == 0) return SIZE_MAX;

    const size_t n = csa_fwd_.size();
    size_t pos_min = SIZE_MAX;

    const std::string rev_p(pattern.rbegin(), pattern.rend());

    // ── Special: P termina en el trailing_char de alguna frase k ─────────────
    {
        size_t sp_rev = 0, ep_rev = n - 1;
        sdsl::backward_search(csa_rev_, 0, n - 1,
                              rev_p.begin(), rev_p.end(),
                              sp_rev, ep_rev);
        if (sp_rev <= ep_rev) {
            const auto sp = grid_.query_special(sp_rev, ep_rev, m);
            if (sp.count > 0 && sp.occ_min_pos < pos_min)
                pos_min = sp.occ_min_pos;
        }
    }

    // ── Crossings: P[0..i-1] | P[i..m-1] con i = 1..m-1 ────────────────────
    for (size_t i = 1; i < m; ++i) {
        size_t sp_r = 0, ep_r = n - 1;
        sdsl::backward_search(csa_fwd_, 0, n - 1,
                              pattern.begin() + i, pattern.end(),
                              sp_r, ep_r);
        if (sp_r > ep_r) continue;

        const std::string left_rev(pattern.rbegin() + (m - i), pattern.rend());
        size_t sp_l = 0, ep_l = n - 1;
        sdsl::backward_search(csa_rev_, 0, n - 1,
                              left_rev.begin(), left_rev.end(),
                              sp_l, ep_l);
        if (sp_l > ep_l) continue;

        const auto r = grid_.query_min_2d(sp_r, ep_r, sp_l, ep_l);
        if (r.count == 0) continue;

        if (r.boundary_min < i) continue;  // invariante: occ_pos = boundary - split >= 0
        const size_t p = r.boundary_min - i;
        if (p < pos_min) pos_min = p;
    }

    return pos_min;
}

// ─────────────────────────────────────────────────────────────────────────────
// locate_max()
// ─────────────────────────────────────────────────────────────────────────────

size_t LZ77Index::locate_max(const std::string& pattern) const {
    const size_t m = pattern.size();
    if (m < 2 || grid_.point_count() == 0) return SIZE_MAX;

    const size_t n = csa_fwd_.size();
    size_t pos_max = 0;
    bool found = false;

    const std::string rev_p(pattern.rbegin(), pattern.rend());

    // ── Special: P termina en el trailing_char de alguna frase k ─────────────
    {
        size_t sp_rev = 0, ep_rev = n - 1;
        sdsl::backward_search(csa_rev_, 0, n - 1,
                              rev_p.begin(), rev_p.end(),
                              sp_rev, ep_rev);
        if (sp_rev <= ep_rev) {
            const auto sp = grid_.query_special(sp_rev, ep_rev, m);
            if (sp.count > 0 && sp.occ_max_pos > pos_max) {
                pos_max = sp.occ_max_pos;
                found = true;
            }
        }
    }

    // ── Crossings: P[0..i-1] | P[i..m-1] con i = 1..m-1 ────────────────────
    for (size_t i = 1; i < m; ++i) {
        size_t sp_r = 0, ep_r = n - 1;
        sdsl::backward_search(csa_fwd_, 0, n - 1,
                              pattern.begin() + i, pattern.end(),
                              sp_r, ep_r);
        if (sp_r > ep_r) continue;

        const std::string left_rev(pattern.rbegin() + (m - i), pattern.rend());
        size_t sp_l = 0, ep_l = n - 1;
        sdsl::backward_search(csa_rev_, 0, n - 1,
                              left_rev.begin(), left_rev.end(),
                              sp_l, ep_l);
        if (sp_l > ep_l) continue;

        const auto r = grid_.query_max_2d(sp_r, ep_r, sp_l, ep_l);
        if (r.count == 0 || r.boundary_max < i) continue;

        const size_t p = r.boundary_max - i;
        if (!found || p > pos_max) {
            pos_max = p;
            found = true;
        }
    }

    return found ? pos_max : SIZE_MAX;
}

}  // namespace lz77tax
