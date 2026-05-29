#include "index.hpp"

#include <algorithm>
#include <cassert>
#include <climits>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>

#include <sdsl/csa_wt.hpp>
#include <sdsl/construct_isa.hpp>
#include <sdsl/construct.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>

#include <sr-index/r_csa.h>
#include <sr-index/construct.h>
#include <sr-index/config.h>
#include <sr-index/index_base.h>

// ─────────────────────────────────────────────────────────────────────────────
// CountOnlyRCSA — subclase de sri::RCSA<> que solo carga/serializa las
// estructuras necesarias para Count() (Alphabet + PsiCoreRLE).
// Los samples/marks/mark-to-sample del Locate se omiten completamente:
//   • No se cargan en loadAllItems() → no consumen RAM.
//   • No se escriben en serialize()  → serialize_bytes() reporta solo
//     los ~25 MB reales de Alphabet + PsiRLE por dirección.
// Locate() queda deshabilitado (stubs); la clase es solo para Count().
//
// Nota: la clase está en file scope (no en anonymous namespace) para que
// sdsl::has_serialize<CountOnlyRCSA>::value sea true y sdsl::size_in_bytes
// funcione correctamente a través del despacho virtual de serialize().
// ─────────────────────────────────────────────────────────────────────────────

class CountOnlyRCSA : public sri::RCSA<> {
    using Base    = sri::RCSA<>;
    using StrType = typename Alphabet::string_type;
 public:
    // serialize es public en la base: override debe también ser public para que
    // sdsl::has_serialize lo detecte y rcsa_bytes() lo pueda llamar externamente.
    size_type serialize(std::ostream& out, sdsl::structure_tree_node* v,
                        const std::string& name) const override {
        auto child = sdsl::structure_tree::add_child(v, name, "CountOnlyRCSA");
        size_type written = 0;
        written += this->template serializeItem<Alphabet>(
            key(ItemKey::ALPHABET), out, child, "alphabet");
        written += this->template serializeItem<sri::PsiCoreRLE<>>(
            key(ItemKey::NAVIGATE), out, child, "psi");
        return written;
    }

 protected:
    void loadAllItems(TSource& t_source) override {
        this->template loadItem<Alphabet>(key(ItemKey::ALPHABET), t_source);
        this->template loadItem<sri::PsiCoreRLE<>>(key(ItemKey::NAVIGATE), t_source, true);
        // Skip: TSample, TBvMark+rank/select, TMarkToSampleIdx — Count no los usa.
    }

    void constructIndex(TSource& t_source) override {
        auto lf      = this->constructLF(t_source);   // también setea this->n_
        auto sym     = this->constructGetSymbol(t_source);
        auto is_mt   = this->constructIsRangeEmpty();
        auto no_init = [](std::size_t) noexcept -> int { return 0; };
        auto no_upd  = [](auto&&...) noexcept {};
        auto no_cmp  = [](auto&&...) noexcept {};
        auto mk_rng  = [](std::size_t n) noexcept -> Range { return {0, n}; };

        using LF_t  = decltype(lf);
        using Sym_t = decltype(sym);
        using I_t   = decltype(no_init);
        using U_t   = decltype(no_upd);
        using C_t   = decltype(no_cmp);
        using R_t   = decltype(mk_rng);
        using E_t   = decltype(is_mt);

        this->index_ = std::make_shared<
            sri::RIndexBase<StrType, LF_t, U_t, C_t, I_t, Sym_t, R_t, E_t>
        >(StrType{}, lf, no_upd, no_cmp, this->n_, no_init, sym, mk_rng, is_mt);
    }
};

namespace lz77tax {

// PIMPL: mantiene los dos RCSA lejos del header (evita ODR por funciones no-inline del sr-index)
struct LZ77Index::RCSAImpl {
    CountOnlyRCSA fwd;
    CountOnlyRCSA rev;
};

LZ77Index::LZ77Index() = default;
LZ77Index::~LZ77Index() = default;

// ─────────────────────────────────────────────────────────────────────────────
// Helpers internos
// ─────────────────────────────────────────────────────────────────────────────

namespace {

// Escribe `text` a un archivo temporal en `tmp_dir` y devuelve su path.
std::filesystem::path write_temp_text(const std::string& text,
                                      const std::filesystem::path& tmp_dir,
                                      const std::string& suffix) {
    auto path = tmp_dir / ("lz77_" + suffix + ".txt");
    std::ofstream f(path, std::ios::binary);
    f.write(text.data(), text.size());
    return path;
}

// Construye un CountOnlyRCSA para `text_path`, extrae ISA desde el cache de sdsl y lo devuelve.
// config tiene delete_files=true → los archivos intermedios se borran al destruir config.
sdsl::int_vector<> build_rcsa_and_get_isa(CountOnlyRCSA& rcsa,
                                           const std::filesystem::path& text_path,
                                           const std::filesystem::path& cache_dir) {
    sri::Config config(text_path, cache_dir, sri::SAAlgo::SDSL_LIBDIVSUFSORT,
                       /*delete_files=*/true);
    sri::construct(rcsa, text_path.string(), config);

    // ISA no se construye en sri::construct — hacerlo aquí.
    sdsl::construct_isa(config);
    sdsl::int_vector<> isa;
    sdsl::load_from_cache(isa, sdsl::conf::KEY_ISA, config);
    return isa;
}

std::filesystem::path make_tmp_dir(const std::string& prefix) {
    std::mt19937_64 rng(std::random_device{}());
    auto tmp = std::filesystem::temp_directory_path()
               / (prefix + "_" + std::to_string(rng()));
    std::filesystem::create_directories(tmp);
    return tmp;
}

}  // namespace

// ─────────────────────────────────────────────────────────────────────────────
// build() — versión desde texto en memoria (small-scale y benchmark)
// ─────────────────────────────────────────────────────────────────────────────

void LZ77Index::build(const std::string& text) {
    const std::string text_s = text + '\0';
    n_ = text_s.size();
    LZ77Parsing phrases = LZ77Parser().parse(text_s);
    z_ = phrases.size();

    alphabet_.reset();
    for (unsigned char c : text) alphabet_.set(c);

    rcsa_ = std::make_unique<RCSAImpl>();

    // Directorios temporales separados para fwd y rev (evitar colisiones de keys).
    const auto tmp_fwd = make_tmp_dir("lz77fwd");
    const auto tmp_rev = make_tmp_dir("lz77rev");

    // Texto SIN centinela para sri::construct (lo añade internamente via sdsl).
    const auto path_fwd = write_temp_text(text, tmp_fwd, "text");
    const auto isa_fwd  = build_rcsa_and_get_isa(rcsa_->fwd, path_fwd, tmp_fwd);
    std::filesystem::remove(path_fwd);
    std::filesystem::remove_all(tmp_fwd);

    const std::string text_rev(text.rbegin(), text.rend());
    const auto path_rev = write_temp_text(text_rev, tmp_rev, "text");
    const auto isa_rev  = build_rcsa_and_get_isa(rcsa_->rev, path_rev, tmp_rev);
    std::filesystem::remove(path_rev);
    std::filesystem::remove_all(tmp_rev);

    grid_.build(phrases, isa_fwd, isa_rev, n_);

    // Liberar la RAM del vector de frases: Grid2D ya extrajo todo lo que necesita.
    phrases = LZ77Parsing{};
}

namespace {
// Cuenta bytes que serialize() escribe llamando directamente al método virtual.
// Evita sdsl::size_in_bytes que requiere has_serialize<T> vía free function.
size_t rcsa_bytes(const CountOnlyRCSA& rcsa) {
    struct NullBuf : std::streambuf {
        std::streamsize xsputn(const char*, std::streamsize n) override { return n; }
        int overflow(int c) override { return c; }
    } buf;
    std::ostream ns(&buf);
    return rcsa.serialize(ns, nullptr, "");
}
}  // namespace

size_t LZ77Index::csa_fwd_bytes() const {
    return rcsa_ ? rcsa_bytes(rcsa_->fwd) : 0;
}

size_t LZ77Index::csa_rev_bytes() const {
    return rcsa_ ? rcsa_bytes(rcsa_->rev) : 0;
}


// ─────────────────────────────────────────────────────────────────────────────
// count()
// ─────────────────────────────────────────────────────────────────────────────

namespace {
// Convierte rango half-open [start, stop) de CountOnlyRCSA::Count a closed [sp, ep].
// Precondición: sub solo contiene caracteres en el alfabeto del índice.
bool rcsa_search(const CountOnlyRCSA& rcsa, const std::string& sub,
                 size_t& sp, size_t& ep) {
    auto [start, stop] = rcsa.Count(sub);
    if (start >= stop) return false;
    sp = start;
    ep = stop - 1;
    return true;
}

// Retorna true si todos los chars de s están en el alfabeto indexado.
bool in_alphabet(const std::string& s, const std::bitset<256>& alpha) {
    for (unsigned char c : s)
        if (!alpha.test(c)) return false;
    return true;
}
}  // namespace

size_t LZ77Index::count(const std::string& pattern) const {
    const size_t m = pattern.size();

    if (m < 2 || grid_.point_count() == 0 || !rcsa_) return 0;
    if (!in_alphabet(pattern, alphabet_)) return 0;

    size_t total = 0;
    const std::string rev_p(pattern.rbegin(), pattern.rend());

    // ── Special: P termina en el trailing_char de alguna frase k ─────────────
    {
        size_t sp_rev, ep_rev;
        if (rcsa_search(rcsa_->rev, rev_p, sp_rev, ep_rev)) {
            const auto sp = grid_.query_special(sp_rev, ep_rev, m);
            total += sp.count;
        }
    }

    // ── Crossings: P[0..i-1] | P[i..m-1] con i = 1..m-1 ────────────────────
    for (size_t i = 1; i < m; ++i) {
        size_t sp_r, ep_r;
        if (!rcsa_search(rcsa_->fwd, pattern.substr(i), sp_r, ep_r)) continue;

        const std::string left_rev(pattern.rbegin() + (m - i), pattern.rend());
        size_t sp_l, ep_l;
        if (!rcsa_search(rcsa_->rev, left_rev, sp_l, ep_l)) continue;

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
    if (!in_alphabet(pattern, alphabet_)) return {SIZE_MAX, 0};

    size_t pos_min = SIZE_MAX, pos_max = 0;
    bool found = false;

    const std::string rev_p(pattern.rbegin(), pattern.rend());

    // ── Special ───────────────────────────────────────────────────────────────
    {
        size_t sp_rev, ep_rev;
        if (rcsa_search(rcsa_->rev, rev_p, sp_rev, ep_rev)) {
            const auto sp = grid_.query_special(sp_rev, ep_rev, m);
            if (sp.count > 0) {
                if (sp.occ_min_pos < pos_min) pos_min = sp.occ_min_pos;
                if (sp.occ_max_pos > pos_max) pos_max = sp.occ_max_pos;
                found = true;
            }
        }
    }

    // ── Crossings ─────────────────────────────────────────────────────────────
    for (size_t i = 1; i < m; ++i) {
        size_t sp_r, ep_r;
        if (!rcsa_search(rcsa_->fwd, pattern.substr(i), sp_r, ep_r)) continue;

        const std::string left_rev(pattern.rbegin() + (m - i), pattern.rend());
        size_t sp_l, ep_l;
        if (!rcsa_search(rcsa_->rev, left_rev, sp_l, ep_l)) continue;

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
    if (!in_alphabet(pattern, alphabet_)) return SIZE_MAX;

    size_t pos_min = SIZE_MAX;
    const std::string rev_p(pattern.rbegin(), pattern.rend());

    // ── Special ───────────────────────────────────────────────────────────────
    {
        size_t sp_rev, ep_rev;
        if (rcsa_search(rcsa_->rev, rev_p, sp_rev, ep_rev)) {
            const auto sp = grid_.query_special(sp_rev, ep_rev, m);
            if (sp.count > 0 && sp.occ_min_pos < pos_min)
                pos_min = sp.occ_min_pos;
        }
    }

    // ── Crossings ─────────────────────────────────────────────────────────────
    for (size_t i = 1; i < m; ++i) {
        size_t sp_r, ep_r;
        if (!rcsa_search(rcsa_->fwd, pattern.substr(i), sp_r, ep_r)) continue;

        const std::string left_rev(pattern.rbegin() + (m - i), pattern.rend());
        size_t sp_l, ep_l;
        if (!rcsa_search(rcsa_->rev, left_rev, sp_l, ep_l)) continue;

        const auto r = grid_.query_min_2d(sp_r, ep_r, sp_l, ep_l);
        if (r.count == 0) continue;

        if (r.boundary_min < i) continue;
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
    if (!in_alphabet(pattern, alphabet_)) return SIZE_MAX;

    size_t pos_max = 0;
    bool found = false;
    const std::string rev_p(pattern.rbegin(), pattern.rend());

    // ── Special ───────────────────────────────────────────────────────────────
    {
        size_t sp_rev, ep_rev;
        if (rcsa_search(rcsa_->rev, rev_p, sp_rev, ep_rev)) {
            const auto sp = grid_.query_special(sp_rev, ep_rev, m);
            if (sp.count > 0 && sp.occ_max_pos > pos_max) {
                pos_max = sp.occ_max_pos;
                found = true;
            }
        }
    }

    // ── Crossings ─────────────────────────────────────────────────────────────
    for (size_t i = 1; i < m; ++i) {
        size_t sp_r, ep_r;
        if (!rcsa_search(rcsa_->fwd, pattern.substr(i), sp_r, ep_r)) continue;

        const std::string left_rev(pattern.rbegin() + (m - i), pattern.rend());
        size_t sp_l, ep_l;
        if (!rcsa_search(rcsa_->rev, left_rev, sp_l, ep_l)) continue;

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
