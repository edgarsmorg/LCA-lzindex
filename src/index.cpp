#include "index.hpp"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <climits>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <string>

#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>

#include "patricia.h"
#include "dfuds.h"
#include "directcodes.h"

namespace lz77tax {

// ─────────────────────────────────────────────────────────────────────────────
// TrieImpl — PIMPL que guarda los tries DFUDS y las DAC de skips
// ─────────────────────────────────────────────────────────────────────────────

struct LZ77Index::TrieImpl {
    lz77index::basics::dfuds* sst     = nullptr;
    lz77index::basics::dfuds* rev     = nullptr;
    lz77index::basics::FTRep* sst_sk  = nullptr;  // skip lengths SST
    lz77index::basics::FTRep* rev_sk  = nullptr;  // skip lengths rev_trie
    // bp_hb::data y static_doublebitmap_s::bitmap_data apuntan aqui —
    // NO se pueden liberar hasta destruir los dfuds.
    unsigned int*              sst_bm  = nullptr;
    unsigned int*              rev_bm  = nullptr;

    TrieImpl() = default;
    ~TrieImpl() {
        delete sst;  delete rev;
        if (sst_sk) lz77index::basics::destroyFT(sst_sk);
        if (rev_sk) lz77index::basics::destroyFT(rev_sk);
        delete[] sst_bm;  delete[] rev_bm;
    }
    TrieImpl(const TrieImpl&) = delete;
    TrieImpl& operator=(const TrieImpl&) = delete;
};

LZ77Index::LZ77Index()  = default;
LZ77Index::~LZ77Index() = default;

// ─────────────────────────────────────────────────────────────────────────────
// Helpers (anonymous namespace)
// ─────────────────────────────────────────────────────────────────────────────

namespace {

// ── Serialization helpers ─────────────────────────────────────────────────────

void write_bitset256(std::ostream& out, const std::bitset<256>& bs) {
    for (int w = 0; w < 4; ++w) {
        uint64_t word = 0;
        for (int b = 0; b < 64; ++b)
            if (bs[w * 64 + b]) word |= (1ULL << b);
        out.write(reinterpret_cast<const char*>(&word), sizeof(word));
    }
}

void read_bitset256(std::istream& in, std::bitset<256>& bs) {
    bs.reset();
    for (int w = 0; w < 4; ++w) {
        uint64_t word = 0;
        in.read(reinterpret_cast<char*>(&word), sizeof(word));
        for (int b = 0; b < 64; ++b)
            if (word & (1ULL << b)) bs.set(w * 64 + b);
    }
}

// ── Query helpers ─────────────────────────────────────────────────────────────

bool in_alphabet(const std::string& s, const std::bitset<256>& alpha) {
    for (unsigned char c : s)
        if (!alpha.test(c)) return false;
    return true;
}

// ── Trie search ───────────────────────────────────────────────────────────────

// Busca `pattern[0..plen-1]` en el SST.
// Retorna [L,R] 1-indexed (L>R → no encontrado).
// dfuds no tiene metodos const → toma referencias no-const.
void search_sst(lz77index::basics::dfuds& sst,
                lz77index::basics::FTRep* skips,
                const unsigned char* pattern, unsigned int plen,
                unsigned int& L, unsigned int& R) {
    unsigned int node = 1;  // root
    unsigned int pos  = 0;
    unsigned int s_pos;
    while (pos < plen && !sst.inspect(node)) {
        node = sst.labeled_child(node, pattern[pos], &s_pos);
        if (node == static_cast<unsigned int>(-1)) { L = 1; R = 0; return; }
        pos += lz77index::basics::accessFT(skips, s_pos + 1);
    }
    L = sst.leftmost_leaf_rank(node);
    R = sst.rightmost_leaf_rank(node);
}

// Busca `pattern[0..plen-1]` en el rev_trie.
void search_rev(lz77index::basics::dfuds& rev,
                lz77index::basics::FTRep* skips,
                const unsigned char* pattern, unsigned int plen,
                unsigned int& L, unsigned int& R) {
    unsigned int node = 1;
    unsigned int pos  = 0;
    unsigned int s_pos;
    while (pos < plen && !rev.inspect(node)) {
        node = rev.labeled_child(node, pattern[pos], &s_pos);
        if (node == static_cast<unsigned int>(-1)) { L = 1; R = 0; return; }
        pos += lz77index::basics::accessFT(skips, s_pos + 1);
    }
    L = rev.leftmost_leaf_rank(node);
    R = rev.rightmost_leaf_rank(node);
}

// ── Patricia verification ──────────────────────────────────────────────────────

// Verifica que T[occ_pos .. occ_pos+m-1] == pattern.
// boundary es el texto con centinela almacenado en text_s_.
bool verify_occurrence(const std::string& text_s,
                       size_t occ_pos, const std::string& pattern) {
    const size_t m = pattern.size();
    if (occ_pos + m > text_s.size()) return false;
    return text_s.compare(occ_pos, m, pattern) == 0;
}

}  // namespace

// ─────────────────────────────────────────────────────────────────────────────
// build() — construye índice desde texto en memoria
// ─────────────────────────────────────────────────────────────────────────────

void LZ77Index::build(const std::string& text) {
    build_core(text);
    // Sub-índice sobre el texto reverso: su ocurrencia primaria más a la izquierda
    // de P^R es la ocurrencia más a la derecha de P en el texto directo.
    rev_index_ = std::make_unique<LZ77Index>();
    rev_index_->build_core(std::string(text.rbegin(), text.rend()));
}

void LZ77Index::build_core(const std::string& text) {
    text_s_ = text + '\0';
    n_ = text_s_.size();
    LZ77Parsing phrases = LZ77Parser().parse(text_s_);
    z_ = phrases.size();

    alphabet_.reset();
    for (unsigned char c : text) alphabet_.set(c);

    const size_t z1 = z_ - 1;
    trie_ = std::make_unique<TrieImpl>();

    if (z1 == 0) {
        // Texto de 1 frase: grilla vacía, tries no necesarios.
        return;
    }

    assert(z1 < static_cast<size_t>(std::numeric_limits<unsigned int>::max()));

    auto* text_u = const_cast<unsigned char*>(
        reinterpret_cast<const unsigned char*>(text_s_.data()));
    const unsigned int tlen = static_cast<unsigned int>(n_);

    // ── 1. Construir SST (Patricia sobre frases 1..z-1) ──────────────────────
    lz77index::patricia pat_sst(text_u, tlen);
    for (size_t k = 0; k + 1 < z_; ++k)
        pat_sst.add(static_cast<unsigned int>(phrases[k + 1].start_pos));

    unsigned int* ids       = pat_sst.leafIds();
    unsigned int  sst_nodes = pat_sst.nodes();
    unsigned int* sst_bm    = pat_sst.dfuds();
    unsigned char* sst_lab  = pat_sst.dfuds_labels();
    unsigned int* sst_sk    = pat_sst.dfuds_skips();
    trie_->sst    = new lz77index::basics::dfuds(sst_bm, sst_nodes, sst_lab);
    trie_->sst_sk = lz77index::basics::createFT(sst_sk, sst_nodes);
    trie_->sst_bm = sst_bm;  // bp_hb::data apunta aqui — no liberar hasta ~TrieImpl
    delete[] sst_sk;          // createFT copia internamente; sst_sk ya no necesario

    // ── 2. Construir rev_trie (Patricia sobre frases 0..z-2, texto invertido) ─
    lz77index::patricia pat_rev(text_u, tlen);
    // Frase 0 es siempre literal: end_0 = 0, width = 1
    pat_rev.addReverse(0u, 1u);
    for (size_t k = 1; k + 1 < z_; ++k) {
        const unsigned int end_k   = static_cast<unsigned int>(
            phrases[k].start_pos + phrases[k].length);
        const unsigned int width_k = static_cast<unsigned int>(phrases[k].length + 1);
        pat_rev.addReverse(end_k, width_k);
    }

    unsigned int* ids_rev   = pat_rev.leafIds();
    unsigned int  rev_nodes = pat_rev.nodes();
    unsigned int* rev_bm    = pat_rev.dfuds();
    unsigned char* rev_lab  = pat_rev.dfuds_labels();
    unsigned int* rev_sk    = pat_rev.dfuds_skips();
    trie_->rev    = new lz77index::basics::dfuds(rev_bm, rev_nodes, rev_lab);
    trie_->rev_sk = lz77index::basics::createFT(rev_sk, rev_nodes);
    trie_->rev_bm = rev_bm;  // bp_hb::data apunta aqui — no liberar hasta ~TrieImpl
    delete[] rev_sk;

    // ── 3. Computar coordenadas de la grilla desde leaf ranks ─────────────────
    // ids[j]     = insertion_id (= k) de la hoja en SST posicion j (0-indexed)
    // ids_rev[j] = insertion_id (= k) de la hoja en rev_trie posicion j
    //
    // inv_ids_rev[k] = rev_trie rank de frase k
    //   = j tal que ids_rev[j] == k
    //   = se obtiene ordenando idx por ids_rev

    // Después del sort: idx[j] = índice preorder en rev_trie de la inserción j.
    // R[i] = posición preorder (0-indexed) de la inserción ids[i] en el rev_trie
    //      = idx[ids[i]].
    std::vector<unsigned int> idx(z1);
    std::iota(idx.begin(), idx.end(), 0u);
    std::sort(idx.begin(), idx.end(), [&](unsigned int a, unsigned int b) {
        return ids_rev[a] < ids_rev[b];
    });

    // R[i] = rev_trie preorder leaf rank de la frase en SST posicion i
    // boundaries[i] = start_{k+1} para la frase ids[i]
    // phrase_lens[i] = total width de la frase ids[i] (= phrases[k].length + 1)
    std::vector<size_t> R(z1), boundaries(z1), phrase_lens(z1);
    for (size_t i = 0; i < z1; ++i) {
        const size_t k = static_cast<size_t>(ids[i]);
        R[i]           = static_cast<size_t>(idx[k]);
        boundaries[i]  = phrases[k + 1].start_pos;
        phrase_lens[i] = phrases[k].length + 1;
    }

    delete[] ids;
    delete[] ids_rev;

    // ── 4. Construir grilla con coordenadas trie-based ────────────────────────
    grid_.build_from_trie_coords(R, boundaries, phrase_lens);

    phrases = LZ77Parsing{};
}

// ─────────────────────────────────────────────────────────────────────────────
// Accesores de tamaño
// ─────────────────────────────────────────────────────────────────────────────

size_t LZ77Index::trie_bytes() const {
    if (!trie_) return 0;
    size_t total = 0;
    if (trie_->sst)   total += trie_->sst->size();
    if (trie_->rev)   total += trie_->rev->size();
    if (trie_->sst_sk) total += lz77index::basics::sizeFT(trie_->sst_sk);
    if (trie_->rev_sk) total += lz77index::basics::sizeFT(trie_->rev_sk);
    return total;
}

size_t LZ77Index::index_bytes() const {
    const auto bd = grid_.size_breakdown();
    size_t total = bd.wm + bd.bv_fwd + bd.bv_rev + bd.rank_fwd + bd.rank_rev
                 + bd.text_pos + bd.phrase_total_len + bd.wm_min_rmq
                 + trie_bytes();
    if (rev_index_) total += rev_index_->index_bytes();
    return total;
}

// ─────────────────────────────────────────────────────────────────────────────
// save() / load()
// ─────────────────────────────────────────────────────────────────────────────

void LZ77Index::save(const std::filesystem::path& prefix) const {
    // 1. Meta: n_, z_, alphabet_
    {
        std::ofstream f(std::filesystem::path(prefix).replace_extension(".meta"),
                        std::ios::binary);
        sdsl::write_member(n_, f);
        sdsl::write_member(z_, f);
        write_bitset256(f, alphabet_);
        const uint8_t has_rev = rev_index_ ? 1 : 0;
        sdsl::write_member(has_rev, f);
    }
    // 2. Grid
    {
        std::ofstream f(std::filesystem::path(prefix).replace_extension(".grid"),
                        std::ios::binary);
        grid_.serialize(f);
    }
    // 3. Tries: usar FILE* API nativa de uiHRDC
    if (trie_) {
        auto trie_path = std::filesystem::path(prefix).replace_extension(".trie");
        FILE* fp = fopen(trie_path.c_str(), "wb");
        if (!fp) throw std::runtime_error("No se pudo abrir " + trie_path.string());
        if (trie_->sst)    trie_->sst->save(fp);
        if (trie_->sst_sk) lz77index::basics::saveFT(trie_->sst_sk, fp);
        if (trie_->rev)    trie_->rev->save(fp);
        if (trie_->rev_sk) lz77index::basics::saveFT(trie_->rev_sk, fp);
        fclose(fp);
    }
    // 4. Sub-índice reverso (con su propio prefijo "<prefix>_rev.*")
    if (rev_index_)
        rev_index_->save(std::filesystem::path(prefix.string() + "_rev"));
}

void LZ77Index::load(const std::filesystem::path& prefix, const std::string& text) {
    text_s_ = text + '\0';
    uint8_t has_rev = 0;
    // 1. Meta
    {
        std::ifstream f(std::filesystem::path(prefix).replace_extension(".meta"),
                        std::ios::binary);
        sdsl::read_member(n_, f);
        sdsl::read_member(z_, f);
        read_bitset256(f, alphabet_);
        sdsl::read_member(has_rev, f);
    }
    // 2. Grid
    {
        std::ifstream f(std::filesystem::path(prefix).replace_extension(".grid"),
                        std::ios::binary);
        grid_.load(f);
    }
    // 3. Tries
    trie_ = std::make_unique<TrieImpl>();
    if (z_ > 1) {
        auto trie_path = std::filesystem::path(prefix).replace_extension(".trie");
        FILE* fp = fopen(trie_path.c_str(), "rb");
        if (!fp) throw std::runtime_error("No se pudo abrir " + trie_path.string());
        trie_->sst    = lz77index::basics::dfuds::load(fp);
        trie_->sst_sk = lz77index::basics::loadFT(fp);
        trie_->rev    = lz77index::basics::dfuds::load(fp);
        trie_->rev_sk = lz77index::basics::loadFT(fp);
        fclose(fp);
    }
    // 4. Sub-índice reverso
    if (has_rev) {
        rev_index_ = std::make_unique<LZ77Index>();
        rev_index_->load(std::filesystem::path(prefix.string() + "_rev"),
                         std::string(text.rbegin(), text.rend()));
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// count()
// ─────────────────────────────────────────────────────────────────────────────

size_t LZ77Index::count(const std::string& pattern) const {
    const size_t m = pattern.size();
    if (m < 2 || grid_.point_count() == 0 || !trie_) return 0;
    if (!in_alphabet(pattern, alphabet_)) return 0;

    const std::string rev_p(pattern.rbegin(), pattern.rend());
    const auto* rev_u  = reinterpret_cast<const unsigned char*>(rev_p.data());
    const auto* pat_u  = reinterpret_cast<const unsigned char*>(pattern.data());

    size_t total = 0;

    // ── Special: patrón completo en rev_trie, buscamos frases que lo contengan ─
    if (trie_->rev) {
        unsigned int Ll, Lr;
        search_rev(*trie_->rev, trie_->rev_sk, rev_u,
                   static_cast<unsigned int>(m), Ll, Lr);
        if (Lr >= Ll) {
            const auto sp = grid_.query_special_direct(Ll - 1, Lr - 1, m);
            if (sp.count > 0) {
                // Verificar una ocurrencia representante
                if (sp.occ_min_pos != SIZE_MAX &&
                    verify_occurrence(text_s_, sp.occ_min_pos, pattern))
                    total += sp.count;
            }
        }
    }

    // ── Crossings: split P[0..i-1] | P[i..m-1] para i=1..m-1 ────────────────
    for (size_t i = 1; i < m; ++i) {
        // Mitad derecha P[i..m-1] en SST
        unsigned int Rl, Rr;
        search_sst(*trie_->sst, trie_->sst_sk,
                   pat_u + i, static_cast<unsigned int>(m - i), Rl, Rr);
        if (Rr < Rl) continue;

        // Mitad izquierda invertida en rev_trie
        const std::string left_rev(pattern.rbegin() + static_cast<std::ptrdiff_t>(m - i),
                                   pattern.rend());
        const auto* lrev_u = reinterpret_cast<const unsigned char*>(left_rev.data());
        unsigned int Ll, Lr;
        search_rev(*trie_->rev, trie_->rev_sk,
                   lrev_u, static_cast<unsigned int>(i), Ll, Lr);
        if (Lr < Ll) continue;

        const auto res = grid_.query_direct(Rl - 1, Rr - 1, Ll - 1, Lr - 1);
        if (res.first == 0) continue;

        // Contar puntos validos: phrase_total_len >= i (ocurrencia empieza en frase k).
        // Patricia: una verificacion basta para confirmar todo el rango.
        bool verified = false;
        size_t valid_count = 0;
        for (const auto& [wt_idx, y_rel] : res.second) {
            if (grid_.phrase_total_len(wt_idx) < i) continue;
            ++valid_count;
            if (!verified) {
                const size_t boundary = grid_.text_pos(wt_idx);
                if (boundary >= i && verify_occurrence(text_s_, boundary - i, pattern))
                    verified = true;
            }
        }
        if (verified) total += valid_count;
    }

    return total;
}

// ─────────────────────────────────────────────────────────────────────────────
// locate_extremal()
// ─────────────────────────────────────────────────────────────────────────────

// Ocurrencia primaria más a la izquierda de `pattern` (crossings + end-aligned),
// o SIZE_MAX si no ocurre. Solo usa el RMQ de mínimo (nunca de máximo).
size_t LZ77Index::locate_leftmost(const std::string& pattern) const {
    const size_t m = pattern.size();
    if (m < 2 || grid_.point_count() == 0 || !trie_) return SIZE_MAX;
    if (!in_alphabet(pattern, alphabet_)) return SIZE_MAX;

    const std::string rev_p(pattern.rbegin(), pattern.rend());
    const auto* rev_u = reinterpret_cast<const unsigned char*>(rev_p.data());
    const auto* pat_u = reinterpret_cast<const unsigned char*>(pattern.data());

    size_t pos_min = SIZE_MAX;

    // ── Special (end-aligned): P completo termina al final de una frase ─────────
    if (trie_->rev) {
        unsigned int Ll, Lr;
        search_rev(*trie_->rev, trie_->rev_sk, rev_u,
                   static_cast<unsigned int>(m), Ll, Lr);
        if (Lr >= Ll) {
            const auto sp = grid_.query_special_direct(Ll - 1, Lr - 1, m);
            if (sp.count > 0 && sp.occ_min_pos != SIZE_MAX &&
                verify_occurrence(text_s_, sp.occ_min_pos, pattern) &&
                sp.occ_min_pos < pos_min) {
                pos_min = sp.occ_min_pos;
            }
        }
    }

    // ── Crossings: split P[0..i-1] | P[i..m-1] para i=1..m-1 ────────────────────
    for (size_t i = 1; i < m; ++i) {
        unsigned int Rl, Rr;
        search_sst(*trie_->sst, trie_->sst_sk,
                   pat_u + i, static_cast<unsigned int>(m - i), Rl, Rr);
        if (Rr < Rl) continue;

        const std::string left_rev(pattern.rbegin() + static_cast<std::ptrdiff_t>(m - i),
                                   pattern.rend());
        const auto* lrev_u = reinterpret_cast<const unsigned char*>(left_rev.data());
        unsigned int Ll, Lr;
        search_rev(*trie_->rev, trie_->rev_sk,
                   lrev_u, static_cast<unsigned int>(i), Ll, Lr);
        if (Lr < Ll) continue;

        const auto rmin = grid_.query_min_direct(Rl - 1, Rr - 1, Ll - 1, Lr - 1);
        if (rmin.count == 0 || rmin.wt_idx == SIZE_MAX) continue;

        // El punto de mínimo boundary del RMQ. Dos casos:
        size_t occ = SIZE_MAX;
        if (grid_.phrase_total_len(rmin.wt_idx) >= i && rmin.boundary_min >= i) {
            // (a) El extremo del RMQ es un crossing válido. UNA verificación
            //     Patricia decide todo el rectángulo: si falla, ninguna ocurrencia
            //     calza con el patrón → no hay que enumerar.
            const size_t cand = rmin.boundary_min - i;
            if (verify_occurrence(text_s_, cand, pattern)) occ = cand;
        } else {
            // (b) El extremo del RMQ no cumple el filtro de largo de frase; puede
            //     haber otro punto (boundary mayor) que sí lo cumpla. Se busca el de
            //     mínimo boundary entre los válidos por largo y se verifica una vez.
            const auto mf = grid_.query_min_filtered(Rl - 1, Rr - 1, Ll - 1, Lr - 1, i);
            if (mf.count > 0 && mf.boundary_min >= i) {
                const size_t cand = mf.boundary_min - i;
                if (verify_occurrence(text_s_, cand, pattern)) occ = cand;
            }
        }

        if (occ != SIZE_MAX && occ < pos_min) pos_min = occ;
    }

    return pos_min;
}

std::pair<size_t, size_t> LZ77Index::locate_extremal(const std::string& pattern) const {
    const size_t p_min = locate_leftmost(pattern);
    if (p_min == SIZE_MAX) return {SIZE_MAX, 0};

    // Más a la derecha en T = más a la izquierda de P^R en T^R, remapeada.
    // Si una ocurrencia de P empieza en s en T (largo m), en T^R (largo L=n_-1)
    // la ocurrencia de P^R empieza en j = L - m - s; luego s = L - m - j.
    size_t p_max = p_min;
    if (rev_index_) {
        const std::string rp(pattern.rbegin(), pattern.rend());
        const size_t j = rev_index_->locate_leftmost(rp);
        if (j != SIZE_MAX) {
            const size_t L = n_ - 1;               // largo del texto sin centinela
            p_max = L - pattern.size() - j;
        }
    }
    return {p_min, p_max};
}

}  // namespace lz77tax
