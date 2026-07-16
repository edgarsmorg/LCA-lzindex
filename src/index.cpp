#include "index.hpp"
#include "lz77/text_oracle.hpp"

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
#include <sdsl/util.hpp>
#include <sdsl/divsufsort.hpp>

namespace lz77tax {

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

void write_byte_vec(std::ostream& out, const std::vector<uint8_t>& v) {
    uint64_t sz = v.size();
    out.write(reinterpret_cast<const char*>(&sz), sizeof(sz));
    if (sz) out.write(reinterpret_cast<const char*>(v.data()), sz);
}

void read_byte_vec(std::istream& in, std::vector<uint8_t>& v) {
    uint64_t sz = 0;
    in.read(reinterpret_cast<char*>(&sz), sizeof(sz));
    v.resize(sz);
    if (sz) in.read(reinterpret_cast<char*>(v.data()), sz);
}

// ── Query helpers ─────────────────────────────────────────────────────────────

bool in_alphabet(const std::string& s, const std::bitset<256>& alpha) {
    for (unsigned char c : s)
        if (!alpha.test(c)) return false;
    return true;
}

// ── Búsqueda binaria sobre sparse suffix arrays (PSA) ───────────────────────────
//
// Reemplaza los tries Patricia: cada PSA es el arreglo de las z-1 posiciones de
// frase ordenadas por su sufijo. La comparación se hace contra el texto real
// (text_s, con centinela '\0' al final) → búsqueda EXACTA, sin falsos positivos
// que verificar. El índice en la PSA ES la coordenada de la grilla (X para
// psa_fwd, Y para psa_rev).

// Los comparadores devuelven el signo Y el LCP alcanzado (cuantos caracteres
// coincidieron), porque la busqueda binaria los necesita para acelerarse.
// <0 si el sufijo es menor que pat, 0 si pat es prefijo del sufijo, >0 si mayor.
// `skip` = caracteres que YA se saben iguales: la comparacion arranca ahi.

// Sufijo forward text_s[p..] contra pat[0..plen-1].
template <class TextT>
inline int cmp_fwd(const TextT& text_s, size_t p,
                   const unsigned char* pat, size_t plen,
                   size_t skip, size_t& lcp) {
    const size_t n = text_s.size();
    for (size_t k = skip; k < plen; ++k) {
        const unsigned char a = (p + k < n)
            ? static_cast<unsigned char>(text_s[p + k]) : 0;  // 0 = centinela
        const unsigned char b = pat[k];
        if (a != b) { lcp = k; return a < b ? -1 : 1; }
    }
    lcp = plen;
    return 0;
}

// Cadena hacia atras en text_s terminando en e (text_s[e], [e-1], …) contra
// reverse(pat[0..plen-1]). Reproduce el orden de sufijos de rev_s = reverse(T)+'\0'
// usado al construir psa_rev. Equivale a "pat termina en la posicion e".
template <class TextT>
inline int cmp_rev(const TextT& text_s, size_t e,
                   const unsigned char* pat, size_t plen,
                   size_t skip, size_t& lcp) {
    for (size_t k = skip; k < plen; ++k) {
        const unsigned char a = (e >= k)
            ? static_cast<unsigned char>(text_s[e - k]) : 0;  // 0 = centinela
        const unsigned char b = pat[plen - 1 - k];
        if (a != b) { lcp = k; return a < b ? -1 : 1; }
    }
    lcp = plen;
    return 0;
}

// Sobrecargas para los accesores: extraen por bloques, amortizando el descenso
// recursivo sobre varios caracteres. Sirven igual al accesor directo y a su
// vista espejada (el sub-indice reverso).
//
// Extraccion PEREZOSA: en la busqueda binaria la comparacion falla muy pronto
// (medido: ~3 caracteres mas alla de `skip`), asi que extraer los m caracteres de
// una vez desperdicia el grueso del trabajo. Se pide primero un bloque corto y
// solo si empatan se pide el resto.
constexpr size_t kCmpChunk = 8;

// Muestra por LIMITE DE FRASE, indexada por coordenada de grilla. Como toda
// comparacion de la busqueda binaria arranca en un limite (X hacia adelante,
// Y hacia atras), guardar los primeros kSmp caracteres alli resuelve el 95 % de
// las comparaciones sin tocar la recursion del accesor (medido). El indice de
// grilla es el rango de la propia busqueda, asi que el acceso es O(1) (no hay
// phrase_of). El resto cae al extract() normal. Costo: 2 kSmp bytes por frase.
// kSmp = caracteres de muestra por limite de frase. 0 (default) = config
// COMPACTO (indice mas chico que el sr). 8 = config VELOCIDAD (locate mas rapido
// que el sr, a costa de +espacio por la muestra). Es un trade espacio/tiempo.
#ifndef LZ_SAMPLE
#define LZ_SAMPLE 0
#endif
constexpr size_t kSmp = LZ_SAMPLE;

// `smp` (o nullptr) = kSmp caracteres de T desde el limite, en orden forward.
// `slen` = cuantos son validos (< kSmp cerca de los bordes del texto).
template <class OracleT>
inline int cmp_fwd_blk(const OracleT& o, size_t p,
                       const unsigned char* pat, size_t plen,
                       size_t skip, size_t& lcp,
                       const uint8_t* smp, size_t slen) {
    uint8_t buf[512];
    const size_t n = o.size();
    size_t done = skip;

    if (smp) {                                  // camino rapido O(1)
        const size_t lim = std::min(plen, slen);
        for (size_t k = done; k < lim; ++k)
            if (smp[k] != pat[k]) { lcp = k; return smp[k] < pat[k] ? -1 : 1; }
        if (plen <= slen) { lcp = plen; return 0; }
        done = slen;
    }

    while (done < plen) {
        const size_t want  = std::min(done == skip ? kCmpChunk : plen - done,
                                      plen - done);
        const size_t avail = (p + done < n) ? std::min(want, n - p - done) : 0;
        if (avail) o.extract(p + done, avail, buf);
        for (size_t k = 0; k < want; ++k) {
            const unsigned char a = (k < avail) ? buf[k] : 0;
            const unsigned char b = pat[done + k];
            if (a != b) { lcp = done + k; return a < b ? -1 : 1; }
        }
        done += want;
    }
    lcp = plen;
    return 0;
}
// `smp` = kSmp caracteres de T terminando en e, en orden forward (smp[slen-1]=T[e]).
template <class OracleT>
inline int cmp_rev_blk(const OracleT& o, size_t e,
                       const unsigned char* pat, size_t plen,
                       size_t skip, size_t& lcp,
                       const uint8_t* smp, size_t slen) {
    uint8_t buf[512];
    size_t done = skip;

    if (smp) {                                  // camino rapido O(1)
        const size_t lim = std::min(plen, slen);
        for (size_t k = done; k < lim; ++k) {
            const unsigned char a = smp[slen - 1 - k];      // T[e-k]
            const unsigned char b = pat[plen - 1 - k];
            if (a != b) { lcp = k; return a < b ? -1 : 1; }
        }
        if (plen <= slen) { lcp = plen; return 0; }
        done = slen;
    }

    while (done < plen) {
        const size_t want  = std::min(done == skip ? kCmpChunk : plen - done,
                                      plen - done);
        const size_t hi    = (e >= done) ? e - done : 0;
        const size_t avail = (e >= done) ? std::min(want, hi + 1) : 0;
        if (avail) o.extract(hi + 1 - avail, avail, buf);
        for (size_t k = 0; k < want; ++k) {
            const unsigned char a = (k < avail) ? buf[avail - 1 - k] : 0;
            const unsigned char b = pat[plen - 1 - done - k];
            if (a != b) { lcp = done + k; return a < b ? -1 : 1; }
        }
        done += want;
    }
    lcp = plen;
    return 0;
}

// Dispatch: el texto plano (string) lee directo e ignora la muestra; los
// accesores usan la muestra + extract. La sobrecarga se elige por el tipo del
// texto en tiempo de compilacion, asi cada rama compila solo lo que su tipo
// soporta (el string no tiene extract()).
template <class TextT>
inline int cmp_fwd_s(const TextT& t, size_t p, const unsigned char* pat, size_t plen,
                     size_t skip, size_t& lcp, const uint8_t*, size_t) {
    return cmp_fwd(t, p, pat, plen, skip, lcp);
}
template <class TextT>
inline int cmp_rev_s(const TextT& t, size_t e, const unsigned char* pat, size_t plen,
                     size_t skip, size_t& lcp, const uint8_t*, size_t) {
    return cmp_rev(t, e, pat, plen, skip, lcp);
}
inline int cmp_fwd_s(const TextOracle& o, size_t p, const unsigned char* pat, size_t plen,
                     size_t skip, size_t& lcp, const uint8_t* smp, size_t slen) {
    return cmp_fwd_blk(o, p, pat, plen, skip, lcp, smp, slen);
}
inline int cmp_rev_s(const TextOracle& o, size_t e, const unsigned char* pat, size_t plen,
                     size_t skip, size_t& lcp, const uint8_t* smp, size_t slen) {
    return cmp_rev_blk(o, e, pat, plen, skip, lcp, smp, slen);
}
inline int cmp_fwd_s(const MirrorOracle& o, size_t p, const unsigned char* pat, size_t plen,
                     size_t skip, size_t& lcp, const uint8_t* smp, size_t slen) {
    return cmp_fwd_blk(o, p, pat, plen, skip, lcp, smp, slen);
}
inline int cmp_rev_s(const MirrorOracle& o, size_t e, const unsigned char* pat, size_t plen,
                     size_t skip, size_t& lcp, const uint8_t* smp, size_t slen) {
    return cmp_rev_blk(o, e, pat, plen, skip, lcp, smp, slen);
}

// Rango [L,R] 0-indexed inclusivo de ranks r en [0,z1) con cmp(pos(r)) == 0.
// R < L (L=1, R=0) si no hay ninguna.
//
// Busqueda binaria ACELERADA (Manber-Myers, variante sin arreglo LCP): se
// mantienen l = lcp(pat, sufijo[lo]) y r = lcp(pat, sufijo[hi]) con lo/hi los
// bordes actuales. Para todo mid entre ellos vale
//     lcp(pat, sufijo[mid]) >= min(l, r)
// (los sufijos estan ordenados y pat queda encajonado entre ambos bordes), asi
// que esos min(l,r) caracteres no se vuelven a comparar. Con el accesor eso se
// traduce en pedirle menos texto. lo/hi arrancan en centinelas virtuales
// (-1 y z1), donde el lcp es 0 por definicion.
// cmp(rank, skip, lcp) compara pat contra el sufijo del rango `rank`, saltando
// los primeros `skip` caracteres (ya sabidos iguales) y devolviendo el lcp.
template <class Cmp>
void bsearch_ranks(size_t z1, size_t plen, Cmp cmp, size_t& L, size_t& R) {
    if (z1 == 0) { L = 1; R = 0; return; }

    // ── Cota inferior: primer rank con cmp >= 0 ───────────────────────────────
    int64_t lo = -1, hi = static_cast<int64_t>(z1);
    size_t  l = 0,  r = 0;
    while (hi - lo > 1) {
        const int64_t mid = lo + (hi - lo) / 2;
        size_t k = 0;
        const int c = cmp(static_cast<size_t>(mid), std::min(l, r), k);
        if (c < 0) { lo = mid; l = k; } else { hi = mid; r = k; }
    }
    const int64_t first = hi;
    // Al salir, r = lcp(pat, sufijo[first]): hay rango solo si coincide entero.
    if (first >= static_cast<int64_t>(z1) || r != plen) { L = 1; R = 0; return; }

    // ── Cota superior: ultimo rank con cmp == 0 ───────────────────────────────
    // sufijo[first] coincide con pat => sirve de borde izquierdo con l = plen.
    lo = first; l = plen;
    hi = static_cast<int64_t>(z1); r = 0;
    while (hi - lo > 1) {
        const int64_t mid = lo + (hi - lo) / 2;
        size_t k = 0;
        const int c = cmp(static_cast<size_t>(mid), std::min(l, r), k);
        if (c <= 0) { lo = mid; l = k; } else { hi = mid; r = k; }
    }
    L = static_cast<size_t>(first);
    R = static_cast<size_t>(lo);
}

// [L,R] = X-ranks i con pat[0..plen-1] como prefijo de text_s[start_{k+1}..].
// La posición del X-rank i ES grid.text_pos(i) (= psa_fwd, sin copia redundante).
// `fwd_smp` (o nullptr) = muestra por punto de grilla en orden X.
template <class TextT>
void search_sst(const TextT& text_s, const Grid2D& grid, size_t z1,
                const unsigned char* pat, size_t plen, size_t& L, size_t& R,
                const uint8_t* fwd_smp, size_t n) {
    bsearch_ranks(z1, plen,
        [&](size_t i, size_t skip, size_t& lcp) {
            const size_t p = grid.text_pos(i);
            const uint8_t* s = fwd_smp ? fwd_smp + i * kSmp : nullptr;
            const size_t sl = fwd_smp ? std::min(kSmp, n - p) : 0;
            return cmp_fwd_s(text_s, p, pat, plen, skip, lcp, s, sl); }, L, R);
}

// [L,R] = Y-ranks j tales que pat[0..plen-1] termina en la posición psa_rev[j].
// `bwd_smp` (o nullptr) = muestra por punto de grilla en orden Y.
template <class TextT>
void search_rev(const TextT& text_s, const sdsl::int_vector<>& psa_rev,
                const unsigned char* pat, size_t plen, size_t& L, size_t& R,
                const uint8_t* bwd_smp) {
    bsearch_ranks(psa_rev.size(), plen,
        [&](size_t j, size_t skip, size_t& lcp) {
            const size_t e = static_cast<size_t>(psa_rev[j]);
            const uint8_t* s = bwd_smp ? bwd_smp + j * kSmp : nullptr;
            const size_t sl = bwd_smp ? std::min(kSmp, e + 1) : 0;
            return cmp_rev_s(text_s, e, pat, plen, skip, lcp, s, sl); }, L, R);
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
    if (z1 == 0) {
        // Texto de 1 frase: grilla y PSA vacías.
        return;
    }

    // ── 1. Ranks de sufijo forward/reverso de las z-1 posiciones de frase ──────
    // X-rank de la frase k = rank de start_{k+1} por sufijo forward de text_s.
    // Y-rank de la frase k = rank de end_k por sufijo hacia atrás
    //   = rank del sufijo forward de rev_s = reverse(T)+'\0' en pos (n-2-end_k).
    // Se construyen SA/ISA con divsufsort (scratch, O(n)) y se descartan.
    // ponytail: el SA forward duplica el que LZ77Parser ya arma internamente;
    // exponer su ISA ahorraría un divsufsort por índice si el build de 1.8 GB molesta.
    std::vector<int32_t> fkey(z1), rkey(z1);
    {
        std::vector<int32_t> sa(n_), isa(n_);
        // forward
        sdsl::divsufsort(reinterpret_cast<const uint8_t*>(text_s_.data()),
                         sa.data(), static_cast<int32_t>(n_));
        for (size_t i = 0; i < n_; ++i) isa[sa[i]] = static_cast<int32_t>(i);
        for (size_t k = 0; k < z1; ++k)
            fkey[k] = isa[phrases[k + 1].start_pos];
        // reverso: rev_s = reverse(T) + '\0'
        std::string rev_s(text.rbegin(), text.rend());
        rev_s.push_back('\0');
        sdsl::divsufsort(reinterpret_cast<const uint8_t*>(rev_s.data()),
                         sa.data(), static_cast<int32_t>(n_));
        for (size_t i = 0; i < n_; ++i) isa[sa[i]] = static_cast<int32_t>(i);
        for (size_t k = 0; k < z1; ++k) {
            const size_t end_k = phrases[k].start_pos + phrases[k].length;
            rkey[k] = isa[n_ - 2 - end_k];
        }
    }

    // ── 2. Ordenar las z-1 frases por X-rank y por Y-rank ──────────────────────
    std::vector<uint32_t> order_fwd(z1), order_rev(z1);
    std::iota(order_fwd.begin(), order_fwd.end(), 0u);
    std::iota(order_rev.begin(), order_rev.end(), 0u);
    std::sort(order_fwd.begin(), order_fwd.end(),
              [&](uint32_t a, uint32_t b) { return fkey[a] < fkey[b]; });
    std::sort(order_rev.begin(), order_rev.end(),
              [&](uint32_t a, uint32_t b) { return rkey[a] < rkey[b]; });
    std::vector<size_t> yrank(z1);
    for (size_t j = 0; j < z1; ++j) yrank[order_rev[j]] = j;

    // ── 3. Coordenadas de la grilla (indexadas por X-rank i) ──────────────────
    std::vector<size_t> R(z1), boundaries(z1), phrase_lens(z1);
    for (size_t i = 0; i < z1; ++i) {
        const size_t k = order_fwd[i];
        R[i]           = yrank[k];
        boundaries[i]  = phrases[k + 1].start_pos;
        phrase_lens[i] = phrases[k].length + 1;
    }

    // ── 4. Construir grilla con coordenadas ya relativas y derivar psa_rev ─────
    // psa_fwd NO existe: la posición del X-rank i es grid_.text_pos(i).
    // psa_rev se deriva de la grilla (no se serializa): ver rebuild_psa_rev().
    grid_.build(R, boundaries, phrase_lens);
    rebuild_psa_rev();

    phrases = LZ77Parsing{};
}

// psa_rev[j] = posición end_k de la frase con Y-rank j, reconstruida desde la
// grilla: para cada X-rank i, su Y-rank es wm[i] y su end es text_pos(i)-1
// (boundaries consecutivos). O(z log σ); evita serializar psa_rev.
void LZ77Index::rebuild_psa_rev() {
    const size_t z1 = grid_.point_count();
    if (z1 == 0) { psa_rev_ = sdsl::int_vector<>(); return; }
    const uint8_t w = static_cast<uint8_t>(sdsl::bits::hi(n_) + 1);
    psa_rev_ = sdsl::int_vector<>(z1, 0, w);
    const auto& wm = grid_.wm();
    for (size_t i = 0; i < z1; ++i)
        psa_rev_[wm[i]] = grid_.text_pos(i) - 1;
}

// ─────────────────────────────────────────────────────────────────────────────
// Accesores de tamaño
// ─────────────────────────────────────────────────────────────────────────────

// Bytes del sustrato de búsqueda de esta dirección: solo psa_rev en RAM
// (psa_fwd no existe — es grid.text_pos; psa_rev NO se serializa, se reconstruye
// en load). Se mantiene el nombre por compatibilidad con las tools de medición.
size_t LZ77Index::trie_bytes() const {
    return sdsl::size_in_bytes(psa_rev_);
}

size_t LZ77Index::index_bytes() const {
    // Cuenta solo lo que se SERIALIZA (mismo criterio que el tamaño en disco):
    // NI phrase_total_len NI psa_rev (ambos derivables de text_pos, reconstruidos
    // al cargar). El accesor SÍ (start_ excluido, derivable de la grilla).
    const auto bd = grid_.size_breakdown();
    size_t total = bd.wm + bd.text_pos + bd.wm_min_rmq
                 + (oracle_ ? oracle_->size_in_bytes() : 0)
                 + fwd_smp_.size() + bwd_smp_.size();
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
    // 3. Accesor (.orac): el índice se vuelve un self-index — se carga sin el
    //    texto. NO serializa start_ (derivable de la grilla) ni el MirrorOracle
    //    del sub-índice reverso (comparte el accesor directo). Sí las muestras si
    //    las hay (config velocidad). Si no hay accesor construido, no se escribe.
    if (oracle_) {
        std::ofstream f(std::filesystem::path(prefix).replace_extension(".orac"),
                        std::ios::binary);
        oracle_->serialize(f);
    }
    {
        std::ofstream f(std::filesystem::path(prefix).replace_extension(".smp"),
                        std::ios::binary);
        write_byte_vec(f, fwd_smp_);
        write_byte_vec(f, bwd_smp_);
    }
    // 4. PSA: NO se serializan. psa_fwd = grid.text_pos (ya en .grid);
    //    psa_rev se reconstruye desde la grilla en load() (rebuild_psa_rev).
    // 5. Sub-índice reverso (con su propio prefijo "<prefix>_rev.*")
    if (rev_index_)
        rev_index_->save(std::filesystem::path(prefix.string() + "_rev"));
}

// [0] + sorted(grid.text_pos): los inicios de frase, en orden de posición.
std::vector<size_t> LZ77Index::derive_starts() const {
    const size_t z1 = grid_.point_count();
    std::vector<size_t> starts;
    starts.reserve(z1 + 1);
    starts.push_back(0);
    for (size_t i = 0; i < z1; ++i) starts.push_back(grid_.text_pos(i));
    std::sort(starts.begin() + 1, starts.end());
    return starts;
}

void LZ77Index::load(const std::filesystem::path& prefix, const std::string& text) {
    // Self-index si existe el accesor serializado (.orac): se carga sin el texto.
    // Si no (índice viejo), modo compatibilidad con el texto.
    const bool self_index = std::filesystem::exists(
        std::filesystem::path(prefix).replace_extension(".orac"));
    load_rec(prefix, self_index ? std::string() : text, self_index, nullptr);
}

void LZ77Index::load_rec(const std::filesystem::path& prefix, const std::string& text,
                         bool self_index, const TextOracle* parent_oracle) {
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
    // 3. PSA: psa_rev se reconstruye desde la grilla (no se serializa).
    rebuild_psa_rev();
    // 4. Fuente de texto para las consultas
    if (self_index) {
        // Muestras (config velocidad); vacías en config compacto.
        {
            const auto smp = std::filesystem::path(prefix).replace_extension(".smp");
            if (std::filesystem::exists(smp)) {
                std::ifstream f(smp, std::ios::binary);
                read_byte_vec(f, fwd_smp_);
                read_byte_vec(f, bwd_smp_);
            }
        }
        if (parent_oracle == nullptr) {
            // Índice directo: carga su accesor; start_ se deriva de la grilla.
            oracle_ = std::make_unique<TextOracle>();
            std::ifstream f(std::filesystem::path(prefix).replace_extension(".orac"),
                            std::ios::binary);
            oracle_->load(f, derive_starts());
        } else {
            // Sub-índice reverso: espejo del accesor del padre (no carga texto).
            mirror_ = std::make_unique<MirrorOracle>(*parent_oracle, n_);
        }
        use_oracle_ = true;
    } else {
        text_s_ = text + '\0';   // compatibilidad: sin .orac se necesita el texto
    }
    // 5. Sub-índice reverso
    if (has_rev) {
        rev_index_ = std::make_unique<LZ77Index>();
        rev_index_->load_rec(std::filesystem::path(prefix.string() + "_rev"),
                             self_index ? std::string()
                                        : std::string(text.rbegin(), text.rend()),
                             self_index,
                             self_index ? oracle_.get() : nullptr);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Accesor al texto (experimental): evita almacenar T
// ─────────────────────────────────────────────────────────────────────────────

void LZ77Index::build_oracle() {
    if (oracle_ || mirror_) return;   // ya cargado/construido (p.ej. self-index)
    if (text_s_.empty()) return;
    std::vector<uint8_t> syms;
    for (int c = 0; c < 256; ++c)
        if (alphabet_.test(c)) syms.push_back(static_cast<uint8_t>(c));
    syms.push_back(0);                       // centinela '\0' de text_s_
    std::sort(syms.begin(), syms.end());
    syms.erase(std::unique(syms.begin(), syms.end()), syms.end());
    LZ77Parsing ph = LZ77Parser().parse(text_s_);
    oracle_ = std::make_unique<TextOracle>();
    oracle_->build(ph, text_s_.size(), syms);
    build_samples();
    // El sub-índice reverso NO construye accesor propio: T^R se lee desde este
    // mismo accesor con el remapeo i -> L-1-i (ver MirrorOracle).
    if (rev_index_) {
        rev_index_->mirror_ = std::make_unique<MirrorOracle>(*oracle_, n_);
        rev_index_->build_samples();
    }
}

// Muestra por punto de grilla: kSmp caracteres de T a cada lado del limite de
// frase (hacia adelante desde el X-punto, hacia atras desde el Y-punto). Se
// llena desde text_s_ (disponible en RAM al construir). En consulta permite el
// camino rapido O(1) sin recursion del accesor.
void LZ77Index::build_samples() {
    if (kSmp == 0 || grid_.point_count() == 0) return;
    const size_t z1 = grid_.point_count();
    const size_t n  = text_s_.size();

    fwd_smp_.assign(z1 * kSmp, 0);
    for (size_t i = 0; i < z1; ++i) {
        const size_t p  = grid_.text_pos(i);
        const size_t fl = std::min(kSmp, n - p);
        for (size_t t = 0; t < fl; ++t)
            fwd_smp_[i * kSmp + t] = static_cast<uint8_t>(text_s_[p + t]);
    }

    bwd_smp_.assign(z1 * kSmp, 0);
    for (size_t j = 0; j < z1; ++j) {
        const size_t e  = static_cast<size_t>(psa_rev_[j]);
        const size_t sl = std::min(kSmp, e + 1);
        for (size_t t = 0; t < sl; ++t)              // orden forward: smp[sl-1]=T[e]
            bwd_smp_[j * kSmp + t] = static_cast<uint8_t>(text_s_[e - sl + 1 + t]);
    }
}

void LZ77Index::use_oracle(bool on) {
    use_oracle_ = on;
    if (rev_index_) rev_index_->use_oracle(on);
}

size_t LZ77Index::oracle_bytes() const {
    // mirror_ no aporta bytes: es una vista del accesor directo. Las muestras SI
    // (cada sub-indice tiene la suya, sobre su propio texto).
    size_t b = oracle_ ? oracle_->size_in_bytes() : 0;
    b += fwd_smp_.size() + bwd_smp_.size();
    if (rev_index_) b += rev_index_->oracle_bytes();
    return b;
}

// ─────────────────────────────────────────────────────────────────────────────
// count()
// ─────────────────────────────────────────────────────────────────────────────

size_t LZ77Index::count(const std::string& pattern) const {
    const size_t m = pattern.size();
    if (m < 2 || grid_.point_count() == 0) return 0;
    if (!in_alphabet(pattern, alphabet_)) return 0;

    // Misma fuente de texto que locate: accesor (self-index) o texto plano.
    const uint8_t* fsmp = (use_oracle_ && !fwd_smp_.empty()) ? fwd_smp_.data() : nullptr;
    const uint8_t* bsmp = (use_oracle_ && !bwd_smp_.empty()) ? bwd_smp_.data() : nullptr;
    if (use_oracle_) {
        if (oracle_) return count_impl(*oracle_, pattern, fsmp, bsmp);
        if (mirror_) return count_impl(*mirror_, pattern, fsmp, bsmp);
    }
    return count_impl(text_s_, pattern, nullptr, nullptr);
}

template <class TextT>
size_t LZ77Index::count_impl(const TextT& text_s, const std::string& pattern,
                             const uint8_t* fsmp, const uint8_t* bsmp) const {
    const size_t m = pattern.size();
    const auto* pat_u = reinterpret_cast<const unsigned char*>(pattern.data());

    size_t total = 0;

    // ── Special (end-aligned): frases que terminan con P completo ──────────────
    // La búsqueda binaria es exacta → todo punto del rango es real; el filtro
    // phrase_total_len >= m (en query_special) exige que P quepa en la frase.
    {
        size_t Ll, Lr;
        search_rev(text_s, psa_rev_, pat_u, m, Ll, Lr, bsmp);
        if (Lr >= Ll)
            total += grid_.query_special(Ll, Lr, m).count;
    }

    // ── Crossings: split P[0..i-1] | P[i..m-1] para i=1..m-1 ────────────────────
    for (size_t i = 1; i < m; ++i) {
        // Mitad derecha P[i..m-1] empieza en un boundary (psa_fwd).
        size_t Rl, Rr;
        search_sst(text_s, grid_, grid_.point_count(), pat_u + i, m - i, Rl, Rr,
                   fsmp, n_);
        if (Rr < Rl) continue;

        // Mitad izquierda P[0..i-1] termina en end_k (psa_rev).
        size_t Ll, Lr;
        search_rev(text_s, psa_rev_, pat_u, i, Ll, Lr, bsmp);
        if (Lr < Ll) continue;

        const auto res = grid_.query(Rl, Rr, Ll, Lr);
        if (res.first == 0) continue;

        // Cada punto con phrase_total_len >= i es una ocurrencia primaria real
        // (la ocurrencia empieza dentro de la frase k = su primer boundary cruzado).
        for (const auto& [wt_idx, y_rel] : res.second) {
            (void)y_rel;
            if (grid_.phrase_total_len(wt_idx) >= i) ++total;
        }
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
    if (m < 2 || grid_.point_count() == 0) return SIZE_MAX;
    if (!in_alphabet(pattern, alphabet_)) return SIZE_MAX;

    const auto* pat_u = reinterpret_cast<const unsigned char*>(pattern.data());

    // El accesor y el texto plano ofrecen la misma interfaz (size/operator[]);
    // la búsqueda es idéntica, solo cambia de dónde se leen los caracteres.
    const uint8_t* fsmp = (use_oracle_ && !fwd_smp_.empty()) ? fwd_smp_.data() : nullptr;
    const uint8_t* bsmp = (use_oracle_ && !bwd_smp_.empty()) ? bwd_smp_.data() : nullptr;
    if (use_oracle_) {
        if (oracle_) return locate_leftmost_impl(*oracle_, pattern, fsmp, bsmp);
        if (mirror_) return locate_leftmost_impl(*mirror_, pattern, fsmp, bsmp);
    }
    return locate_leftmost_impl(text_s_, pattern, nullptr, nullptr);
}

template <class TextT>
size_t LZ77Index::locate_leftmost_impl(const TextT& text_s,
                                       const std::string& pattern,
                                       const uint8_t* fsmp,
                                       const uint8_t* bsmp) const {
    const size_t m = pattern.size();
    const auto* pat_u = reinterpret_cast<const unsigned char*>(pattern.data());

    size_t pos_min = SIZE_MAX;

    // ── Special (end-aligned): P completo termina al final de una frase ─────────
    // Búsqueda binaria exacta → occ_min_pos es real, sin verificación.
    {
        size_t Ll, Lr;
        search_rev(text_s, psa_rev_, pat_u, m, Ll, Lr, bsmp);
        if (Lr >= Ll) {
            const auto sp = grid_.query_special(Ll, Lr, m);
            if (sp.count > 0 && sp.occ_min_pos != SIZE_MAX && sp.occ_min_pos < pos_min)
                pos_min = sp.occ_min_pos;
        }
    }

    // ── Crossings: split P[0..i-1] | P[i..m-1] para i=1..m-1 ────────────────────
    for (size_t i = 1; i < m; ++i) {
        size_t Rl, Rr;
        search_sst(text_s, grid_, grid_.point_count(), pat_u + i, m - i, Rl, Rr,
                   fsmp, n_);
        if (Rr < Rl) continue;

        size_t Ll, Lr;
        search_rev(text_s, psa_rev_, pat_u, i, Ll, Lr, bsmp);
        if (Lr < Ll) continue;

        const auto rmin = grid_.query_min(Rl, Rr, Ll, Lr);
        if (rmin.count == 0 || rmin.wt_idx == SIZE_MAX) continue;

        // El punto de mínimo boundary del RMQ. Búsqueda binaria exacta → cualquier
        // punto con phrase_total_len >= i es una ocurrencia real (empieza en la
        // frase k). Dos casos:
        size_t occ = SIZE_MAX;
        if (grid_.phrase_total_len(rmin.wt_idx) >= i && rmin.boundary_min >= i) {
            // (a) El extremo del RMQ ya cumple el filtro de largo de frase.
            occ = rmin.boundary_min - i;
        } else {
            // (b) El extremo del RMQ no cumple el filtro; puede haber otro punto
            //     (boundary mayor) que sí lo cumpla. Se busca el de mínimo boundary
            //     entre los válidos por largo.
            const auto mf = grid_.query_min_filtered(Rl, Rr, Ll, Lr, i);
            if (mf.count > 0 && mf.boundary_min >= i)
                occ = mf.boundary_min - i;
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
