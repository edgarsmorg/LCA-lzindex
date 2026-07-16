#pragma once

// Accesor al texto derivado del parsing LZ77: permite leer T sin almacenarlo.
// Reemplaza el texto plano (8n bits) por ~z(log n + log |dict|) bits.
//
// Mecanismo (el del LZ-index de Kreft-Navarro): cada frase k guarda el delta
// hacia su copia previa. Para leer T[p] se localiza la frase que contiene p; si
// p es el caracter explicito del final de la frase, se devuelve; si cae en la
// region copiada, se recurre a source_k + (p - start_k). La profundidad de esa
// recursion (h) crece muy lentamente con n (medido: 26/32/37 a 10/100/500 MB).
//
// Los deltas start_k - source_k se repiten masivamente en colecciones de copias
// (84 163 copias -> 10 491 deltas distintos a 100 MB), por lo que se codifican
// con un diccionario: log|dict| bits por frase en vez de log n.

#include "phrase.hpp"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <unordered_map>
#include <vector>

#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>

#include <iosfwd>

namespace lz77tax {

class TextOracle {
public:
    /// Construye el accesor desde el parsing. `syms` = simbolos presentes en T.
    void build(const LZ77Parsing& ph, size_t n, const std::vector<uint8_t>& syms) {
        n_ = n;
        z_ = ph.size();
        for (size_t i = 0; i < syms.size(); ++i) {
            map_[syms[i]]  = static_cast<uint8_t>(i);
            unmap_[i]      = syms[i];
        }
        const size_t sb = syms.size() > 1
                        ? sdsl::bits::hi(syms.size() - 1) + 1 : 1;

        start_ = sdsl::int_vector<>(z_, 0, sdsl::bits::hi(n) + 1);
        for (size_t k = 0; k < z_; ++k) start_[k] = ph[k].start_pos;

        // Diccionario de deltas: el codigo 0 marca "frase literal" (sin source).
        std::unordered_map<size_t, size_t> id;
        std::vector<size_t> dict;
        std::vector<size_t> code(z_, 0);
        for (size_t k = 0; k < z_; ++k) {
            if (ph[k].length == 0) continue;
            const size_t d = ph[k].start_pos - ph[k].source;
            auto it = id.find(d);
            if (it == id.end()) { id.emplace(d, dict.size()); dict.push_back(d); }
            code[k] = id[d] + 1;
        }
        dict_ = sdsl::int_vector<>(dict.size(), 0, sdsl::bits::hi(n) + 1);
        for (size_t i = 0; i < dict.size(); ++i) dict_[i] = dict[i];
        const size_t cw = dict.size() ? sdsl::bits::hi(dict.size()) + 1 : 1;
        code_ = sdsl::int_vector<>(z_, 0, cw);
        for (size_t k = 0; k < z_; ++k) code_[k] = code[k];
        endc_ = sdsl::int_vector<>(z_, 0, sb);
        for (size_t k = 0; k < z_; ++k) endc_[k] = map_[ph[k].next_char];
    }

    size_t size() const { return n_; }

    size_t size_in_bytes() const {
        // start_ NO se cuenta: es el conjunto de inicios de frase, identico a
        // grid.text_pos ordenado (ya serializado en la grilla). Se mantiene en
        // RAM pero es derivable de la grilla, igual que psa_rev_ → no se cobra.
        return sdsl::size_in_bytes(dict_)
             + sdsl::size_in_bytes(code_)  + sdsl::size_in_bytes(endc_);
    }

    // ── Serializacion (para el self-index: cargar sin el texto) ────────────────
    // NO serializa start_: se reconstruye desde la grilla en load() (inicios de
    // frase = grid.text_pos ordenado). unmap_ basta para extract (map_ solo se
    // usa al construir).
    size_t serialize(std::ostream& out) const {
        size_t w = 0;
        w += sdsl::write_member(n_, out);
        w += sdsl::write_member(z_, out);
        out.write(reinterpret_cast<const char*>(unmap_), 256); w += 256;
        w += sdsl::serialize(dict_, out);
        w += sdsl::serialize(code_, out);
        w += sdsl::serialize(endc_, out);
        return w;
    }

    /// Carga el accesor. `starts` = inicios de frase (derivados de la grilla).
    void load(std::istream& in, const std::vector<size_t>& starts) {
        sdsl::read_member(n_, in);
        sdsl::read_member(z_, in);
        in.read(reinterpret_cast<char*>(unmap_), 256);
        sdsl::load(dict_, in);
        sdsl::load(code_, in);
        sdsl::load(endc_, in);
        start_ = sdsl::int_vector<>(z_, 0, n_ > 0 ? sdsl::bits::hi(n_) + 1 : 1);
        for (size_t k = 0; k < z_ && k < starts.size(); ++k) start_[k] = starts[k];
    }

    /// Escribe T[p..p+len-1] en out. Recursion acotada por la profundidad h.
    void extract(size_t p, size_t len, uint8_t* out) const {
        while (len > 0) {
            const size_t k   = phrase_of(p);
            const size_t st  = start_[k];
            const size_t end = (k + 1 < z_) ? size_t(start_[k + 1]) - 1 : n_ - 1;
            if (p == end) { *out++ = unmap_[endc_[k]]; ++p; --len; continue; }
            const size_t take = std::min(len, end - p);
            extract(st - dict_[code_[k] - 1] + (p - st), take, out);
            out += take; p += take; len -= take;
        }
    }

    /// Acceso a un solo caracter (mas lento por caracter que extract en bloque).
    uint8_t operator[](size_t p) const { uint8_t c; extract(p, 1, &c); return c; }

private:
    /// Ultima frase con start <= p.
    size_t phrase_of(size_t p) const {
        size_t lo = 0, hi = z_ - 1;
        while (lo < hi) {
            const size_t mid = (lo + hi + 1) / 2;
            if (start_[mid] <= p) lo = mid; else hi = mid - 1;
        }
        return lo;
    }

    sdsl::int_vector<> start_, dict_, code_, endc_;
    size_t  n_ = 0, z_ = 0;
    uint8_t map_[256] = {0}, unmap_[256] = {0};
};

// El sub-indice sobre T^R no necesita accesor propio: leer T^R es leer T al
// reves. Con T = T[0..L-1] y centinela en L, ambos textos miden n = L+1 y
//     T^R_s[i] = T_s[L-1-i]   para i < L,   T^R_s[L] = '\0'
// de modo que un unico accesor sobre T sirve a los dos sub-indices y el espacio
// del accesor se reduce a la mitad.
class MirrorOracle {
public:
    MirrorOracle(const TextOracle& base, size_t n) : base_(&base), n_(n) {}

    size_t size() const { return n_; }

    void extract(size_t p, size_t len, uint8_t* out) const {
        const size_t L = n_ - 1;
        // Tramo espejado: posiciones p..min(p+len,L)-1 vienen de T[L-1-i].
        const size_t c = (p < L) ? std::min(len, L - p) : 0;
        if (c) {
            base_->extract(L - p - c, c, out);   // T[L-p-c .. L-1-p]
            std::reverse(out, out + c);          // out[j] = T[L-1-p-j]
        }
        // Centinela (i == L) y cualquier sobrante fuera de rango.
        for (size_t k = c; k < len; ++k) out[k] = 0;
    }

    uint8_t operator[](size_t p) const { uint8_t c; extract(p, 1, &c); return c; }

    /// El espejo no almacena nada propio: solo apunta al accesor directo.
    size_t size_in_bytes() const { return 0; }

private:
    const TextOracle* base_;
    size_t            n_;
};

}  // namespace lz77tax
