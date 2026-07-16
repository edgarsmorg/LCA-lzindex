#include "parser.hpp"

#include <algorithm>
#include <vector>

#include <sdsl/divsufsort.hpp>

namespace lz77tax {

// ── Kasai's algorithm para LCP array ─────────────────────────────────────────
// lcp[i] = longitud del prefijo común más largo entre los sufijos SA[i-1] y SA[i].
// lcp[0] = 0 por convención.
static void build_lcp_kasai(const std::string& text,
                             const std::vector<int32_t>& sa,
                             const std::vector<int32_t>& isa,
                             std::vector<int32_t>& lcp) {
    const size_t n = text.size();
    lcp.assign(n, 0);
    int32_t h = 0;
    for (size_t i = 0; i < n; i++) {
        if (isa[i] > 0) {
            int32_t j = sa[isa[i] - 1];
            while (static_cast<size_t>(i + h) < n &&
                   static_cast<size_t>(j + h) < n &&
                   text[i + h] == text[j + h])
                ++h;
            lcp[isa[i]] = h;
            if (h > 0) --h;
        }
    }
}

// ── Parser LZ77 greedy de izquierda a derecha ─────────────────────────────────
//
// Para cada posición i buscamos el factor previo más largo: el mayor len tal que
// T[i..i+len-1] aparece en T[0..i-1]. La búsqueda usa SA + LCP:
//
//  - En el SA, T[i..] está en el rank r = ISA[i].
//  - Sus vecinos con SA[k] < i son las ocurrencias previas.
//  - El LCP entre el rank k y r es min(LCP[k+1..r]) (o LCP[r+1..k] si k > r).
//  - Escaneamos izquierda y derecha desde r, manteniendo un mínimo corriente.
//
// Por qué parar en el primer SA[k] < i encontrado:
//   El mínimo corriente min(LCP[k+1..r]) es NO CRECIENTE al alejar k de r.
//   Por tanto la primera fuente válida hallada en cada dirección garantiza el
//   mayor LCP posible en esa dirección; seguir escaneando solo puede empeorar.
//
// Complejidad: O(n log n) para construir SA con divsufsort + O(n) para LCP
// (Kasai). El greedy tiene z frases con un scan por frase; en texto repetitivo
// z ≪ n y el scan termina rápido (min_lcp cae a 0), pero el peor caso es
// O(n·z) ≈ O(n²). Para los corpus objetivo (z ~ n/log n) es eficiente en práctica.
//
// Límite: usa int32_t para SA/ISA/LCP → texto máximo 2^31 - 1 ≈ 2 GB.
//
LZ77Parsing LZ77Parser::parse(const std::string& text) {
    if (text.empty()) {
        phrase_count_ = 0;
        return {};
    }

    const size_t n = text.size();

    // 1. Suffix Array con divsufsort (header-only en sdsl duscob fork)
    std::vector<int32_t> sa(n);
    sdsl::divsufsort(reinterpret_cast<const uint8_t*>(text.data()),
                     sa.data(), static_cast<int32_t>(n));

    // 2. Inverse Suffix Array
    std::vector<int32_t> isa(n);
    for (size_t i = 0; i < n; i++)
        isa[sa[i]] = static_cast<int32_t>(i);

    // 3. LCP array (Kasai)
    std::vector<int32_t> lcp(n, 0);
    build_lcp_kasai(text, sa, isa, lcp);

    // 4. Greedy LZ77 parsing
    LZ77Parsing phrases;
    size_t i = 0;

    while (i < n) {
        const size_t r = static_cast<size_t>(isa[i]);
        size_t best_len = 0;
        size_t best_src = 0;

        // Escaneo izquierdo: desde rank r-1 bajando a 0.
        // Invariante: al entrar en iteración k, min_lcp = min(LCP[k+1..r])
        //             = LCP(SA[k], SA[r]) = largo del match entre T[SA[k]..] y T[i..].
        // Al bajar k, min_lcp es no-creciente → la primera fuente SA[k]<i da el
        // mejor match posible hacia la izquierda; no se necesita seguir.
        if (r > 0) {
            size_t min_lcp = static_cast<size_t>(lcp[r]);  // LCP(SA[r-1], SA[r])
            for (int32_t k = static_cast<int32_t>(r) - 1; k >= 0; --k) {
                if (min_lcp == 0) break;
                if (static_cast<size_t>(sa[k]) < i) {
                    if (min_lcp > best_len) {
                        best_len = min_lcp;
                        best_src = static_cast<size_t>(sa[k]);
                    }
                    break;
                }
                // Prepara min_lcp para la siguiente iteración (k-1):
                // LCP(SA[k-1], SA[r]) = min(LCP[k..r]) = min(min_lcp, lcp[k])
                if (k > 0)
                    min_lcp = std::min(min_lcp, static_cast<size_t>(lcp[k]));
                else
                    break;  // k=0 sin fuente válida: no hay más candidatos
            }
        }

        // Escaneo derecho: desde rank r+1 subiendo a n-1.
        // Invariante: al entrar en iteración k, min_lcp = min(LCP[r+1..k])
        //             = LCP(SA[r], SA[k]) = largo del match entre T[i..] y T[SA[k]..].
        if (r + 1 < n) {
            size_t min_lcp = static_cast<size_t>(lcp[r + 1]);  // LCP(SA[r], SA[r+1])
            for (size_t k = r + 1; k < n; ++k) {
                if (min_lcp == 0) break;
                if (static_cast<size_t>(sa[k]) < i) {
                    if (min_lcp > best_len) {
                        best_len = min_lcp;
                        best_src = static_cast<size_t>(sa[k]);
                    }
                    break;
                }
                // Prepara min_lcp para la siguiente iteración (k+1):
                // LCP(SA[r], SA[k+1]) = min(LCP[r+1..k+1]) = min(min_lcp, lcp[k+1])
                if (k + 1 < n)
                    min_lcp = std::min(min_lcp, static_cast<size_t>(lcp[k + 1]));
                else
                    break;  // k=n-1 sin fuente válida: no hay más candidatos
            }
        }

        // Garantizamos que siempre existe un next_char (T[i+len] en [0..n-1]).
        if (best_len > 0 && i + best_len >= n)
            best_len = n - i - 1;

        if (best_len == 0) {
            // Frase literal: solo el carácter T[i]
            phrases.push_back({i, 0, static_cast<uint8_t>(text[i]), 0});
            i += 1;
        } else {
            // Frase de copia: copia best_len chars + agrega T[i+best_len]
            phrases.push_back({i, best_len,
                               static_cast<uint8_t>(text[i + best_len]),
                               best_src});
            i += best_len + 1;
        }
    }

    phrase_count_ = phrases.size();
    return phrases;
}

}  // namespace lz77tax
