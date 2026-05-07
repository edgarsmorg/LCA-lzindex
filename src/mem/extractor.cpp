#include "extractor.hpp"

#include <vector>

namespace lz77tax {

MEMList MEMExtractor::extract(const std::string& query, size_t min_length) const {
    const size_t m = query.size();
    const size_t n = csa_.size();
    if (m == 0 || min_length == 0) return {};

    // lb[j] = menor i tal que Q[i..j] ocurre en T.
    //         m+1 (sentinel) si Q[j] como carácter aislado no ocurre en T.
    std::vector<size_t> lb(m, m + 1);
    std::vector<size_t> sps(m, 0);  // intervalo SA [sps[j], eps[j]] = ocurrencias de Q[lb[j]..j]

    for (size_t j = 0; j < m; ++j) {
        size_t sp = 0, ep = n - 1;
        for (size_t k = j + 1; k-- > 0; ) {
            size_t new_sp, new_ep;
            sdsl::backward_search(csa_, sp, ep,
                                  static_cast<unsigned char>(query[k]),
                                  new_sp, new_ep);
            if (new_sp > new_ep) {
                lb[j] = k + 1;   // Q[k..j] no ocurre; longest match es Q[k+1..j]
                break;
            }
            sp = new_sp;
            ep = new_ep;
            if (k == 0) lb[j] = 0;  // Q[0..j] ocurre completo
        }
        sps[j] = sp;
    }

    MEMList result;
    for (size_t j = 0; j < m; ++j) {
        const size_t i = lb[j];
        if (i > j) continue;  // ningún carácter Q[j] existe en T

        const size_t len = j - i + 1;
        if (len < min_length) continue;

        // Right-maximal: Q[i..j+1] no ocurre en T.
        // lb[j+1] > i  ⟺  el match más largo terminando en j+1 empieza después de i
        //               ⟺  Q[i..j+1] no ocurre.
        const bool right_max = (j == m - 1) || (lb[j + 1] > i);
        if (!right_max) continue;

        result.push_back({i, csa_[sps[j]], len});
    }
    return result;
}

}  // namespace lz77tax
