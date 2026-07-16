#include <gtest/gtest.h>

#include "lz77/grid.hpp"

#include <algorithm>
#include <numeric>
#include <random>
#include <vector>

using namespace lz77tax;

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Grilla aleatoria de z1 puntos: R es una permutación de {0..z1-1} (invariante
/// real de la grilla: cada frase aporta exactamente una coordenada Y), los
/// boundaries son posiciones de texto crecientes y los largos de frase varían.
struct RandomGrid {
    std::vector<size_t> R, boundaries, phrase_lens;
};

static RandomGrid make_random_grid(size_t z1, std::mt19937_64& rng) {
    RandomGrid g;
    g.R.resize(z1);
    std::iota(g.R.begin(), g.R.end(), 0);
    std::shuffle(g.R.begin(), g.R.end(), rng);

    g.boundaries.resize(z1);
    g.phrase_lens.resize(z1);
    size_t pos = 1;
    for (size_t j = 0; j < z1; ++j) {
        const size_t len = 1 + rng() % 8;
        pos += len;
        g.boundaries[j]  = pos;
        g.phrase_lens[j] = len;
    }
    // Los boundaries se generaron crecientes en j; barajarlos rompe la
    // correlación entre la coordenada X y la posición de texto, que es lo que
    // hace no trivial al RMQ.
    std::shuffle(g.boundaries.begin(), g.boundaries.end(), rng);
    return g;
}

/// Mínimo boundary por enumeración directa de los puntos del rectángulo.
static Grid2D::MinResult brute_min(const RandomGrid& g,
                                   size_t x_lo, size_t x_hi,
                                   size_t y_lo, size_t y_hi,
                                   size_t min_phrase_len = 0) {
    size_t count = 0, best = SIZE_MAX, best_idx = SIZE_MAX;
    for (size_t x = x_lo; x <= x_hi && x < g.R.size(); ++x) {
        if (g.R[x] < y_lo || g.R[x] > y_hi) continue;
        if (g.phrase_lens[x] < min_phrase_len) continue;
        ++count;
        if (g.boundaries[x] < best) { best = g.boundaries[x]; best_idx = x; }
    }
    if (count == 0) return {0, SIZE_MAX, SIZE_MAX};
    return {count, best, best_idx};
}

// ─────────────────────────────────────────────────────────────────────────────
// Construcción
// ─────────────────────────────────────────────────────────────────────────────

TEST(Grid2D_Build, PointCountMatchesInput) {
    Grid2D grid;
    grid.build(/*R=*/{2, 0, 3, 1}, /*boundaries=*/{40, 10, 50, 20},
               /*phrase_lens=*/{3, 1, 4, 2});
    EXPECT_EQ(grid.point_count(), 4u);
}

TEST(Grid2D_Build, EmptyGrid) {
    Grid2D grid;
    grid.build({}, {}, {});
    EXPECT_EQ(grid.point_count(), 0u);

    const auto mn = grid.query_min(0, 10, 0, 10);
    EXPECT_EQ(mn.count, 0u);
    EXPECT_EQ(mn.boundary_min, SIZE_MAX);
}

TEST(Grid2D_Build, AccessorsReturnStoredValues) {
    Grid2D grid;
    grid.build({1, 0}, {30, 10}, {5, 2});
    EXPECT_EQ(grid.text_pos(0), 30u);
    EXPECT_EQ(grid.text_pos(1), 10u);
    EXPECT_EQ(grid.phrase_total_len(0), 5u);
    EXPECT_EQ(grid.phrase_total_len(1), 2u);
}

// ─────────────────────────────────────────────────────────────────────────────
// Consultas degeneradas
// ─────────────────────────────────────────────────────────────────────────────

TEST(Grid2D_Query, EmptyOrInvalidRanges) {
    Grid2D grid;
    grid.build({2, 0, 3, 1}, {40, 10, 50, 20}, {3, 1, 4, 2});

    // Rango invertido → vacío
    EXPECT_EQ(grid.query(3, 1, 0, 3).first, 0u);
    EXPECT_EQ(grid.query_min(3, 1, 0, 3).count, 0u);
    // Rango fuera de la grilla → vacío, sin desbordar
    EXPECT_EQ(grid.query(0, 99, 0, 3).first, 0u);
    EXPECT_EQ(grid.query_min(0, 3, 0, 99).count, 0u);
}

// ─────────────────────────────────────────────────────────────────────────────
// Invariante central: query_min == mínimo por enumeración
// ─────────────────────────────────────────────────────────────────────────────

TEST(Grid2D_MinRmq, MinMatchesBruteForce_Randomized) {
    std::mt19937_64 rng(20260715);
    for (size_t z1 : {1u, 2u, 3u, 5u, 8u, 16u, 17u, 64u, 100u}) {
        const auto g = make_random_grid(z1, rng);
        Grid2D grid;
        grid.build(g.R, g.boundaries, g.phrase_lens);
        ASSERT_EQ(grid.point_count(), z1);

        for (int rep = 0; rep < 200; ++rep) {
            size_t x_lo = rng() % z1, x_hi = rng() % z1;
            size_t y_lo = rng() % z1, y_hi = rng() % z1;
            if (x_lo > x_hi) std::swap(x_lo, x_hi);
            if (y_lo > y_hi) std::swap(y_lo, y_hi);

            const auto expected = brute_min(g, x_lo, x_hi, y_lo, y_hi);
            const auto actual   = grid.query_min(x_lo, x_hi, y_lo, y_hi);

            ASSERT_EQ(actual.count, expected.count)
                << "z1=" << z1 << " x=[" << x_lo << "," << x_hi
                << "] y=[" << y_lo << "," << y_hi << "]";
            if (expected.count > 0) {
                ASSERT_EQ(actual.boundary_min, expected.boundary_min)
                    << "z1=" << z1 << " x=[" << x_lo << "," << x_hi
                    << "] y=[" << y_lo << "," << y_hi << "]";
            }
        }
    }
}

TEST(Grid2D_Query, RangeSearchMatchesBruteForce_Randomized) {
    std::mt19937_64 rng(987654321);
    for (size_t z1 : {1u, 4u, 16u, 33u}) {
        const auto g = make_random_grid(z1, rng);
        Grid2D grid;
        grid.build(g.R, g.boundaries, g.phrase_lens);

        for (int rep = 0; rep < 100; ++rep) {
            size_t x_lo = rng() % z1, x_hi = rng() % z1;
            size_t y_lo = rng() % z1, y_hi = rng() % z1;
            if (x_lo > x_hi) std::swap(x_lo, x_hi);
            if (y_lo > y_hi) std::swap(y_lo, y_hi);

            const auto expected = brute_min(g, x_lo, x_hi, y_lo, y_hi);
            const auto [count, pts] = grid.query(x_lo, x_hi, y_lo, y_hi);
            ASSERT_EQ(count, expected.count);
            // Todo punto devuelto debe caer realmente en el rectángulo.
            for (const auto& [x, y] : pts) {
                EXPECT_GE(x, x_lo); EXPECT_LE(x, x_hi);
                EXPECT_GE(y, y_lo); EXPECT_LE(y, y_hi);
                EXPECT_EQ(g.R[x], y);
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Filtro por largo de frase
// ─────────────────────────────────────────────────────────────────────────────

TEST(Grid2D_FilteredExtremal, MinSkipsShortPhraseExtremal) {
    Grid2D grid;
    grid.build(/*R=*/{0, 1}, /*boundaries=*/{10, 20}, /*phrase_lens=*/{1, 5});

    const auto fast = grid.query_min(0, 1, 0, 1);
    ASSERT_EQ(fast.count, 2u);
    EXPECT_EQ(fast.boundary_min, 10u);
    EXPECT_LT(grid.phrase_total_len(fast.wt_idx), 3u);

    const auto filtered = grid.query_min_filtered(0, 1, 0, 1, 3);
    ASSERT_EQ(filtered.count, 1u);
    EXPECT_EQ(filtered.boundary_min, 20u);
    EXPECT_GE(grid.phrase_total_len(filtered.wt_idx), 3u);
}

TEST(Grid2D_FilteredExtremal, MatchesBruteForce_Randomized) {
    std::mt19937_64 rng(555);
    for (size_t z1 : {4u, 16u, 40u}) {
        const auto g = make_random_grid(z1, rng);
        Grid2D grid;
        grid.build(g.R, g.boundaries, g.phrase_lens);

        for (int rep = 0; rep < 100; ++rep) {
            size_t x_lo = rng() % z1, x_hi = rng() % z1;
            size_t y_lo = rng() % z1, y_hi = rng() % z1;
            if (x_lo > x_hi) std::swap(x_lo, x_hi);
            if (y_lo > y_hi) std::swap(y_lo, y_hi);
            const size_t min_len = rng() % 10;

            const auto expected = brute_min(g, x_lo, x_hi, y_lo, y_hi, min_len);
            const auto actual   = grid.query_min_filtered(x_lo, x_hi, y_lo, y_hi, min_len);
            ASSERT_EQ(actual.count, expected.count);
            if (expected.count > 0)
                ASSERT_EQ(actual.boundary_min, expected.boundary_min);
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Ocurrencias end-aligned
// ─────────────────────────────────────────────────────────────────────────────

TEST(Grid2D_Special, MatchesBruteForce_Randomized) {
    std::mt19937_64 rng(4242);
    for (size_t z1 : {4u, 16u, 40u}) {
        const auto g = make_random_grid(z1, rng);
        Grid2D grid;
        grid.build(g.R, g.boundaries, g.phrase_lens);

        for (int rep = 0; rep < 100; ++rep) {
            size_t y_lo = rng() % z1, y_hi = rng() % z1;
            if (y_lo > y_hi) std::swap(y_lo, y_hi);
            const size_t plen = 1 + rng() % 6;

            size_t count = 0, occ_min = SIZE_MAX, occ_max = 0;
            for (size_t x = 0; x < z1; ++x) {
                if (g.R[x] < y_lo || g.R[x] > y_hi) continue;
                if (g.phrase_lens[x] < plen) continue;
                ++count;
                const size_t occ = g.boundaries[x] - plen;
                occ_min = std::min(occ_min, occ);
                occ_max = std::max(occ_max, occ);
            }

            const auto sp = grid.query_special(y_lo, y_hi, plen);
            ASSERT_EQ(sp.count, count) << "y=[" << y_lo << "," << y_hi
                                       << "] plen=" << plen;
            if (count > 0) {
                EXPECT_EQ(sp.occ_min_pos, occ_min);
                EXPECT_EQ(sp.occ_max_pos, occ_max);
            }
        }
    }
}
