#include <gtest/gtest.h>

#include "wavelet/wm_rmq_min.hpp"
#include "wavelet/wt_rmq_min.hpp"  // oracle para comparación side-by-side

#include <algorithm>
#include <climits>
#include <numeric>
#include <random>
#include <vector>

#include <sdsl/int_vector.hpp>

using namespace lz77tax;

static sdsl::int_vector<> to_iv(const std::vector<size_t>& v) {
    sdsl::int_vector<> iv(v.size(), 0, 64);
    for (size_t i = 0; i < v.size(); ++i) iv[i] = v[i];
    return iv;
}

// ─────────────────────────────────────────────────────────────────────────────
// Brute-force oracle
// ─────────────────────────────────────────────────────────────────────────────

struct BfResult { size_t count; size_t argmin; };

static BfResult brute_force(const std::vector<size_t>& ys,
                             const std::vector<size_t>& vals,
                             size_t x_lo, size_t x_hi,
                             size_t y_lo, size_t y_hi) {
    size_t count = 0, best_v = SIZE_MAX, best_i = SIZE_MAX;
    for (size_t x = x_lo; x <= x_hi; ++x) {
        if (ys[x] < y_lo || ys[x] > y_hi) continue;
        count++;
        if (vals[x] < best_v) { best_v = vals[x]; best_i = x; }
    }
    return {count, best_i};
}

/// Verifica WmMinRmq contra brute-force en todos los rectángulos de n×sigma.
static void check_all_rectangles(const std::vector<size_t>& ys,
                                  const std::vector<size_t>& vals,
                                  size_t sigma) {
    const size_t n = ys.size();
    if (n == 0) return;

    const sdsl::int_vector<> vals_iv = to_iv(vals);
    WmMinRmq wm;
    wm.build(ys, vals_iv, sigma);

    for (size_t x0 = 0; x0 < n; ++x0) {
        for (size_t x1 = x0; x1 < n; ++x1) {
            for (size_t y0 = 0; y0 < sigma; ++y0) {
                for (size_t y1 = y0; y1 < sigma; ++y1) {
                    const auto got = wm.range_argmin_2d(x0, x1, y0, y1, vals_iv);
                    const auto exp = brute_force(ys, vals, x0, x1, y0, y1);

                    ASSERT_EQ(got.count, exp.count)
                        << "count mismatch x=[" << x0 << "," << x1
                        << "] y=[" << y0 << "," << y1 << "]";

                    if (exp.count > 0) {
                        ASSERT_EQ(vals[got.argmin_global], vals[exp.argmin]);
                    }
                }
            }
        }
    }
}

/// Verifica WmMinRmq == WtMinRmq (misma respuesta, no necesariamente mismo índice en empate).
static void check_matches_wt(const std::vector<size_t>& ys,
                               const std::vector<size_t>& vals,
                               size_t sigma,
                               size_t n_queries = 1000,
                               uint64_t seed = 42) {
    const size_t n = ys.size();
    const sdsl::int_vector<> vals_iv = to_iv(vals);

    WmMinRmq wm;
    WtMinRmq wt;
    wm.build(ys, vals_iv, sigma);
    wt.build(ys, vals_iv, sigma);

    std::mt19937 rng(seed);
    std::uniform_int_distribution<size_t> xd(0, n - 1);
    std::uniform_int_distribution<size_t> yd(0, sigma - 1);

    for (size_t trial = 0; trial < n_queries; ++trial) {
        size_t x0 = xd(rng), x1 = xd(rng);
        size_t y0 = yd(rng), y1 = yd(rng);
        if (x0 > x1) std::swap(x0, x1);
        if (y0 > y1) std::swap(y0, y1);

        const auto r_wm = wm.range_argmin_2d(x0, x1, y0, y1, vals_iv);
        const auto r_wt = wt.range_argmin_2d(x0, x1, y0, y1, vals_iv);

        ASSERT_EQ(r_wm.count, r_wt.count) << "trial=" << trial;
        if (r_wt.count > 0) {
            ASSERT_EQ(vals[r_wm.argmin_global], vals[r_wt.argmin_global])
                << "trial=" << trial;
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests de construcción
// ─────────────────────────────────────────────────────────────────────────────

TEST(WmMinRmq_Build, Empty) {
    WmMinRmq wm;
    sdsl::int_vector<> empty_iv;
    wm.build({}, empty_iv, 0);
    EXPECT_EQ(wm.size(), 0u);
    const auto r = wm.range_argmin_2d(0, 0, 0, 0, empty_iv);
    EXPECT_EQ(r.count, 0u);
    EXPECT_EQ(r.argmin_global, SIZE_MAX);
}

TEST(WmMinRmq_Build, SingleElement) {
    WmMinRmq wm;
    const sdsl::int_vector<> vals_iv = to_iv({42});
    wm.build({0}, vals_iv, 1);
    EXPECT_EQ(wm.size(), 1u);
    const auto r = wm.range_argmin_2d(0, 0, 0, 0, vals_iv);
    EXPECT_EQ(r.count, 1u);
    EXPECT_EQ(r.argmin_global, 0u);
}

TEST(WmMinRmq_Build, TwoElements) {
    // Puntos: (x=0,y=1,v=10), (x=1,y=0,v=5)
    WmMinRmq wm;
    const std::vector<size_t> ys   = {1, 0};
    const std::vector<size_t> vals = {10, 5};
    const sdsl::int_vector<> vals_iv = to_iv(vals);
    wm.build(ys, vals_iv, 2);
    EXPECT_EQ(wm.size(), 2u);

    const auto r = wm.range_argmin_2d(0, 1, 0, 1, vals_iv);
    EXPECT_EQ(r.count, 2u);
    EXPECT_EQ(vals[r.argmin_global], 5u);

    const auto r2 = wm.range_argmin_2d(0, 1, 0, 0, vals_iv);
    EXPECT_EQ(r2.count, 1u);
    EXPECT_EQ(r2.argmin_global, 1u);

    const auto r3 = wm.range_argmin_2d(0, 1, 1, 1, vals_iv);
    EXPECT_EQ(r3.count, 1u);
    EXPECT_EQ(r3.argmin_global, 0u);
}

// ─────────────────────────────────────────────────────────────────────────────
// Argmin vs brute-force — exhaustivo en n pequeño
// ─────────────────────────────────────────────────────────────────────────────

TEST(WmMinRmq_Query, ArgminMatchesBrute_SmallFixed) {
    const std::vector<size_t> ys   = {3, 1, 6, 0, 5, 2, 7, 4};
    const std::vector<size_t> vals = {100, 20, 80, 50, 10, 60, 30, 70};
    check_all_rectangles(ys, vals, 8);
}

TEST(WmMinRmq_Query, ArgminMatchesBrute_AllSameY) {
    const size_t n = 6;
    const std::vector<size_t> ys(n, 0);
    const std::vector<size_t> vals = {5, 3, 8, 1, 7, 2};
    const sdsl::int_vector<> vals_iv = to_iv(vals);
    WmMinRmq wm;
    wm.build(ys, vals_iv, 1);

    for (size_t x0 = 0; x0 < n; ++x0) {
        for (size_t x1 = x0; x1 < n; ++x1) {
            const auto got = wm.range_argmin_2d(x0, x1, 0, 0, vals_iv);
            const auto exp = brute_force(ys, vals, x0, x1, 0, 0);
            ASSERT_EQ(got.count, exp.count);
            if (exp.count > 0)
                ASSERT_EQ(vals[got.argmin_global], vals[exp.argmin]);
        }
    }
}

TEST(WmMinRmq_Query, ArgminMatchesBrute_Random_Medium) {
    // n=200, σ=200, permutación aleatoria seed fijo — 1000 rectángulos aleatorios
    const size_t n = 200;
    std::vector<size_t> ys(n), vals(n);
    std::iota(ys.begin(),   ys.end(),   0);
    std::iota(vals.begin(), vals.end(), 0);
    std::mt19937 rng(12345);
    std::shuffle(ys.begin(),   ys.end(),   rng);
    std::shuffle(vals.begin(), vals.end(), rng);

    const sdsl::int_vector<> vals_iv = to_iv(vals);
    WmMinRmq wm;
    wm.build(ys, vals_iv, n);

    std::uniform_int_distribution<size_t> dis(0, n - 1);
    for (int trial = 0; trial < 1000; ++trial) {
        size_t x0 = dis(rng), x1 = dis(rng);
        size_t y0 = dis(rng), y1 = dis(rng);
        if (x0 > x1) std::swap(x0, x1);
        if (y0 > y1) std::swap(y0, y1);
        const auto got = wm.range_argmin_2d(x0, x1, y0, y1, vals_iv);
        const auto exp = brute_force(ys, vals, x0, x1, y0, y1);
        ASSERT_EQ(got.count, exp.count) << "trial=" << trial;
        if (exp.count > 0)
            ASSERT_EQ(vals[got.argmin_global], vals[exp.argmin]) << "trial=" << trial;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Rectángulo vacío
// ─────────────────────────────────────────────────────────────────────────────

TEST(WmMinRmq_Query, EmptyRectangle) {
    const std::vector<size_t> ys   = {0, 1, 2};
    const std::vector<size_t> vals = {10, 20, 30};
    const sdsl::int_vector<> vals_iv = to_iv(vals);
    WmMinRmq wm;
    wm.build(ys, vals_iv, 3);

    const auto r = wm.range_argmin_2d(0, 2, 5, 10, vals_iv);
    EXPECT_EQ(r.count, 0u);
    EXPECT_EQ(r.argmin_global, SIZE_MAX);
}

// ─────────────────────────────────────────────────────────────────────────────
// Validez del índice global
// ─────────────────────────────────────────────────────────────────────────────

TEST(WmMinRmq_Query, ArgminGlobalIsValidIndex) {
    const size_t n = 16;
    std::vector<size_t> ys(n), vals(n);
    std::iota(ys.begin(),   ys.end(),   0);
    std::iota(vals.begin(), vals.end(), 0);
    std::mt19937 rng(99999);
    std::shuffle(ys.begin(),   ys.end(),   rng);
    std::shuffle(vals.begin(), vals.end(), rng);

    const sdsl::int_vector<> vals_iv = to_iv(vals);
    WmMinRmq wm;
    wm.build(ys, vals_iv, n);

    for (size_t x0 = 0; x0 < n; ++x0) {
        for (size_t x1 = x0; x1 < n; ++x1) {
            const auto r = wm.range_argmin_2d(x0, x1, 0, n - 1, vals_iv);
            if (r.count > 0) {
                ASSERT_LT(r.argmin_global, n) << "índice global fuera de [0,n)";
                ASSERT_GE(r.argmin_global, x0);
                ASSERT_LE(r.argmin_global, x1);
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Comparación WmMinRmq == WtMinRmq (side-by-side oracle)
// ─────────────────────────────────────────────────────────────────────────────

TEST(WmMinRmq_VsWt, SmallRandom) {
    const size_t n = 50;
    std::vector<size_t> ys(n), vals(n);
    std::iota(ys.begin(),   ys.end(),   0);
    std::iota(vals.begin(), vals.end(), 0);
    std::mt19937 rng(7);
    std::shuffle(ys.begin(),   ys.end(),   rng);
    std::shuffle(vals.begin(), vals.end(), rng);
    check_matches_wt(ys, vals, n, /*n_queries=*/2000, 7);
}

TEST(WmMinRmq_VsWt, MediumRandom) {
    const size_t n = 1000;
    std::vector<size_t> ys(n), vals(n);
    std::iota(ys.begin(),   ys.end(),   0);
    std::iota(vals.begin(), vals.end(), 1000);
    std::mt19937 rng(2024);
    std::shuffle(ys.begin(),   ys.end(),   rng);
    std::shuffle(vals.begin(), vals.end(), rng);
    check_matches_wt(ys, vals, n, /*n_queries=*/5000, 2024);
}

// ─────────────────────────────────────────────────────────────────────────────
// Métricas de espacio: deben ser > 0 para n > 0
// ─────────────────────────────────────────────────────────────────────────────

TEST(WmMinRmq_Size, SizeBytesPositive) {
    const std::vector<size_t> ys   = {2, 0, 1, 3};
    const std::vector<size_t> vals = {5, 8, 2, 9};
    const sdsl::int_vector<> vals_iv = to_iv(vals);
    WmMinRmq wm;
    wm.build(ys, vals_iv, 4);

    EXPECT_GT(wm.size_in_bytes(), 0u);
    const auto bd = wm.size_breakdown();
    EXPECT_GT(bd.bv,  0u);
    EXPECT_GT(bd.rmq, 0u);
}
