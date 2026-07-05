#include <gtest/gtest.h>

#include "wavelet/wm_rmq_min.hpp"
#include "test_helpers.hpp"

#include <algorithm>
#include <climits>
#include <numeric>
#include <random>
#include <vector>

#include <sdsl/bits.hpp>
#include <sdsl/construct.hpp>
#include <sdsl/int_vector.hpp>

using namespace lz77tax;

static sdsl::int_vector<> to_iv(const std::vector<size_t>& v) {
    sdsl::int_vector<> iv(v.size(), 0, 64);
    for (size_t i = 0; i < v.size(); ++i) iv[i] = v[i];
    return iv;
}

static SharedWm build_shared_wm(const std::vector<size_t>& ys, size_t sigma) {
    SharedWm wm;
    if (ys.empty() || sigma == 0) return wm;
    sdsl::int_vector<> yv(ys.size(), 0, sdsl::bits::hi(sigma - 1) + 1);
    for (size_t i = 0; i < ys.size(); ++i) yv[i] = ys[i];
    sdsl::construct_im(static_cast<sdsl::wm_int<>&>(wm), yv);
    return wm;
}

// brute_force_range_2d() viene de test_helpers.hpp (namespace lz77tax::test).
using lz77tax::test::brute_force_range_2d;

/// Verifica WmMinRmq contra brute-force en todos los rectángulos de n×sigma.
static void check_all_rectangles(const std::vector<size_t>& ys,
                                  const std::vector<size_t>& vals,
                                  size_t sigma) {
    const size_t n = ys.size();
    if (n == 0) return;

    const sdsl::int_vector<> vals_iv = to_iv(vals);
    SharedWm shared_wm = build_shared_wm(ys, sigma);
    WmMinRmq wm;
    wm.build(shared_wm, vals_iv);

    for (size_t x0 = 0; x0 < n; ++x0) {
        for (size_t x1 = x0; x1 < n; ++x1) {
            for (size_t y0 = 0; y0 < sigma; ++y0) {
                for (size_t y1 = y0; y1 < sigma; ++y1) {
                    const auto got = wm.range_argmin_2d(shared_wm, x0, x1, y0, y1, vals_iv);
                    const auto exp = brute_force_range_2d(ys, vals, x0, x1, y0, y1);

                    ASSERT_EQ(got.count, exp.count)
                        << "count mismatch x=[" << x0 << "," << x1
                        << "] y=[" << y0 << "," << y1 << "]";

                    if (exp.count > 0) {
                        ASSERT_EQ(vals[got.argmin_global], exp.argmin_val);
                    }
                }
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests de construcción
// ─────────────────────────────────────────────────────────────────────────────

TEST(WmMinRmq_Build, Empty) {
    WmMinRmq wm;
    sdsl::int_vector<> empty_iv;
    SharedWm shared_wm;
    wm.build(shared_wm, empty_iv);
    EXPECT_EQ(wm.size(), 0u);
    const auto r = wm.range_argmin_2d(shared_wm, 0, 0, 0, 0, empty_iv);
    EXPECT_EQ(r.count, 0u);
    EXPECT_EQ(r.argmin_global, SIZE_MAX);
}

TEST(WmMinRmq_Build, SingleElement) {
    WmMinRmq wm;
    const sdsl::int_vector<> vals_iv = to_iv({42});
    SharedWm shared_wm = build_shared_wm({0}, 1);
    wm.build(shared_wm, vals_iv);
    EXPECT_EQ(wm.size(), 1u);
    const auto r = wm.range_argmin_2d(shared_wm, 0, 0, 0, 0, vals_iv);
    EXPECT_EQ(r.count, 1u);
    EXPECT_EQ(r.argmin_global, 0u);
}

TEST(WmMinRmq_Build, TwoElements) {
    // Puntos: (x=0,y=1,v=10), (x=1,y=0,v=5)
    WmMinRmq wm;
    const std::vector<size_t> ys   = {1, 0};
    const std::vector<size_t> vals = {10, 5};
    const sdsl::int_vector<> vals_iv = to_iv(vals);
    SharedWm shared_wm = build_shared_wm(ys, 2);
    wm.build(shared_wm, vals_iv);
    EXPECT_EQ(wm.size(), 2u);

    const auto r = wm.range_argmin_2d(shared_wm, 0, 1, 0, 1, vals_iv);
    EXPECT_EQ(r.count, 2u);
    EXPECT_EQ(vals[r.argmin_global], 5u);

    const auto r2 = wm.range_argmin_2d(shared_wm, 0, 1, 0, 0, vals_iv);
    EXPECT_EQ(r2.count, 1u);
    EXPECT_EQ(r2.argmin_global, 1u);

    const auto r3 = wm.range_argmin_2d(shared_wm, 0, 1, 1, 1, vals_iv);
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
    SharedWm shared_wm = build_shared_wm(ys, 1);
    WmMinRmq wm;
    wm.build(shared_wm, vals_iv);

    for (size_t x0 = 0; x0 < n; ++x0) {
        for (size_t x1 = x0; x1 < n; ++x1) {
            const auto got = wm.range_argmin_2d(shared_wm, x0, x1, 0, 0, vals_iv);
            const auto exp = brute_force_range_2d(ys, vals, x0, x1, 0, 0);
            ASSERT_EQ(got.count, exp.count);
            if (exp.count > 0)
                ASSERT_EQ(vals[got.argmin_global], exp.argmin_val);
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
    SharedWm shared_wm = build_shared_wm(ys, n);
    WmMinRmq wm;
    wm.build(shared_wm, vals_iv);

    std::uniform_int_distribution<size_t> dis(0, n - 1);
    for (int trial = 0; trial < 1000; ++trial) {
        size_t x0 = dis(rng), x1 = dis(rng);
        size_t y0 = dis(rng), y1 = dis(rng);
        if (x0 > x1) std::swap(x0, x1);
        if (y0 > y1) std::swap(y0, y1);
        const auto got = wm.range_argmin_2d(shared_wm, x0, x1, y0, y1, vals_iv);
        const auto exp = brute_force_range_2d(ys, vals, x0, x1, y0, y1);
        ASSERT_EQ(got.count, exp.count) << "trial=" << trial;
        if (exp.count > 0)
            ASSERT_EQ(vals[got.argmin_global], exp.argmin_val) << "trial=" << trial;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Rectángulo vacío
// ─────────────────────────────────────────────────────────────────────────────

TEST(WmMinRmq_Query, EmptyRectangle) {
    const std::vector<size_t> ys   = {0, 1, 2};
    const std::vector<size_t> vals = {10, 20, 30};
    const sdsl::int_vector<> vals_iv = to_iv(vals);
    SharedWm shared_wm = build_shared_wm(ys, 3);
    WmMinRmq wm;
    wm.build(shared_wm, vals_iv);

    const auto r = wm.range_argmin_2d(shared_wm, 0, 2, 5, 10, vals_iv);
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
    SharedWm shared_wm = build_shared_wm(ys, n);
    WmMinRmq wm;
    wm.build(shared_wm, vals_iv);

    for (size_t x0 = 0; x0 < n; ++x0) {
        for (size_t x1 = x0; x1 < n; ++x1) {
            const auto r = wm.range_argmin_2d(shared_wm, x0, x1, 0, n - 1, vals_iv);
            if (r.count > 0) {
                ASSERT_LT(r.argmin_global, n) << "índice global fuera de [0,n)";
                ASSERT_GE(r.argmin_global, x0);
                ASSERT_LE(r.argmin_global, x1);
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Métricas de espacio: deben ser > 0 para n > 0
// ─────────────────────────────────────────────────────────────────────────────

TEST(WmMinRmq_Size, SizeBytesPositive) {
    const std::vector<size_t> ys   = {2, 0, 1, 3};
    const std::vector<size_t> vals = {5, 8, 2, 9};
    const sdsl::int_vector<> vals_iv = to_iv(vals);
    SharedWm shared_wm = build_shared_wm(ys, 4);
    WmMinRmq wm;
    wm.build(shared_wm, vals_iv);

    EXPECT_GT(wm.size_in_bytes(), 0u);
    const auto bd = wm.size_breakdown();
    EXPECT_EQ(bd.bv,  0u);
    EXPECT_EQ(bd.rank_sel, 0u);
    EXPECT_GT(bd.rmq, 0u);
}
