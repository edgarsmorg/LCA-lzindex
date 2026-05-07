#include <gtest/gtest.h>

#include "wavelet/wt_rmq_min.hpp"

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
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Brute-force argmin 2D: escanea todos los puntos en el rectángulo.
static WtMinRmq::ArgminResult brute_force_argmin(
        const std::vector<size_t>& ys,
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

/// Construye un WtMinRmq y lo llama; compara count y valor min con brute-force.
static void check_all_rectangles(const std::vector<size_t>& ys,
                                  const std::vector<size_t>& vals,
                                  size_t sigma) {
    const size_t n = ys.size();
    if (n == 0) return;

    const sdsl::int_vector<> vals_iv = to_iv(vals);
    WtMinRmq wt;
    wt.build(ys, vals_iv, sigma);

    for (size_t x0 = 0; x0 < n; ++x0) {
        for (size_t x1 = x0; x1 < n; ++x1) {
            for (size_t y0 = 0; y0 < sigma; ++y0) {
                for (size_t y1 = y0; y1 < sigma; ++y1) {
                    const auto got = wt.range_argmin_2d(x0, x1, y0, y1, vals_iv);
                    const auto exp = brute_force_argmin(ys, vals, x0, x1, y0, y1);

                    ASSERT_EQ(got.count, exp.count)
                        << "count mismatch x=[" << x0 << "," << x1
                        << "] y=[" << y0 << "," << y1 << "]";

                    if (exp.count > 0) {
                        // En empate el RMQ puede elegir cualquier índice;
                        // comparamos el VALOR mínimo, no el índice.
                        ASSERT_EQ(vals[got.argmin_global], vals[exp.argmin_global])
                            << "min_value mismatch x=[" << x0 << "," << x1
                            << "] y=[" << y0 << "," << y1 << "]";
                    }
                }
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests de construcción
// ─────────────────────────────────────────────────────────────────────────────

TEST(WtMinRmq_Build, Empty) {
    WtMinRmq wt;
    sdsl::int_vector<> empty_iv;
    wt.build({}, empty_iv, 0);
    EXPECT_EQ(wt.size(), 0u);
    const auto r = wt.range_argmin_2d(0, 0, 0, 0, empty_iv);
    EXPECT_EQ(r.count, 0u);
    EXPECT_EQ(r.argmin_global, SIZE_MAX);
}

TEST(WtMinRmq_Build, SingleElement) {
    WtMinRmq wt;
    const sdsl::int_vector<> vals_iv = to_iv({42});
    wt.build({0}, vals_iv, 1);
    EXPECT_EQ(wt.size(), 1u);

    const auto r = wt.range_argmin_2d(0, 0, 0, 0, vals_iv);
    EXPECT_EQ(r.count, 1u);
    EXPECT_EQ(r.argmin_global, 0u);
}

TEST(WtMinRmq_Build, TwoElements) {
    // Puntos: (x=0,y=1,v=10), (x=1,y=0,v=5)
    // ys = {1, 0}, vals = {10, 5}, sigma = 2
    WtMinRmq wt;
    const std::vector<size_t> ys   = {1, 0};
    const std::vector<size_t> vals = {10, 5};
    const sdsl::int_vector<> vals_iv = to_iv(vals);
    wt.build(ys, vals_iv, 2);
    EXPECT_EQ(wt.size(), 2u);

    // Todo el rectángulo: argmin = índice 1 (val=5)
    const auto r = wt.range_argmin_2d(0, 1, 0, 1, vals_iv);
    EXPECT_EQ(r.count, 2u);
    EXPECT_EQ(vals[r.argmin_global], 5u);

    // Solo y=0: solo el punto x=1
    const auto r2 = wt.range_argmin_2d(0, 1, 0, 0, vals_iv);
    EXPECT_EQ(r2.count, 1u);
    EXPECT_EQ(r2.argmin_global, 1u);

    // Solo y=1: solo el punto x=0
    const auto r3 = wt.range_argmin_2d(0, 1, 1, 1, vals_iv);
    EXPECT_EQ(r3.count, 1u);
    EXPECT_EQ(r3.argmin_global, 0u);
}

// ─────────────────────────────────────────────────────────────────────────────
// Argmin vs brute-force — exhaustivo en n pequeño
// ─────────────────────────────────────────────────────────────────────────────

TEST(WtMinRmq_Query, ArgminMatchesBrute_SmallFixed) {
    // n=8, σ=8, permutación de ys y vals arbitrarios
    const std::vector<size_t> ys   = {3, 1, 6, 0, 5, 2, 7, 4};
    const std::vector<size_t> vals = {100, 20, 80, 50, 10, 60, 30, 70};
    check_all_rectangles(ys, vals, 8);
}

TEST(WtMinRmq_Query, ArgminMatchesBrute_AllSameY) {
    // Todos los puntos con y=0 — solo un nodo hoja
    const size_t n = 6;
    const std::vector<size_t> ys(n, 0);
    const std::vector<size_t> vals = {5, 3, 8, 1, 7, 2};
    const sdsl::int_vector<> vals_iv = to_iv(vals);
    WtMinRmq wt;
    wt.build(ys, vals_iv, 1);

    for (size_t x0 = 0; x0 < n; ++x0) {
        for (size_t x1 = x0; x1 < n; ++x1) {
            const auto got = wt.range_argmin_2d(x0, x1, 0, 0, vals_iv);
            const auto exp = brute_force_argmin(ys, vals, x0, x1, 0, 0);
            ASSERT_EQ(got.count, exp.count);
            if (exp.count > 0) {
                ASSERT_EQ(vals[got.argmin_global], vals[exp.argmin_global]);
            }
        }
    }
}

TEST(WtMinRmq_Query, ArgminMatchesBrute_Random_Medium) {
    // n=200, σ=200, permutación aleatoria seed fijo
    const size_t n = 200;
    std::vector<size_t> ys(n), vals(n);
    std::iota(ys.begin(),   ys.end(),   0);
    std::iota(vals.begin(), vals.end(), 0);
    std::mt19937 rng(12345);
    std::shuffle(ys.begin(),   ys.end(),   rng);
    std::shuffle(vals.begin(), vals.end(), rng);

    const sdsl::int_vector<> vals_iv = to_iv(vals);
    WtMinRmq wt;
    wt.build(ys, vals_iv, n);

    // 500 rectángulos aleatorios
    std::uniform_int_distribution<size_t> dis(0, n - 1);
    for (int trial = 0; trial < 500; ++trial) {
        size_t x0 = dis(rng), x1 = dis(rng);
        size_t y0 = dis(rng), y1 = dis(rng);
        if (x0 > x1) std::swap(x0, x1);
        if (y0 > y1) std::swap(y0, y1);

        const auto got = wt.range_argmin_2d(x0, x1, y0, y1, vals_iv);
        const auto exp = brute_force_argmin(ys, vals, x0, x1, y0, y1);

        ASSERT_EQ(got.count, exp.count) << "trial=" << trial;
        if (exp.count > 0) {
            ASSERT_EQ(vals[got.argmin_global], vals[exp.argmin_global])
                << "trial=" << trial;
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Rectángulo vacío
// ─────────────────────────────────────────────────────────────────────────────

TEST(WtMinRmq_Query, EmptyRectangle) {
    const std::vector<size_t> ys   = {0, 1, 2};
    const std::vector<size_t> vals = {10, 20, 30};
    const sdsl::int_vector<> vals_iv = to_iv(vals);
    WtMinRmq wt;
    wt.build(ys, vals_iv, 3);

    // y_lo > y_hi implícito via rango vacío: y fuera de rango
    const auto r = wt.range_argmin_2d(0, 2, 5, 10, vals_iv);
    EXPECT_EQ(r.count, 0u);
    EXPECT_EQ(r.argmin_global, SIZE_MAX);
}

// ─────────────────────────────────────────────────────────────────────────────
// Unwind: verificar que argmin_global es realmente un índice válido en raíz
// ─────────────────────────────────────────────────────────────────────────────

TEST(WtMinRmq_Query, ArgminGlobalIsValidIndex) {
    const size_t n = 16;
    std::vector<size_t> ys(n), vals(n);
    std::iota(ys.begin(),   ys.end(),   0);
    std::iota(vals.begin(), vals.end(), 0);
    std::mt19937 rng(99999);
    std::shuffle(ys.begin(),   ys.end(),   rng);
    std::shuffle(vals.begin(), vals.end(), rng);

    const sdsl::int_vector<> vals_iv = to_iv(vals);
    WtMinRmq wt;
    wt.build(ys, vals_iv, n);

    for (size_t x0 = 0; x0 < n; ++x0) {
        for (size_t x1 = x0; x1 < n; ++x1) {
            const auto r = wt.range_argmin_2d(x0, x1, 0, n - 1, vals_iv);
            if (r.count > 0) {
                ASSERT_LT(r.argmin_global, n) << "índice global fuera de [0,n)";
                // el índice global debe estar en el rango X consultado
                ASSERT_GE(r.argmin_global, x0);
                ASSERT_LE(r.argmin_global, x1);
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Métricas de espacio: deben ser > 0 para n > 0
// ─────────────────────────────────────────────────────────────────────────────

TEST(WtMinRmq_Size, SizeBytesPositive) {
    const std::vector<size_t> ys   = {2, 0, 1, 3};
    const std::vector<size_t> vals = {5, 8, 2, 9};
    const sdsl::int_vector<> vals_iv = to_iv(vals);
    WtMinRmq wt;
    wt.build(ys, vals_iv, 4);

    EXPECT_GT(wt.size_in_bytes(), 0u);
    const auto bd = wt.size_breakdown();
    EXPECT_GT(bd.rmq, 0u);  // siempre hay RMQs
}
