#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <vector>

#include "riskengine/core/rng/random_stream.hpp"
#include "riskengine/risk/var.hpp"

using namespace riskengine;

TEST_CASE("Normal VaR and ES reproduce the reference table", "[var]") {
    // Analytic normal VaR/ES table, checked independently against SciPy (docs/model_risk_report.md,
    // Appendix C).
    struct Row {
        double alpha, var, es;
    };
    for (const Row& row : {Row{0.95, 1.6449, 2.0627}, Row{0.975, 1.9600, 2.3378}, Row{0.99, 2.3263, 2.6652}}) {
        const RiskMeasures m = normal_var_es(0.0, 1.0, row.alpha);
        CHECK(std::abs(m.var - row.var) < 5e-5);
        CHECK(std::abs(m.es - row.es) < 5e-5);
    }
    // FRTB's point: ES 97.5 % sits just above VaR 99 % under normality.
    CHECK(normal_var_es(0, 1, 0.975).es > normal_var_es(0, 1, 0.99).var);
    const RiskMeasures shifted = normal_var_es(2.0, 3.0, 0.99);
    CHECK(std::abs(shifted.var - (2.0 + 3.0 * 2.3263478740408408)) < 1e-12);
}

TEST_CASE("Empirical VaR and ES converge to the normal values", "[var]") {
    std::vector<double> losses(1'000'000);
    RandomStream rng(SeedKey{31}, 0);
    for (double& l : losses) l = rng.normal();
    for (double alpha : {0.95, 0.99}) {
        const RiskMeasures e = empirical_var_es(losses, alpha);
        const RiskMeasures n = normal_var_es(0, 1, alpha);
        // Quantile standard error sqrt(alpha (1 - alpha) / N) / n(q): about 0.004 at 99 %.
        CHECK(std::abs(e.var - n.var) < 0.02);
        CHECK(std::abs(e.es - n.es) < 0.02);
    }
}

TEST_CASE("Empirical conventions: index ceil(alpha N) and the tail mean from it upwards", "[var]") {
    std::vector<double> losses;
    for (int i = 1; i <= 100; ++i) losses.push_back(double(i)); // losses 1..100, shuffled below
    std::vector<double> shuffled;
    for (int i = 0; i < 100; ++i) shuffled.push_back(losses[(i * 37) % 100]);
    const RiskMeasures m = empirical_var_es(shuffled, 0.95);
    CHECK(m.var == 95.0);                                       // the 95th smallest loss
    CHECK(m.es == (95.0 + 96 + 97 + 98 + 99 + 100) / 6.0);      // mean of the losses >= VaR
}

TEST_CASE("VaR is not subadditive; ES is", "[var][subadditivity]") {
    // Two independent bonds, each defaulting with probability 4 % for a loss of 100 (Artzner et al.
    // 1999). The distributions are represented exactly by 10,000 equally weighted scenarios.
    std::vector<double> a, b, both;
    for (int i = 0; i < 100; ++i)
        for (int j = 0; j < 100; ++j) {
            const double la = i < 4 ? 100.0 : 0.0, lb = j < 4 ? 100.0 : 0.0;
            a.push_back(la);
            b.push_back(lb);
            both.push_back(la + lb);
        }
    const RiskMeasures ra = empirical_var_es(a, 0.95), rb = empirical_var_es(b, 0.95), rab = empirical_var_es(both, 0.95);
    CHECK(ra.var == 0.0);
    CHECK(rb.var == 0.0);
    CHECK(rab.var == 100.0); // diversification *increases* VaR: 100 > 0 + 0
    CHECK(rab.es <= ra.es + rb.es);
}

TEST_CASE("Cornish-Fisher reduces to the normal quantile without skew or excess kurtosis", "[var]") {
    for (double alpha : {0.9, 0.99, 0.999}) CHECK(cornish_fisher_quantile(alpha, 0, 0) == norm_icdf(alpha));
    // Positive skew fattens the right tail: a higher upper quantile.
    CHECK(cornish_fisher_quantile(0.99, 1.0, 0.0) > norm_icdf(0.99));
    CHECK(cornish_fisher_quantile(0.99, 0.0, 3.0) > norm_icdf(0.99));
}
