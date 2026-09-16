// Report 4.3: delta, gamma and theta read off the nodes of an extended binomial tree (rooted two
// steps before today, Pelsser and Vorst 1994), against Black-Scholes for European options and
// against the same estimator on a fine tree (n = 2^15) for the American put, which has no closed
// form. CRR and BBS lattices, n and n + 1 for every n0 so that the even/odd effect is visible.

#include <cmath>

#include "harness/experiment.hpp"
#include "riskengine/methods/analytic/black_scholes.hpp"
#include "riskengine/methods/tree/binomial.hpp"

using namespace riskengine;

int main(int argc, char** argv) {
    return harness::run("tree_greeks", argc, argv, [](harness::Experiment& exp) {
        const MarketState market{Spot{100}, Rate{0.05}, Rate{0.0}, Vol{0.2}};
        const VanillaOption call{Strike{100}, Maturity{1}, OptionType::Call};
        const VanillaOption put{Strike{100}, Maturity{1}, OptionType::Put};
        constexpr unsigned kReferenceSteps = 1u << 15;
        const TreeGreeks american_ref = binomial_greeks(TreeMethod::Bbs, put, Exercise::American, market, kReferenceSteps);

        exp.param("market", "S = 100, r = 0.05, q = 0, sigma = 0.2, T = 1");
        exp.param("reference", "Black-Scholes (European); extended BBS tree at n = 32768 (American put)");
        exp.param("american_put_reference_delta", american_ref.delta);
        exp.param("american_put_reference_gamma", american_ref.gamma);
        exp.param("american_put_reference_theta", american_ref.theta);

        auto csv = exp.csv({"problem", "method", "n", "greek", "value", "reference", "relative_error"});
        struct Problem {
            const char* name;
            VanillaOption option;
            Exercise exercise;
        };
        for (const Problem& p : {Problem{"european_call_k100", call, Exercise::European},
                                 Problem{"american_put_k100", put, Exercise::American}}) {
            const bool european = p.exercise == Exercise::European;
            const Greeks bs = black_scholes_greeks(p.option, market);
            const double ref[3] = {european ? bs.delta : american_ref.delta, european ? bs.gamma : american_ref.gamma,
                                   european ? bs.theta : american_ref.theta};
            for (TreeMethod m : {TreeMethod::Crr, TreeMethod::Bbs})
                for (unsigned n0 : {25u, 50u, 100u, 200u, 500u, 1000u, 2000u, 5000u})
                    for (unsigned n : {n0, n0 + 1}) {
                        const TreeGreeks g = binomial_greeks(m, p.option, p.exercise, market, n);
                        const double value[3] = {g.delta, g.gamma, g.theta};
                        const char* names[3] = {"delta", "gamma", "theta"};
                        for (int i = 0; i < 3; ++i)
                            csv.row(p.name, to_string(m), n, names[i], value[i], ref[i], value[i] / ref[i] - 1.0);
                    }
        }
    });
}
