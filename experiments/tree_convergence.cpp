// Report 4.1-4.2: convergence of binomial trees, on consecutive n and on their error envelope.
//
// Three problems at S = 100, r = 5 %, q = 0, sigma = 20 %, T = 1: an at-the-money European call
// (a node falls on the strike for even n, the strike sits midway between nodes for odd n), a
// European call struck at 110 (the strike's position between nodes moves irregularly with n) and an
// at-the-money American put (no closed form: the reference is BBS-Richardson at n = 2^16, checked
// against Leisen-Reimer at n = 2^16 + 1).
//
// tree_convergence.csv: price and error of every method for every valid n from 10 to 400, which
// is what exposes the oscillation (sampling only even n hides it). tree_convergence_envelope.csv:
// for n0 = 50 ... 10,000, the largest absolute error over the 8 consecutive valid n from n0, an
// envelope that does not depend on where an oscillating error happens to cross zero; the orders of
// convergence are fitted on it.

#include <cmath>

#include "harness/experiment.hpp"
#include "riskengine/methods/analytic/black_scholes.hpp"
#include "riskengine/methods/tree/binomial.hpp"

using namespace riskengine;

namespace {

const MarketState kMarket{Spot{100}, Rate{0.05}, Rate{0.0}, Vol{0.2}};
constexpr TreeMethod kMethods[] = {TreeMethod::Crr, TreeMethod::CrrAveraged, TreeMethod::LeisenReimer,
                                   TreeMethod::Bbs, TreeMethod::BbsRichardson};

bool valid(TreeMethod m, unsigned n) {
    if (m == TreeMethod::LeisenReimer) return n % 2 == 1;
    if (m == TreeMethod::BbsRichardson) return n % 2 == 0;
    return true;
}

struct Problem {
    const char* name;
    VanillaOption option;
    Exercise exercise;
    double reference;
};

} // namespace

int main(int argc, char** argv) {
    return harness::run("tree_convergence", argc, argv, [](harness::Experiment& exp) {
        const VanillaOption atm_call{Strike{100}, Maturity{1}, OptionType::Call};
        const VanillaOption otm_call{Strike{110}, Maturity{1}, OptionType::Call};
        const VanillaOption atm_put{Strike{100}, Maturity{1}, OptionType::Put};
        constexpr unsigned kReferenceSteps = 1u << 16;
        const double american_reference =
            binomial_price(TreeMethod::BbsRichardson, atm_put, Exercise::American, kMarket, kReferenceSteps);
        const double american_check =
            binomial_price(TreeMethod::LeisenReimer, atm_put, Exercise::American, kMarket, kReferenceSteps + 1);
        const Problem problems[] = {
            {"european_call_k100", atm_call, Exercise::European, black_scholes_price(atm_call, kMarket)},
            {"european_call_k110", otm_call, Exercise::European, black_scholes_price(otm_call, kMarket)},
            {"american_put_k100", atm_put, Exercise::American, american_reference},
        };

        exp.param("market", "S = 100, r = 0.05, q = 0, sigma = 0.2, T = 1");
        exp.param("reference", "Black-Scholes for the European calls; BBS-Richardson at n = 65536 for the American put");
        exp.param("american_put_reference", american_reference);
        exp.param("american_put_check_leisen_reimer_65537", american_check);
        exp.param("consecutive_n", "10 .. 400, every n valid for the method (Leisen-Reimer odd, BBS-Richardson even)");

        auto csv = exp.csv({"problem", "method", "n", "price", "error"});
        for (const Problem& p : problems)
            for (TreeMethod m : kMethods)
                for (unsigned n = 10; n <= 400; ++n) {
                    if (!valid(m, n)) continue;
                    const double v = binomial_price(m, p.option, p.exercise, kMarket, n);
                    csv.row(p.name, to_string(m), n, v, v - p.reference);
                }

        harness::Experiment envelope("tree_convergence_envelope", exp.out_dir());
        envelope.param("source", "tree_convergence (same run, same references)");
        envelope.param("envelope", "max |error| over the 8 consecutive valid n starting at n0");
        auto env = envelope.csv({"problem", "method", "n0", "max_abs_error", "n_at_max"});
        for (const Problem& p : problems)
            for (TreeMethod m : kMethods)
                for (unsigned n0 : {50u, 100u, 200u, 500u, 1000u, 2000u, 5000u, 10000u}) {
                    double worst = -1.0;
                    unsigned at = 0;
                    for (unsigned n = n0, taken = 0; taken < 8; ++n) {
                        if (!valid(m, n)) continue;
                        const double e = std::abs(binomial_price(m, p.option, p.exercise, kMarket, n) - p.reference);
                        if (e > worst) {
                            worst = e;
                            at = n;
                        }
                        ++taken;
                    }
                    env.row(p.name, to_string(m), n0, worst, at);
                }
        envelope.write_metadata();
    });
}
