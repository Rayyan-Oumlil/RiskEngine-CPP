include(FetchContent)
FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.7.1
)
FetchContent_MakeAvailable(Catch2)

add_executable(riskengine_tests
    test_smoke.cpp
    test_black_scholes.cpp
    test_implied_vol.cpp
    test_analytic_exotics.cpp
    test_rng.cpp
    test_sobol.cpp
    test_simulation.cpp
    test_monte_carlo.cpp
    test_brownian_bridge.cpp
    test_dual.cpp
    test_mc_greeks.cpp
    test_var.cpp
    test_time_series.cpp
    test_portfolio_risk.cpp
    test_backtest.cpp
    test_trees.cpp
    test_models.cpp
)
target_link_libraries(riskengine_tests PRIVATE riskengine Catch2::Catch2WithMain)
target_compile_definitions(riskengine_tests PRIVATE RISKENGINE_DATA_DIR="${PROJECT_SOURCE_DIR}/data/raw")

include(CTest)
include(Catch)
catch_discover_tests(riskengine_tests)
