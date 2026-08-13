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
)
target_link_libraries(riskengine_tests PRIVATE riskengine Catch2::Catch2WithMain)

include(CTest)
include(Catch)
catch_discover_tests(riskengine_tests)
