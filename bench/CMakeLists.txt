# Micro-benchmarks for report 9 (docs/riskengine_research.md 8). Off by default
# (-DRISKENGINE_BUILD_BENCH=ON); run on a Release build, see bench/README.md.

include(FetchContent)
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
    benchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG v1.9.1
)
FetchContent_MakeAvailable(benchmark)

add_executable(riskengine_bench benchmarks.cpp)
target_link_libraries(riskengine_bench PRIVATE riskengine benchmark::benchmark)
# Record the build in the benchmark's JSON context, as the experiments do in their metadata.
if(TARGET riskengine_harness)
  target_link_libraries(riskengine_bench PRIVATE riskengine_harness)
  target_compile_definitions(riskengine_bench PRIVATE RISKENGINE_HAVE_HARNESS=1)
endif()
