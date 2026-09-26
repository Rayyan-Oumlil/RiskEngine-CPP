cmake_minimum_required(VERSION 3.25)
project(RiskEngineCPP LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()

option(RISKENGINE_ENABLE_SANITIZERS "Build with ASan/UBSan" OFF)
option(RISKENGINE_BUILD_TESTS "Build test suite" ON)
option(RISKENGINE_BUILD_BENCH "Build benchmarks" OFF)
option(RISKENGINE_BUILD_EXPERIMENTS "Build the report experiments" ON)
option(RISKENGINE_BUILD_EXAMPLES "Build the examples" ON)
option(RISKENGINE_BUILD_SERVER "Build the local dashboard server (fetches cpp-httplib, nlohmann_json)" OFF)
option(RISKENGINE_ENABLE_AVX2 "Vectorize random-normal generation with AVX2 (results are bit-identical either way)" OFF)
option(RISKENGINE_ENABLE_AVX512 "Vectorize random-normal generation with AVX-512F, 8 lanes (results are bit-identical either way)" OFF)

add_library(riskengine INTERFACE)
target_include_directories(riskengine INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_compile_features(riskengine INTERFACE cxx_std_20)

find_package(Threads REQUIRED)
target_link_libraries(riskengine INTERFACE Threads::Threads)

# Reproducibility (docs/model_risk_report.md 2.2):
# - never let the compiler fuse a*b+c into an FMA on its own, since that changes results between
#   compilers and target ISAs;
# - never let it evaluate exp/log/erfc/pow of compile-time constants itself: GCC does so with
#   correctly rounded MPFR arithmetic while the run-time libm may differ in the last ulp, so the
#   same expression would give different bits depending on what the optimizer inlines. (sqrt is
#   correctly rounded by IEEE either way and stays a builtin.)
if(MSVC)
  target_compile_options(riskengine INTERFACE /W4 /WX /permissive- /fp:precise)
else()
  target_compile_options(riskengine INTERFACE -Wall -Wextra -Werror -ffp-contract=off
    -fno-builtin-exp -fno-builtin-log -fno-builtin-erfc -fno-builtin-pow)
endif()

# AVX2 / AVX-512F only enable vector kernels that reproduce the scalar bits exactly (core/simd.hpp).
# -ffp-contract=off above forbids fusing into FMA whatever the target. AVX-512F implies AVX2.
if(RISKENGINE_ENABLE_AVX512)
  if(MSVC)
    target_compile_options(riskengine INTERFACE /arch:AVX512)
  else()
    target_compile_options(riskengine INTERFACE -mavx512f)
  endif()
elseif(RISKENGINE_ENABLE_AVX2)
  if(MSVC)
    target_compile_options(riskengine INTERFACE /arch:AVX2)
  else()
    target_compile_options(riskengine INTERFACE -mavx2)
  endif()
endif()

# MinGW GCC keeps the Windows stack only 16-byte aligned but may still spill 32/64-byte vector
# registers with aligned moves (GCC bug 54412): an intermittent crash, depending on where the stack
# happens to sit. Have the assembler emit unaligned moves instead, which cost nothing when aligned.
if(MINGW AND (RISKENGINE_ENABLE_AVX2 OR RISKENGINE_ENABLE_AVX512))
  target_compile_options(riskengine INTERFACE -Wa,-muse-unaligned-vector-move)
endif()

if(RISKENGINE_ENABLE_SANITIZERS AND NOT MSVC)
  target_compile_options(riskengine INTERFACE -fsanitize=address,undefined -fno-omit-frame-pointer)
  target_link_options(riskengine INTERFACE -fsanitize=address,undefined)
endif()

if(RISKENGINE_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

if(RISKENGINE_BUILD_EXPERIMENTS)
  add_subdirectory(experiments)
endif()

if(RISKENGINE_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

if(RISKENGINE_BUILD_BENCH)
  add_subdirectory(bench)
endif()

if(RISKENGINE_BUILD_SERVER)
  add_subdirectory(server)
endif()
