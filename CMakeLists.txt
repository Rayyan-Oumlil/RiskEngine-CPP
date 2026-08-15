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

add_library(riskengine INTERFACE)
target_include_directories(riskengine INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_compile_features(riskengine INTERFACE cxx_std_20)

find_package(Threads REQUIRED)
target_link_libraries(riskengine INTERFACE Threads::Threads)

# Reproducibility: never let the compiler fuse a*b+c into an FMA on its own, since that changes
# results between compilers and target ISAs (docs/riskengine_research.md 10).
if(MSVC)
  target_compile_options(riskengine INTERFACE /W4 /WX /permissive- /fp:precise)
else()
  target_compile_options(riskengine INTERFACE -Wall -Wextra -Werror -ffp-contract=off)
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

if(RISKENGINE_BUILD_BENCH)
  add_subdirectory(bench)
endif()
