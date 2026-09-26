# Python bindings: `import riskengine` (python/riskengine_py.cpp). Off by default
# (-DRISKENGINE_BUILD_PYTHON=ON). Needs Python >= 3.9 with its development headers; nanobind is
# fetched. The module lands in <build>/python, so `PYTHONPATH=<build>/python python ...` finds it.

find_package(Python 3.9 REQUIRED COMPONENTS Interpreter Development.Module)

include(FetchContent)
FetchContent_Declare(
    nanobind
    GIT_REPOSITORY https://github.com/wjakob/nanobind.git
    GIT_TAG v2.15.0
)
FetchContent_MakeAvailable(nanobind)

nanobind_add_module(riskengine_python NB_STATIC riskengine_py.cpp)
# $<1:...> keeps multi-config generators (MSVC) from adding a Release/ or Debug/ subdirectory, so
# the module is always at <build>/python.
set_target_properties(riskengine_python PROPERTIES
    OUTPUT_NAME riskengine
    LIBRARY_OUTPUT_DIRECTORY $<1:${CMAKE_CURRENT_BINARY_DIR}>)
target_link_libraries(riskengine_python PRIVATE riskengine)

# `ctest -R python` runs the pytest suite against the module just built.
if(RISKENGINE_BUILD_TESTS)
  add_test(NAME python_bindings
           COMMAND ${Python_EXECUTABLE} -m pytest -q ${CMAKE_CURRENT_SOURCE_DIR}/tests)
  set_tests_properties(python_bindings PROPERTIES
      ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}")
endif()
