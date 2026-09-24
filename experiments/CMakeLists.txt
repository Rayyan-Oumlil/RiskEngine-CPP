# One executable per figure or table of docs/model_risk_report.md. Each writes
# data/results/<id>.csv and <id>.meta.json (see harness/experiment.hpp).

set(git_info_header "${CMAKE_CURRENT_BINARY_DIR}/generated/riskengine/git_info.hpp")
add_custom_target(riskengine_git_info
  COMMAND ${CMAKE_COMMAND}
    -DSOURCE_DIR=${PROJECT_SOURCE_DIR}
    -DINPUT_FILE=${CMAKE_CURRENT_SOURCE_DIR}/harness/git_info.hpp.in
    -DOUTPUT_FILE=${git_info_header}
    -P ${PROJECT_SOURCE_DIR}/cmake/GitInfo.cmake
  BYPRODUCTS ${git_info_header}
  COMMENT "Recording git commit for experiment metadata")

add_library(riskengine_harness INTERFACE)
add_dependencies(riskengine_harness riskengine_git_info)
target_include_directories(riskengine_harness INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}/generated)
# Generator expressions, so multi-config generators (MSVC) record the configuration actually built.
# The riskengine target's own options (-ffp-contract=off, sanitizers) are recorded too.
set(config_flags "")
foreach(config IN ITEMS Debug Release RelWithDebInfo MinSizeRel)
  string(TOUPPER ${config} config_upper)
  string(APPEND config_flags "$<$<CONFIG:${config}>:${CMAKE_CXX_FLAGS_${config_upper}}>")
endforeach()
target_compile_definitions(riskengine_harness INTERFACE
  "RISKENGINE_COMPILER=\"${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}\""
  "RISKENGINE_BUILD_TYPE=\"$<CONFIG>\""
  "RISKENGINE_CXX_FLAGS=\"${CMAKE_CXX_FLAGS} ${config_flags} $<JOIN:$<TARGET_PROPERTY:riskengine,INTERFACE_COMPILE_OPTIONS>, >\"")
target_link_libraries(riskengine_harness INTERFACE riskengine)

set(RISKENGINE_EXPERIMENTS
  fd_vcurve
  iv_roundtrip
  mc_convergence
  mc_coverage
  mc_efficiency
  qmc_convergence
  greeks_vs_h
  greeks_vs_n
  greeks_matrix
  var_straddle
  var_backtest
  stress_scenarios
  tree_convergence
  tree_greeks
  heston_discretization
  model_risk_exotics
)

foreach(experiment IN LISTS RISKENGINE_EXPERIMENTS)
  add_executable(${experiment} ${experiment}.cpp)
  target_link_libraries(${experiment} PRIVATE riskengine_harness)
endforeach()

# `cmake --build <dir> --target report`: runs every experiment from the repository root (some read
# data/raw) and renders every figure, i.e. regenerates data/results and docs/figures. Results meant
# for the repository must come from a Release build of a clean tree (report Appendix D).
find_package(Python3 COMPONENTS Interpreter)
set(report_commands "")
foreach(experiment IN LISTS RISKENGINE_EXPERIMENTS)
  list(APPEND report_commands COMMAND $<TARGET_FILE:${experiment}>)
endforeach()
if(Python3_Interpreter_FOUND)
  list(APPEND report_commands COMMAND ${Python3_EXECUTABLE} tools/make_figures.py)
else()
  list(APPEND report_commands COMMAND ${CMAKE_COMMAND} -E echo
       "Python 3 not found: results regenerated, figures not rendered (run tools/make_figures.py)")
endif()
add_custom_target(report ${report_commands}
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  DEPENDS ${RISKENGINE_EXPERIMENTS}
  COMMENT "Regenerating the report's results and figures"
  USES_TERMINAL
  VERBATIM)
