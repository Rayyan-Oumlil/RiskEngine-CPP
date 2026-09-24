# Writes the current git commit and dirty state into a header. Runs at build time (see
# experiments/CMakeLists.txt) so the recorded commit is never stale; configure_file only touches
# the output when its content changes, so an unchanged commit does not trigger rebuilds. "Dirty"
# includes untracked, non-ignored files: an unadded source file also breaks reproducibility.
execute_process(
  COMMAND git rev-parse HEAD
  WORKING_DIRECTORY "${SOURCE_DIR}"
  OUTPUT_VARIABLE RISKENGINE_GIT_COMMIT
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
  RESULT_VARIABLE git_result)
if(NOT git_result EQUAL 0)
  set(RISKENGINE_GIT_COMMIT "unknown")
endif()

execute_process(
  COMMAND git status --porcelain
  WORKING_DIRECTORY "${SOURCE_DIR}"
  OUTPUT_VARIABLE git_status
  ERROR_QUIET)
if(git_status STREQUAL "")
  set(RISKENGINE_GIT_DIRTY false)
else()
  set(RISKENGINE_GIT_DIRTY true)
endif()

configure_file("${INPUT_FILE}" "${OUTPUT_FILE}" @ONLY)
