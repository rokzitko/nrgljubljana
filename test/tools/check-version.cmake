foreach(version_option IN ITEMS -V --version)
  execute_process(
    COMMAND "${TOOL}" "${version_option}"
    RESULT_VARIABLE result
    OUTPUT_VARIABLE output
    ERROR_VARIABLE error
  )

  set(expected "${TOOL_NAME} ${EXPECTED_VERSION}\n")
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "${TOOL_NAME} ${version_option} exited with ${result}: ${error}")
  endif()
  if(NOT output STREQUAL expected)
    message(FATAL_ERROR "${TOOL_NAME} ${version_option} output was [${output}], expected [${expected}]")
  endif()
  if(NOT error STREQUAL "")
    message(FATAL_ERROR "${TOOL_NAME} ${version_option} wrote to stderr: ${error}")
  endif()

  execute_process(
    COMMAND "${TOOL}" -v "${version_option}" ignored-positional-argument
    RESULT_VARIABLE combined_result
    OUTPUT_VARIABLE combined_output
    ERROR_VARIABLE combined_error
  )
  if(NOT combined_result EQUAL 0 OR NOT combined_output STREQUAL expected OR NOT combined_error STREQUAL "")
    message(FATAL_ERROR "${TOOL_NAME} did not short-circuit a combined ${version_option} request")
  endif()
endforeach()
