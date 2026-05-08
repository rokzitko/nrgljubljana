foreach(test IN LISTS SIMPLE_TESTS)
  set(gmon_out "${CURRENT_BINARY_DIR}/${test}/gmon.out")
  if(NOT EXISTS "${gmon_out}")
    continue()
  endif()

  if(EXISTS "${GMON_SUM}")
    execute_process(
      COMMAND "${GPROF_EXECUTABLE}" "${NRG_EXECUTABLE}" -s -z -c "${GMON_SUM}" "${gmon_out}"
      WORKING_DIRECTORY "${BINARY_DIR}"
      RESULT_VARIABLE result)
  else()
    execute_process(
      COMMAND "${CMAKE_COMMAND}" -E copy "${gmon_out}" "${GMON_SUM}"
      RESULT_VARIABLE result)
  endif()

  if(NOT result EQUAL 0)
    message(FATAL_ERROR "Failed to collect gprof output from ${gmon_out}")
  endif()
endforeach()
