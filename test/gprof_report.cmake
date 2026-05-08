execute_process(
  COMMAND "${GPROF_EXECUTABLE}" -z -c "${NRG_EXECUTABLE}" "${GMON_SUM}"
  OUTPUT_FILE "${OUTPUT_FILE}"
  RESULT_VARIABLE result)

if(NOT result EQUAL 0)
  message(FATAL_ERROR "Failed to generate gprof report ${OUTPUT_FILE}")
endif()
