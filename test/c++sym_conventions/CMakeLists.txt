# Cross-symmetry convention tests: the same physical problem solved with different symmetry types must give the same
# expectation values and FDM thermodynamics (shared references). See README.md.

file(STRINGS simple_tests simple_tests)

foreach(test ${simple_tests})
  configure_file(${test}/param ${test}/param COPYONLY)
  configure_file(${test}/data ${test}/data COPYONLY)
  add_test(NAME ${test} COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/runtest ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${test}
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test})
  set_tests_properties(${test} PROPERTIES ENVIRONMENT "${TEST_BUILD_ENV}" LABELS conventions)
endforeach()
