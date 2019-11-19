find_package(Perl REQUIRED)

# Regression tests

file(STRINGS simple_tests simple_tests)

foreach(test ${simple_tests})
  configure_file(${test}/param ${test}/param COPYONLY)
  configure_file(${test}/matrix ${test}/matrix COPYONLY)
  file(GLOB dat_files "${test}/*.dat")
  file(COPY ${dat_files} DESTINATION ${test})
  add_test(NAME ${test} COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/${test}/run ${PROJECT_BINARY_DIR}/tools; ${CMAKE_CURRENT_SOURCE_DIR}/compare
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test})
endforeach()
