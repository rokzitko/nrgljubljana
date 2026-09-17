find_package(Perl REQUIRED)

# Regression tests. Each case copies only its run script and parameter file; the run script takes any further
# input from its source directory, which it receives as the third argument.

file(STRINGS simple_tests simple_tests)

foreach(test ${simple_tests})
  configure_file(${test}/run ${test}/run COPYONLY) # also ensures the test dir is created
  configure_file(${test}/param ${test}/param COPYONLY)
  add_test(NAME ${test} COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/${test}/run ${PROJECT_BINARY_DIR}/tools ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${test}
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test})
  set_tests_properties(${test} PROPERTIES ENVIRONMENT "${TEST_BUILD_ENV}")
endforeach()
