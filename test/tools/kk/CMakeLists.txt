find_package(Perl REQUIRED)

# Regression tests

file(STRINGS simple_tests simple_tests)

foreach(test ${simple_tests})
  configure_file(${test}/run ${test}/run COPYONLY) # also ensures the test dir is created
  add_test(NAME ${test} COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/${test}/run ${PROJECT_BINARY_DIR}/tools ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${test}
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test})
endforeach()

configure_file(kk1/imaw.dat kk1/imaw.dat COPYONLY)
