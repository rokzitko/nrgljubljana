find_package(Perl REQUIRED)

# Regression tests

file(STRINGS simple_tests simple_tests)

foreach(test ${simple_tests})
  configure_file(${test}/param ${test}/param COPYONLY)
  add_test(NAME ${test} COMMAND ${PROJECT_SOURCE_DIR}/nrginit/nrginit ${Mathematica_KERNEL_EXECUTABLE}; ${CMAKE_CURRENT_SOURCE_DIR}/compare
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test})
endforeach()
