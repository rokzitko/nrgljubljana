# Regression tests

file(STRINGS simple_tests simple_tests)

foreach(test ${simple_tests})
  configure_file(${test}/param ${test}/param COPYONLY)
  configure_file(${test}/data ${test}/data COPYONLY)
  add_test(NAME ${test} COMMAND ${PROJECT_BINARY_DIR}/c++/nrg ; ${CMAKE_CURRENT_SOURCE_DIR}/compare
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test})
endforeach()
