find_package(Perl REQUIRED)

# Regression tests

file(STRINGS simple_tests simple_tests)

foreach(test ${simple_tests})
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${test}/spec.bin)
    configure_file(${test}/spec.bin ${test}/spec.bin COPYONLY)
  endif()
  add_test(NAME ${test} COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/${test}/run ${PROJECT_BINARY_DIR}/tools ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${test}
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test})
endforeach()
