find_package(Perl REQUIRED)

# Only the two nrgchain comparisons select a scalar backend; mixchain's own method is unchanged.

file(STRINGS simple_tests simple_tests)

foreach(test ${simple_tests})
  if(test STREQUAL "mixchain1_flat" OR test STREQUAL "mixchain2_nrgchain")
    foreach(backend IN LISTS CHAIN_TEST_BACKENDS)
      add_chain_test(${test} ${backend} MODE tool SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/${test}
        INPUTS param param-nrgchain PARAM_FILES param-nrgchain
        COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/${test}/run ${PROJECT_BINARY_DIR}/tools ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${test})
    endforeach()
    continue()
  endif()
  configure_file(${test}/run ${test}/run COPYONLY) # also ensures the test dir is created
  configure_file(${test}/param ${test}/param COPYONLY)
  add_test(NAME ${test} COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/${test}/run ${PROJECT_BINARY_DIR}/tools ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${test}
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test})
  set_tests_properties(${test} PROPERTIES ENVIRONMENT "${TEST_BUILD_ENV}")
endforeach()
