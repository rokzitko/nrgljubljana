# Regression tests

file(STRINGS simple_tests simple_tests)

foreach(test ${simple_tests})
  if(test STREQUAL "legacy_discretization_dmft_asym_U1-polarized")
    foreach(backend ${CHAIN_TEST_BACKENDS})
      add_chain_test(${test} ${backend} MODE runtime SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/${test}
        INPUTS param data PARAM_FILES param
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/runtest ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${test}
        LABELS chain-fixed-seed)
    endforeach()
    continue()
  endif()
  configure_file(${test}/param ${test}/param COPYONLY)
  configure_file(${test}/data ${test}/data COPYONLY)
  add_test(NAME ${test} COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/runtest ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${test}
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test})
  set_tests_properties(${test} PROPERTIES ENVIRONMENT "${TEST_BUILD_ENV}")
endforeach()
