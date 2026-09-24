# Regression tests

file(STRINGS all_tests all_tests)
file(STRINGS inst1_inputs INST1_INPUTS)
set(registered_tests)

foreach(test ${all_tests})
  foreach(backend IN LISTS CHAIN_TEST_BACKENDS)
    add_chain_test(${test} ${backend} MODE tool SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/${test}
      INPUTS param ${INST1_INPUTS}
      COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/runtest ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${test})
    list(APPEND registered_tests ${test}_${backend})
  endforeach()
endforeach()

if(GPROF)
  gprof("${registered_tests}" ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
endif()
