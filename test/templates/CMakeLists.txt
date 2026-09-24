# Regression tests

file(STRINGS simple_tests simple_tests)
file(STRINGS template50_inputs TEMPLATE50_INPUTS)
set(profile_tests)

foreach(test ${simple_tests})
  if(test STREQUAL "template49_ref")
    # Prepared data is backend-neutral; no reconstruction takes place.
    configure_file(${test}/param ${test}/param COPYONLY)
    configure_file(${test}/data ${test}/data COPYONLY)
    add_test(NAME ${test} COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/runtest ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${test} ${Mathematica_KERNEL_EXECUTABLE}
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test})
    set_tests_properties(${test} PROPERTIES ENVIRONMENT "${TEST_BUILD_ENV}")
    list(APPEND profile_tests ${test})
  else()
    foreach(backend IN LISTS CHAIN_TEST_BACKENDS)
      add_chain_test(${test} ${backend} MODE initializer SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/${test}
                     INPUTS param ${TEMPLATE50_INPUTS}
                     COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/runtest ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${test} ${Mathematica_KERNEL_EXECUTABLE})
      list(APPEND profile_tests ${test}_${backend})
    endforeach()
  endif()
endforeach()

if(GPROF)
  gprof("${profile_tests}" ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
endif()
