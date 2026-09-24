file(STRINGS simple_tests simple_tests)

foreach(test ${simple_tests})
  set(inputs param)
  foreach(file template data.in ham_0)
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${test}/${file})
      list(APPEND inputs ${file})
    endif()
  endforeach()
  foreach(backend IN LISTS CHAIN_TEST_BACKENDS)
    add_chain_test(${test} ${backend} MODE tool SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/${test}
      INPUTS ${inputs}
      COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/${test}/run ${PROJECT_BINARY_DIR}/tools ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${test})
  endforeach()
endforeach()
