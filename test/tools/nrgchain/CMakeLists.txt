find_package(Perl REQUIRED)

# Regression tests

file(STRINGS simple_tests simple_tests)

foreach(test ${simple_tests})
  set(param_files param)
  if(test STREQUAL "nrgchain19_steffen")
    list(APPEND param_files load.param linear.param conflict.param no-tridiag.param)
  endif()
  foreach(backend IN LISTS CHAIN_TEST_BACKENDS)
    add_chain_test(${test} ${backend} MODE tool SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/${test}
      INPUTS ${param_files} Delta.dat PARAM_FILES ${param_files}
      COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/${test}/run ${PROJECT_BINARY_DIR}/tools ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${test})
  endforeach()
endforeach()
