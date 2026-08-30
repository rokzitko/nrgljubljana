find_package(Perl REQUIRED)

# Regression tests

file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/simple_tests" simple_tests)

foreach(test IN LISTS simple_tests)
  set(test_source_dir "${CMAKE_CURRENT_SOURCE_DIR}/${test}")
  set(test_binary_dir "${CMAKE_CURRENT_BINARY_DIR}/${test}")
  file(MAKE_DIRECTORY "${test_binary_dir}")
  file(GLOB test_fixtures CONFIGURE_DEPENDS RELATIVE "${test_source_dir}" "${test_source_dir}/*.dat")
  foreach(fixture IN LISTS test_fixtures)
    configure_file("${test_source_dir}/${fixture}" "${test_binary_dir}/${fixture}" COPYONLY)
  endforeach()
  add_test(NAME ${test}
           COMMAND sh "${test_source_dir}/run" "${PROJECT_BINARY_DIR}/tools" "${PROJECT_SOURCE_DIR}"
                   "${test_source_dir}" "${PERL_EXECUTABLE}"
           WORKING_DIRECTORY "${test_binary_dir}")
endforeach()
