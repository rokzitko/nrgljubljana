find_package(Perl REQUIRED)

# Regression tests

file(STRINGS simple_tests simple_tests)

foreach(test ${simple_tests})
  configure_file(${test}/run ${test}/run COPYONLY) # also ensures the test dir is created
  add_test(NAME ${test} COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/${test}/run ${PROJECT_BINARY_DIR}/tools ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${test}
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test})
endforeach()

configure_file(hilbert2_DOS/DOS.dat hilbert2_DOS/DOS.dat COPYONLY)
configure_file(hilbert3_DOS_fine/DOS.dat hilbert3_DOS_fine/DOS.dat COPYONLY)
configure_file(hilbert7_file/points.dat hilbert7_file/points.dat COPYONLY)
configure_file(hilbert8_G/DOS.dat hilbert8_G/DOS.dat COPYONLY)
configure_file(hilbert8_G/resigma.dat hilbert8_G/resigma.dat COPYONLY)
configure_file(hilbert8_G/imsigma.dat hilbert8_G/imsigma.dat COPYONLY)
configure_file(hilbert9_rescale/points.dat hilbert9_rescale/points.dat COPYONLY)
configure_file(hilbert10_bandwidth/points.dat hilbert10_bandwidth/points.dat COPYONLY)
