find_package(Perl REQUIRED)

# Regression tests

file(STRINGS simple_tests simple_tests)

foreach(test ${simple_tests})
  configure_file(${test}/run ${test}/run COPYONLY) # also ensures the test dir is created
  add_test(NAME ${test} COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/${test}/run ${PROJECT_BINARY_DIR}/tools ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${test}
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test})
endforeach()

configure_file(hilb2_DOS/DOS.dat hilb2_DOS/DOS.dat COPYONLY)
configure_file(hilb3_DOS_fine/DOS.dat hilb3_DOS_fine/DOS.dat COPYONLY)
configure_file(hilb7_file/points.dat hilb7_file/points.dat COPYONLY)
configure_file(hilb8_G/DOS.dat hilb8_G/DOS.dat COPYONLY)
configure_file(hilb8_G/resigma.dat hilb8_G/resigma.dat COPYONLY)
configure_file(hilb8_G/imsigma.dat hilb8_G/imsigma.dat COPYONLY)
configure_file(hilb9_rescale/points.dat hilb9_rescale/points.dat COPYONLY)
configure_file(hilb10_bandwidth/points.dat hilb10_bandwidth/points.dat COPYONLY)
