find_package(Perl REQUIRED)

# Regression tests

file(STRINGS simple_tests simple_tests)

foreach(test ${simple_tests})
  configure_file(${test}/param ${test}/param COPYONLY)
  add_test(NAME ${test} COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/runtest ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${test} ${Mathematica_KERNEL_EXECUTABLE}
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test})
endforeach()

configure_file(nrginit10_calcopq_d2/SIAM.m nrginit10_calcopq_d2/SIAM.m COPYONLY)
configure_file(nrginit11_calcopq_d3/SIAM.m nrginit11_calcopq_d3/SIAM.m COPYONLY)
configure_file(nrginit12_calcopq_d4/DQD.m nrginit12_calcopq_d4/DQD.m COPYONLY)
