# Regression tests

file(STRINGS simple_tests simple_tests)

foreach(test ${simple_tests})
  foreach(backend IN LISTS CHAIN_TEST_BACKENDS)
    add_chain_test(${test} ${backend} MODE initializer SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/${test} INPUTS param
                   COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/runtest ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${test} ${Mathematica_KERNEL_EXECUTABLE})
  endforeach()
endforeach()

add_test(NAME zero_u1_matsubara_contracts
         COMMAND ${PERL_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/zero_u1/test_validate_matsubara.pl)
