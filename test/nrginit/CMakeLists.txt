find_package(Perl REQUIRED)

# Regression tests

file(STRINGS simple_tests simple_tests)

foreach(test ${simple_tests})
  file(GLOB inputs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/${test} ${test}/*.m)
  if(test STREQUAL "nrginit2_generate_template" OR test STREQUAL "nrginit47_GENERATE_TEMPLATE")
    # These generate symbolic templates, not numerical Wilson chains.
    configure_file(${test}/param ${test}/param COPYONLY)
    add_test(NAME ${test} COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/runtest ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${test} ${Mathematica_KERNEL_EXECUTABLE}
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test})
    set_tests_properties(${test} PROPERTIES ENVIRONMENT "${TEST_BUILD_ENV}")
  else()
    foreach(backend IN LISTS CHAIN_TEST_BACKENDS)
      add_chain_test(${test} ${backend} MODE initializer SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/${test}
                     INPUTS param ${inputs}
                     COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/runtest ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${test} ${Mathematica_KERNEL_EXECUTABLE})
    endforeach()
  endif()
endforeach()

add_test(NAME nrginit_bool_parser
         COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/test_bool_parser ${PROJECT_SOURCE_DIR} ${Mathematica_KERNEL_EXECUTABLE})

add_test(NAME nrginit_floquet_primitives
         COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/test_floquet_primitives ${PROJECT_SOURCE_DIR} ${Mathematica_KERNEL_EXECUTABLE})

add_test(NAME nrginit_chain_defaults
         COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/test_rkpw ${PROJECT_SOURCE_DIR} ${Mathematica_KERNEL_EXECUTABLE} chain_defaults)

add_test(NAME nrginit_runner_contracts
          COMMAND ${PERL_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_runner_contracts.pl ${PROJECT_SOURCE_DIR})

foreach(backend IN LISTS CHAIN_TEST_BACKENDS)
  add_test(NAME nrginit_mapping_${backend}
           COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/test_rkpw ${PROJECT_SOURCE_DIR} ${Mathematica_KERNEL_EXECUTABLE} chain_mapping)
  set_tests_properties(nrginit_mapping_${backend} PROPERTIES
    ENVIRONMENT "${TEST_BUILD_ENV};CHAIN_TEST_BACKEND=${backend}" LABELS chain-${backend})
  foreach(suite dmft_empty hardgap_mapping)
    add_test(NAME nrginit_${suite}_${backend}
             COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/test_rkpw ${PROJECT_SOURCE_DIR} ${Mathematica_KERNEL_EXECUTABLE} ${suite})
    set_tests_properties(nrginit_${suite}_${backend} PROPERTIES
      ENVIRONMENT "${TEST_BUILD_ENV};CHAIN_TEST_BACKEND=${backend}" LABELS chain-${backend})
  endforeach()
  add_test(NAME nrginit_serialization_${backend}
           COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/test_chain_serialization ${PROJECT_SOURCE_DIR} ${Mathematica_KERNEL_EXECUTABLE})
  set_tests_properties(nrginit_serialization_${backend} PROPERTIES
    ENVIRONMENT "${TEST_BUILD_ENV};CHAIN_TEST_BACKEND=${backend}" LABELS chain-${backend})
endforeach()

if(TEST_CHAIN_LEGACY)
  add_test(NAME nrginit_chain_legacy
           COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/test_rkpw ${PROJECT_SOURCE_DIR} ${Mathematica_KERNEL_EXECUTABLE} chain_legacy)
  set_tests_properties(nrginit_chain_legacy PROPERTIES LABELS chain-legacy)
endif()

if(TEST_CHAIN_RKPW)
  add_test(NAME nrginit_rkpw
           COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/test_rkpw ${PROJECT_SOURCE_DIR} ${Mathematica_KERNEL_EXECUTABLE})
  add_test(NAME nrginit_pipeline_rkpw
           COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/test_rkpw_pipeline ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR}
                        ${Mathematica_KERNEL_EXECUTABLE})
  set_tests_properties(nrginit_rkpw nrginit_pipeline_rkpw PROPERTIES LABELS chain-rkpw ENVIRONMENT "${TEST_BUILD_ENV}")
endif()

option(TEST_CHAIN_CROSSCHECK "Compare initializer legacy and RKPW implementations" OFF)
if(TEST_CHAIN_CROSSCHECK AND TEST_CHAIN_LEGACY AND TEST_CHAIN_RKPW)
  add_test(NAME nrginit_chain_crosscheck
           COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/test_rkpw ${PROJECT_SOURCE_DIR} ${Mathematica_KERNEL_EXECUTABLE} chain_crosscheck)
  set_tests_properties(nrginit_chain_crosscheck PROPERTIES LABELS "chain-crosscheck;chain-uses-legacy;chain-uses-rkpw")
endif()
