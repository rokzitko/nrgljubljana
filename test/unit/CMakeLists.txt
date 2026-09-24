file(GLOB_RECURSE all_tests RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

set(TEST_REFERENCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/reference)
file(MAKE_DIRECTORY ${TEST_REFERENCE_DIR})
foreach(reference_file polynomials.tsv moments.tsv cauchy_transforms.tsv kk_nonlinear.tsv)
  configure_file(${PROJECT_SOURCE_DIR}/test/reference/generated/${reference_file}
                 ${TEST_REFERENCE_DIR}/${reference_file} COPYONLY)
endforeach()
configure_file(include/cmake_configure.hpp.in include/cmake_configure.hpp)

option(GCOV "gcov" OFF)

if(GCOV)
  find_package(Perl REQUIRED)
  find_program(LCOV_EXECUTABLE lcov REQUIRED)
  find_program(GENHTML_EXECUTABLE genhtml REQUIRED)
  message(STATUS "*** gcov instrumentalization enabled")
  set(gtest_test_name "gtest_")
  add_custom_target(gcov
                    COMMAND ${CMAKE_COMMAND} -E rm -rf ${CMAKE_BINARY_DIR}/test/unit
                    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR}
                    COMMAND ${CMAKE_CTEST_COMMAND} -j2 -R gtest --output-on-failure
                    COMMAND ${CMAKE_COMMAND} -E rm -rf ${CMAKE_BINARY_DIR}/gcoverage
                    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/gcoverage
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR} )

  add_custom_command(TARGET gcov POST_BUILD
                     COMMAND ${CMAKE_COMMAND} -E env
                       "LCOV_EXECUTABLE=${LCOV_EXECUTABLE}"
                       "GENHTML_EXECUTABLE=${GENHTML_EXECUTABLE}"
                       ${CMAKE_CURRENT_SOURCE_DIR}/gcov/gcov
                      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/gcoverage )
endif()

foreach(test ${all_tests})
  get_filename_component(test_name ${test} NAME_WE)
  get_filename_component(test_dir ${test} DIRECTORY)
  string(REPLACE "/" "_" test_resource_lock "unit_${test_dir}")
  add_executable(${test_name} ${test})
  target_compile_options(${test_name} PRIVATE $<$<CXX_COMPILER_ID:GNU>:-fconcepts>
                                               $<$<BOOL:${GCOV}>:--coverage>)
  nrgljubljana_apply_strict_floating_point(${test_name})
  set_property(TARGET ${test_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})
  set_property(TARGET ${test_name} PROPERTY BUILD_RPATH "${PROJECT_BINARY_DIR}/c++")
  target_link_libraries(${test_name} nrgljubljana_c gtest_main project_warnings gsl Threads::Threads)
  target_link_options(${test_name} PRIVATE  $<$<BOOL:${GCOV}>:--coverage>)
  target_include_directories(${test_name} PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${PROJECT_BINARY_DIR}/tools/generated)
  if(test_name STREQUAL "scalar_tridiag")
    add_test(NAME ${gtest_test_name}scalar_tridiag_params COMMAND scalar_tridiag --gtest_filter=ScalarTridiagParams.*
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})
    set_tests_properties(${gtest_test_name}scalar_tridiag_params PROPERTIES
      ENVIRONMENT "${TEST_BUILD_ENV}" RESOURCE_LOCK "${test_resource_lock}")
    foreach(backend ${CHAIN_TEST_BACKENDS})
      add_test(NAME ${gtest_test_name}scalar_tridiag_${backend}
        COMMAND scalar_tridiag --gtest_filter=ScalarTridiag/*.*_${backend}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})
      set_tests_properties(${gtest_test_name}scalar_tridiag_${backend} PROPERTIES
        ENVIRONMENT "${TEST_BUILD_ENV}" RESOURCE_LOCK "${test_resource_lock}" LABELS "chain-${backend}")
    endforeach()
    if(TEST_CHAIN_CROSSCHECK AND TEST_CHAIN_LEGACY AND TEST_CHAIN_RKPW)
      add_test(NAME ${gtest_test_name}scalar_tridiag_crosscheck
        COMMAND scalar_tridiag --gtest_filter=ScalarTridiag/*.*_crosscheck
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})
      set_tests_properties(${gtest_test_name}scalar_tridiag_crosscheck PROPERTIES
        ENVIRONMENT "${TEST_BUILD_ENV}" RESOURCE_LOCK "${test_resource_lock}"
        LABELS "chain-crosscheck;chain-uses-legacy;chain-uses-rkpw")
    endif()
  elseif(test_name STREQUAL "star_to_chain")
    if(TEST_CHAIN_RKPW)
      add_test(NAME ${gtest_test_name}star_to_chain_rkpw COMMAND star_to_chain
               WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})
      set_tests_properties(${gtest_test_name}star_to_chain_rkpw PROPERTIES
        ENVIRONMENT "${TEST_BUILD_ENV}" RESOURCE_LOCK "${test_resource_lock}" LABELS chain-rkpw)
    endif()
  else()
    add_test(NAME ${gtest_test_name}${test_name} COMMAND ${test_name} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})
    set_tests_properties(${gtest_test_name}${test_name} PROPERTIES ENVIRONMENT "${TEST_BUILD_ENV}" RESOURCE_LOCK "${test_resource_lock}")
  endif()

  file(GLOB_RECURSE txtfiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${test_dir}/*.txt)
  foreach(file ${txtfiles})
    configure_file(${file} ${file} COPYONLY)
  endforeach()

endforeach()
