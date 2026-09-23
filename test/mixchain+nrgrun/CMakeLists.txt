# End-to-end tests: a Wilson chain from mixchain, staged into a SIAM template, run by nrg.
#
# Two kinds live here. Most compare against a stored ref/ and share the 'runtest' script. A few instead build the
# same physics a second way and compare the two runs, with no stored numbers; those carry their own 'run' script and
# are listed in the same files. Which tests can run depends on the symmetry types the build has: QS is always
# compiled, SPSU2 needs SYM_MORE and U1 needs SYM_ALL (c++/CMakeLists.txt).

find_package(Perl REQUIRED)

set(mixchain_nrgrun_tests)
file(STRINGS simple_tests simple_tests)
list(APPEND mixchain_nrgrun_tests ${simple_tests})

if(SYM_MORE)
  file(STRINGS more_tests more_tests)
  list(APPEND mixchain_nrgrun_tests ${more_tests})
endif()

if(SYM_ALL)
  file(STRINGS all_tests all_tests)
  list(APPEND mixchain_nrgrun_tests ${all_tests})
endif()

foreach(test ${mixchain_nrgrun_tests})
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${test}/run)
    set(command ${CMAKE_CURRENT_SOURCE_DIR}/${test}/run)
  else()
    set(command ${CMAKE_CURRENT_SOURCE_DIR}/runtest)
  endif()
  add_test(NAME ${test}
           COMMAND ${command} ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${test}
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test})
  set_tests_properties(${test} PROPERTIES ENVIRONMENT "${TEST_BUILD_ENV}")
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test})
endforeach()

# The Mathematica route as a second opinion. It needs a kernel, which is passed as a fourth argument, and SPSU2.
if(Mathematica_FOUND AND SYM_MORE)
  file(STRINGS mathematica_tests mathematica_tests)
  foreach(test ${mathematica_tests})
    add_test(NAME ${test}
             COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/${test}/run ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${test} ${Mathematica_KERNEL_EXECUTABLE}
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test})
    set_tests_properties(${test} PROPERTIES ENVIRONMENT "${TEST_BUILD_ENV}")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test})
  endforeach()
  list(APPEND mixchain_nrgrun_tests ${mathematica_tests})
endif()

if(GPROF)
  gprof("${mixchain_nrgrun_tests}" ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
endif()
