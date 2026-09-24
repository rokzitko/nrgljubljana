# Regression tests

file(STRINGS simple_tests simple_tests)
set(registered_tests)

set(chain_runtime_tests
  test40_tri test41_tri_mMAX40 test42_model_1orb test43_model_2orb
  test44_model_2orb_U test45_model_2orb_JU test46_adapt
  legacy_discretization_dmft_asym_qs legacy_discretization_ferromagnetic1
  legacy_spectral_gtls-qsz-T1e-10)

foreach(test ${simple_tests})
  if(test IN_LIST chain_runtime_tests)
    foreach(backend ${CHAIN_TEST_BACKENDS})
      add_chain_test(${test} ${backend} MODE runtime SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/${test}
        INPUTS param data PARAM_FILES param
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/runtest ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${test}
        LABELS chain-fixed-seed)
      list(APPEND registered_tests ${test}_${backend})
    endforeach()
    continue()
  endif()
  configure_file(${test}/param ${test}/param COPYONLY)
  if(test STREQUAL "floquet_spu1_runtime")
    configure_file(${PROJECT_SOURCE_DIR}/test/nrginit+nrgrun/floquet_spu1/data ${test}/data COPYONLY)
  else()
    configure_file(${test}/data ${test}/data COPYONLY)
  endif()
  add_test(NAME ${test} COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/runtest ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${test}
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test} )
  set_tests_properties(${test} PROPERTIES ENVIRONMENT "${TEST_BUILD_ENV}")
  list(APPEND registered_tests ${test})
endforeach()

# Regenerate the inputs independently of the fixed-seed runtime tests above.
if(Mathematica_FOUND)
  set(chain_prepare_tests
    test46_adapt test47_nrgchain test48_adapt_tri=old test49_template1
    test50_ref1 test51_template2 test52_ref2 test53_template3 test54_ref3)
  foreach(test ${chain_prepare_tests})
    set(mode initializer)
    set(inputs param Delta.dat prep)
    if(test MATCHES "^test(47|49|51|53)_")
      set(mode tool)
    endif()
    if(test MATCHES "^test(49|51|53)_")
      list(APPEND inputs prepare.pl)
    endif()
    if(test STREQUAL "test53_template3")
      # Imported only while generating symbolic templates, then discarded.
      # prepare.pl and the solver use freshly reconstructed coefficients.
      list(APPEND inputs xi.dat zeta.dat)
    endif()
    foreach(backend ${CHAIN_TEST_BACKENDS})
      add_chain_test(${test}_prepare ${backend} MODE ${mode} SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/${test}
        INPUTS ${inputs} PARAM_FILES param
        COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/runprepare ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR}
          ${CMAKE_CURRENT_SOURCE_DIR}/${test} ${Mathematica_KERNEL_EXECUTABLE}
        LABELS chain-preparation)
      list(APPEND registered_tests ${test}_prepare_${backend})
    endforeach()
  endforeach()
endif()

configure_file(test_resume/param test_resume/param COPYONLY)
configure_file(test0_minimal_Nmax3/data test_resume/data COPYONLY)
add_test(NAME test_resume COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/resume_test $<TARGET_FILE:nrg> 2
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_resume)
set_tests_properties(test_resume PROPERTIES ENVIRONMENT "${TEST_BUILD_ENV}")

configure_file(test_resume_fdm/param test_resume_fdm/param COPYONLY)
configure_file(test_fdm_only/data test_resume_fdm/data COPYONLY)
add_test(NAME test_resume_fdm COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/resume_test $<TARGET_FILE:nrg> 3 density
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_resume_fdm)
set_tests_properties(test_resume_fdm PROPERTIES ENVIRONMENT "${TEST_BUILD_ENV}")

if(GPROF)
  gprof("${registered_tests}" ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
endif()

if(NRG_MPI)
  file(STRINGS mpi_tests mpi_tests)

  macro(add_mpi_tests NumProc)
    foreach(test ${mpi_tests})
      configure_file(${test}/param ${test}_np${NumProc}/param COPYONLY)
      configure_file(${test}/data ${test}_np${NumProc}/data COPYONLY)
      add_test(NAME ${test}_np${NumProc}
              COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/runmpitest
                ${PROJECT_SOURCE_DIR}
                ${CMAKE_CURRENT_SOURCE_DIR}/${test}
                ${MPIEXEC_EXECUTABLE}
                ${MPIEXEC_PREFLAGS}
                ${MPIEXEC_NUMPROC_FLAG}
                ${NumProc}
                $<TARGET_FILE:nrg>
                ${MPIEXEC_POSTFLAGS}
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test}_np${NumProc})
      set_tests_properties(${test}_np${NumProc} PROPERTIES ENVIRONMENT "${TEST_BUILD_ENV}" PROCESSORS ${NumProc})
    endforeach()
  endmacro()

  add_mpi_tests(1)
  add_mpi_tests(2)
  add_mpi_tests(3)
  add_mpi_tests(4)
endif()
