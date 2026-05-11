set(TEST_BUILD_ENV
  "LD_LIBRARY_PATH=${PROJECT_BINARY_DIR}/c++:$ENV{LD_LIBRARY_PATH}"
  "DYLD_FALLBACK_LIBRARY_PATH=${PROJECT_BINARY_DIR}/c++:$ENV{DYLD_FALLBACK_LIBRARY_PATH}")

if(ASAN)
  set(TEST_ASAN_OPTIONS "$ENV{ASAN_OPTIONS}")
  if(TEST_ASAN_OPTIONS)
    string(APPEND TEST_ASAN_OPTIONS ":")
  endif()
  # MPI/PMIx/libfabric stacks commonly retain shutdown allocations that are
  # outside nrgljubljana; keep ASan checks but disable LeakSanitizer for CTest.
  string(APPEND TEST_ASAN_OPTIONS "detect_leaks=0")
  string(REPLACE ";" "\\;" TEST_ASAN_OPTIONS "${TEST_ASAN_OPTIONS}")
  list(APPEND TEST_BUILD_ENV "ASAN_OPTIONS=${TEST_ASAN_OPTIONS}")
endif()

if(GPROF)
  message(STATUS "*** gprof instrumentalization enabled")
  find_program(GPROF_EXECUTABLE gprof REQUIRED)
  set(NRGLJUBLJANA_GPROF_SUM_SCRIPT ${CMAKE_CURRENT_LIST_DIR}/gprof_sum.cmake)
  set(NRGLJUBLJANA_GPROF_REPORT_SCRIPT ${CMAKE_CURRENT_LIST_DIR}/gprof_report.cmake)
  set(GPROF_DEPEND)
  add_custom_target(gprof_start
                    COMMAND ${CMAKE_COMMAND} -E echo "gprof start"
                    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR}
                    COMMAND ${CMAKE_COMMAND} -E rm -f ${CMAKE_BINARY_DIR}/gmon.sum
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR} )

  macro(gprof simple_tests current_source_dir current_binary_dir)
    get_filename_component(dir ${current_source_dir} NAME_WE)
    add_custom_target(gprof_${dir}
                      COMMAND ${CMAKE_COMMAND} -E echo "gprof_${dir} start"
                      COMMAND ${CMAKE_CTEST_COMMAND}
                      WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
    set(GPROF_DEPEND ${GPROF_DEPEND} "gprof_${dir}" PARENT_SCOPE)
    add_custom_command(TARGET gprof_${dir} POST_BUILD
                       COMMAND ${CMAKE_COMMAND}
                         "-DSIMPLE_TESTS=${simple_tests}"
                        "-DCURRENT_BINARY_DIR=${current_binary_dir}"
                        "-DBINARY_DIR=${CMAKE_BINARY_DIR}"
                        "-DNRG_EXECUTABLE=$<TARGET_FILE:nrg>"
                        "-DGMON_SUM=${CMAKE_BINARY_DIR}/gmon.sum"
                        "-DGPROF_EXECUTABLE=${GPROF_EXECUTABLE}"
                        -P ${NRGLJUBLJANA_GPROF_SUM_SCRIPT}
                      WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
  endmacro()
endif()

add_subdirectory(unit)

add_subdirectory(c++)
add_subdirectory(c++sym_basic)

if(SYM_MORE)
  add_subdirectory(c++sym_more)
  add_subdirectory(c++sym_spsu2)
endif()

if(SYM_ALL)
  add_subdirectory(c++sym_all)
  add_subdirectory(complex)
endif()

if(TEST_LONG AND SYM_ALL)
  add_subdirectory(test_long)
endif()

add_subdirectory(tools)
add_subdirectory(nrgspawn)

if(Mathematica_FOUND AND SYM_ALL)
  add_subdirectory(nrginit)
  add_subdirectory(nrginit_spsu2)
  add_subdirectory(nrginit+nrgrun)
  add_subdirectory(templates)
  if (TEST_LONG)
    add_subdirectory(nrginit+nrgrun_long)
  endif()
  if (SYM_ALL)
      add_subdirectory(models)
      if (TEST_LONG)
        add_subdirectory(models_long)
      endif()
  endif()
endif()

if(GPROF)
  add_custom_target(gprof DEPENDS ${GPROF_DEPEND}
                    COMMAND ${CMAKE_COMMAND} -E echo "gprof end"
                    COMMAND ${CMAKE_COMMAND}
                      "-DNRG_EXECUTABLE=$<TARGET_FILE:nrg>"
                      "-DGMON_SUM=${CMAKE_BINARY_DIR}/gmon.sum"
                      "-DOUTPUT_FILE=${CMAKE_BINARY_DIR}/gprof.result"
                      "-DGPROF_EXECUTABLE=${GPROF_EXECUTABLE}"
                      -P ${NRGLJUBLJANA_GPROF_REPORT_SCRIPT}
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )
endif()
