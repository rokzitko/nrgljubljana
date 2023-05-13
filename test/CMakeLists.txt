if(GPROF)
  message(STATUS "*** gprof instrumentalization enabled")
  set(GPROF_DEPEND)
  add_custom_target(gprof_start
                    COMMAND echo "gprof start"
                    COMMAND make || exit 1
                    COMMAND rm -rf gmon.sum
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR} )

  macro(gprof simple_tests current_source_dir current_binary_dir)
    get_filename_component(dir ${current_source_dir} NAME_WE)
    add_custom_target(gprof_${dir}
                      COMMAND echo "gprof_${dir} start"
                      COMMAND ctest
                      WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
    set(GPROF_DEPEND ${GPROF_DEPEND} "gprof_${dir}" PARENT_SCOPE)
    add_custom_command(TARGET gprof_${dir}
                      COMMAND for test in ${simple_tests}\; \\
                      COMMAND do \\
                      COMMAND if [ -f \"${CMAKE_BINARY_DIR}/gmon.sum\" ] \; \\
                      COMMAND then \\
                      COMMAND   gprof ${CMAKE_BINARY_DIR}/c++/nrg -s -z -c ${CMAKE_BINARY_DIR}/gmon.sum ${current_binary_dir}/$$test/gmon.out\; \\
                      COMMAND else \\
                      COMMAND   cp ${current_binary_dir}/$$test/gmon.out  ${CMAKE_BINARY_DIR}/gmon.sum\; \\
                      COMMAND fi \\
                      COMMAND done
                        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                        )
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

if(Mathematica_FOUND AND SYM_ALL)
  add_subdirectory(nrginit)
  add_subdirectory(nrginit_spsu2)
  add_subdirectory(nrginit+nrgrun)
  if (TEST_LONG)
    add_subdirectory(nrginit+nrgrun_long)
  endif()
endif()

if(GPROF)
  add_custom_target(gprof DEPENDS ${GPROF_DEPEND}
                    COMMAND echo "gprof end"
                    COMMAND gprof -z -c ${CMAKE_BINARY_DIR}/c++/nrg ${CMAKE_BINARY_DIR}/gmon.sum > gprof.result
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )
endif()
