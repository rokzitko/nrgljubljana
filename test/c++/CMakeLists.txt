# Regression tests

file(STRINGS simple_tests simple_tests)

foreach(test ${simple_tests})
  configure_file(${test}/param ${test}/param COPYONLY)
  configure_file(${test}/data ${test}/data COPYONLY)
  add_test(NAME ${test} COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/runtest ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${test}
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test} )
endforeach()

if(GPROF)
  gprof("${simple_tests}" ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
endif()

if (NRG_MPI)
 file(STRINGS mpi_tests mpi_tests)

 macro(add_mpi_tests NumProc)
   foreach(test ${mpi_tests})
     configure_file(${test}/param ${test}_np${NumProc}/param COPYONLY)
     configure_file(${test}/data ${test}_np${NumProc}/data COPYONLY)
     add_test(NAME ${test}_np${NumProc} COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/runmpitest ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${test} ${NumProc}
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test}_np${NumProc})
   endforeach()
 endmacro()

 add_mpi_tests(1)
 add_mpi_tests(2)
 add_mpi_tests(3)
 add_mpi_tests(4)
endif()