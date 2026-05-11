# Regression tests

file(STRINGS all_tests all_tests)

foreach(test ${all_tests})
  configure_file(${test}/param ${test}/param COPYONLY)
  add_test(NAME ${test} COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/runtest ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${test}
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test})
endforeach()

file(STRINGS inst1_inputs INST1_INPUTS)
foreach(file IN LISTS INST1_INPUTS)
  configure_file(inst1_Anderson/${file} inst1_Anderson/${file} COPYONLY)
endforeach()

if(GPROF)
  gprof("${all_tests}" ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
endif()
