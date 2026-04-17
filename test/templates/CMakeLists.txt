# Regression tests

file(STRINGS simple_tests simple_tests)

foreach(test ${simple_tests})
  configure_file(${test}/param ${test}/param COPYONLY)
  add_test(NAME ${test} COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/runtest ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${test} ${Mathematica_KERNEL_EXECUTABLE}
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test})
endforeach()

configure_file(template49_ref/data template49_ref/data COPYONLY)

file(STRINGS template50_inputs TEMPLATE50_INPUTS)
foreach(file IN LISTS TEMPLATE50_INPUTS)
  configure_file(template50_template/${file} template50_template/${file} COPYONLY)
endforeach()

if(GPROF)
  gprof("${simple_tests}" ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
endif()
