# Regression tests

file(STRINGS simple_tests simple_tests)

foreach(test ${simple_tests})
  configure_file(${test}/param ${test}/param COPYONLY)
  configure_file(${test}/model.m ${test}/model.m COPYONLY)
  configure_file(${test}/modeloperators.m ${test}/modeloperators.m COPYONLY)
  add_test(NAME ${test} COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/runtest ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${test} ${Mathematica_KERNEL_EXECUTABLE}
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test})
endforeach()
