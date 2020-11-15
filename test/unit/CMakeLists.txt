add_subdirectory(gtest EXCLUDE_FROM_ALL)

file(GLOB_RECURSE all_tests RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

configure_file(nextline.txt nextline.txt COPYONLY)
configure_file(block.txt block.txt COPYONLY)

foreach(test ${all_tests})
  get_filename_component(test_name ${test} NAME_WE)
  get_filename_component(test_dir ${test} DIRECTORY)
  add_executable(${test_name} ${test})
  target_compile_options(${test_name} PRIVATE -std=c++17 $<$<CXX_COMPILER_ID:GNU>:-fconcepts>)
  set_property(TARGET ${test_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})
  target_link_libraries(${test_name} nrgljubljana_c gtest_main project_warnings)
  add_test(NAME ${test_name} COMMAND ${test_name} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})
endforeach()
