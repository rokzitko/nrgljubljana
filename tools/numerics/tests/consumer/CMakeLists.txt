cmake_minimum_required(VERSION 3.25)
project(numerics_consumer LANGUAGES CXX)

if(NRGLJUBLJANA_NUMERICS_SOURCE_DIR)
  add_subdirectory("${NRGLJUBLJANA_NUMERICS_SOURCE_DIR}" numerics)
else()
  find_package(nrgljubljana_numerics 1 CONFIG REQUIRED)
endif()

option(NRGLJUBLJANA_NUMERICS_CHECK_SOLVER_PACKAGE "Also import the installed solver package" OFF)
if(NRGLJUBLJANA_NUMERICS_CHECK_SOLVER_PACKAGE)
  find_package(nrgljubljana CONFIG REQUIRED)
endif()

add_executable(consumer main.cpp)
target_link_libraries(consumer PRIVATE nrgljubljana::numerics)
enable_testing()
add_test(NAME numerics.consumer COMMAND consumer)
