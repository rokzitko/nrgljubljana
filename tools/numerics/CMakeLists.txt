cmake_minimum_required(VERSION 3.25...3.31)
project(nrgljubljana_numerics VERSION 1.0.0 LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

find_package(GSL 2.1 REQUIRED)
find_package(Threads REQUIRED)
if(NOT TARGET Boost::headers)
  find_package(Boost REQUIRED)
endif()

add_library(nrgljubljana_numerics src/numerics.cpp)
add_library(nrgljubljana::numerics ALIAS nrgljubljana_numerics)
target_compile_features(nrgljubljana_numerics PUBLIC cxx_std_20)
target_include_directories(nrgljubljana_numerics PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(nrgljubljana_numerics PRIVATE
  GSL::gsl GSL::gslcblas Threads::Threads $<BUILD_INTERFACE:Boost::headers>)
# Finite checks and the compensated/exact fallback kernels require strict FP.
target_compile_options(nrgljubljana_numerics PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-fno-fast-math;-ffp-contract=off>
  $<$<CXX_COMPILER_ID:MSVC>:/fp:strict>)
set_target_properties(nrgljubljana_numerics PROPERTIES
  EXPORT_NAME numerics
  POSITION_INDEPENDENT_CODE ON
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR})

if(PROJECT_IS_TOP_LEVEL)
  include(CTest)
  set(numerics_test_default ${BUILD_TESTING})
elseif(NRGLJUBLJANA_BUILD_TESTS_ENABLED)
  set(numerics_test_default ON)
else()
  set(numerics_test_default OFF)
endif()
option(NRGLJUBLJANA_NUMERICS_BUILD_TESTS "Build the focused numerics tests" ${numerics_test_default})
if(NRGLJUBLJANA_NUMERICS_BUILD_TESTS)
  enable_testing()
  add_executable(nrgljubljana_numerics_tests tests/numerics.cpp)
  target_link_libraries(nrgljubljana_numerics_tests PRIVATE nrgljubljana::numerics GSL::gsl Threads::Threads)
  foreach(test IN ITEMS integration resample kk lattice validation gsl_handler)
    add_test(NAME numerics.${test} COMMAND nrgljubljana_numerics_tests ${test})
  endforeach()
endif()

set(numerics_package_dir "${CMAKE_INSTALL_LIBDIR}/cmake/nrgljubljana_numerics")
configure_package_config_file(cmake/nrgljubljana_numericsConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/nrgljubljana_numericsConfig.cmake"
  INSTALL_DESTINATION "${numerics_package_dir}")
write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/nrgljubljana_numericsConfigVersion.cmake"
  VERSION ${PROJECT_VERSION} COMPATIBILITY SameMajorVersion)
install(TARGETS nrgljubljana_numerics EXPORT nrgljubljana_numericsTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(FILES include/nrgljubljana/numerics.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nrgljubljana)
install(EXPORT nrgljubljana_numericsTargets NAMESPACE nrgljubljana::
  DESTINATION "${numerics_package_dir}")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/nrgljubljana_numericsConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/nrgljubljana_numericsConfigVersion.cmake"
  DESTINATION "${numerics_package_dir}")
