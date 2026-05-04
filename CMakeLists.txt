# ##############################################################################
#
# NRG Ljubljana
#
# Copyright (C) Rok Zitko, rok.zitko@ijs.si
#
# NRG Ljubljana is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# NRG Ljubljana is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# nrgljubljana (in the file COPYING.txt in this directory). If not, see
# <http://www.gnu.org/licenses/>.
#
# ##############################################################################

cmake_minimum_required(VERSION 3.25...3.31 FATAL_ERROR)
if(POLICY CMP0173) # CMakeFindFrameworks Availability
  cmake_policy(SET CMP0173 OLD)
endif()

# ############
# Define Project
project(NRGLJUBLJANA VERSION 2026.04 LANGUAGES C CXX)
get_directory_property(IS_SUBPROJECT PARENT_DIRECTORY)

# Get the git hash & print status
execute_process(COMMAND git rev-parse HEAD WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_VARIABLE PROJECT_GIT_HASH OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "${PROJECT_NAME} version : ${PROJECT_VERSION}")
message(STATUS "${PROJECT_NAME} Git hash: ${PROJECT_GIT_HASH}")

# Assert that Install directory is given and valid.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT OR (NOT IS_ABSOLUTE ${CMAKE_INSTALL_PREFIX}))
  message(FATAL_ERROR "No install prefix given (or invalid)")
endif()
message(STATUS "-------- CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX} --------")

# ############
# Options

# Make additional Find Modules available
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/share/cmake/Modules)

# FindMathematica (optional, default OFF only on aarch64; overridable from CLI)
string(TOLOWER "${CMAKE_SYSTEM_PROCESSOR}" NRGLJUBLJANA_SYSTEM_PROCESSOR_LOWER)
set(NRGLJUBLJANA_ENABLE_MATHEMATICA_DEFAULT ON)
if(NRGLJUBLJANA_SYSTEM_PROCESSOR_LOWER STREQUAL "aarch64")
  set(NRGLJUBLJANA_ENABLE_MATHEMATICA_DEFAULT OFF)
endif()
option(NRGLJUBLJANA_ENABLE_MATHEMATICA "Enable Mathematica detection via FindMathematica" ${NRGLJUBLJANA_ENABLE_MATHEMATICA_DEFAULT})
if(NRGLJUBLJANA_ENABLE_MATHEMATICA)
  list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake/Mathematica)
else()
  message(STATUS "Mathematica detection disabled (NRGLJUBLJANA_ENABLE_MATHEMATICA=OFF)")
endif()
option(NRGLJUBLJANA_USE_SYSTEM_DEPS "Use dependencies supplied by the system or package manager instead of CPM" OFF)
option(NRGLJUBLJANA_INSTALL_NRGINIT "Install nrginit Mathematica scripts" ON)
option(NRGLJUBLJANA_ENABLE_APP_OPENMP "Enable application-level OpenMP parallel regions; BLAS/LAPACK remains threaded independently" OFF)
option(Build_Tests "Build tests" ON)
set(NRGLJUBLJANA_BUILD_TESTS_ENABLED OFF)
if(Build_Tests AND NOT IS_SUBPROJECT)
  set(NRGLJUBLJANA_BUILD_TESTS_ENABLED ON)
endif()

message(STATUS "Application OpenMP parallel regions: ${NRGLJUBLJANA_ENABLE_APP_OPENMP}")

# Default to Release build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Type of build" FORCE)
endif()
message(STATUS "-------- BUILD-TYPE: ${CMAKE_BUILD_TYPE} --------")

# Export the list of compile-commands into compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Global compiler options
add_compile_options(
  $<$<CONFIG:Debug>:-Og>
  $<$<CONFIG:Debug>:-fno-var-tracking>
  $<$<CONFIG:Debug>:-ggdb3>
)

# Create an Interface target for compiler warnings
if(NOT TARGET project_warnings)
  add_library(project_warnings INTERFACE)
endif()
install(TARGETS project_warnings EXPORT nrgljubljana-targets)
target_compile_options(project_warnings
  INTERFACE
    -Wall
    -Wextra
    -Wnon-virtual-dtor
#    -Wold-style-cast  # because Eigen
    -Wcast-align
    -Wunused
    -Woverloaded-virtual
    -Wpedantic
#    -Wconversion # USEFUL, but too many false flags..
#    -Wsign-conversion # because {fmt}, Eigen
    $<$<CXX_COMPILER_ID:GNU>:-Wduplicated-cond>
    $<$<CXX_COMPILER_ID:GNU>:-Wlogical-op>
    -Wno-sign-compare
    -Wno-unused-parameter # USEFUL, occasionally...
    -Wno-unused-function   # issue in matrix tool
    -Wno-deprecated-declarations  # HighFive using boost/ublas
    $<$<CXX_COMPILER_ID:GNU>:-Wshadow=local>
    $<$<CXX_COMPILER_ID:Clang>:-Wshadow>
    $<$<CXX_COMPILER_ID:Clang>:-Wno-gcc-compat>
)

# ###########
# MPI support

message(STATUS "-------- MPI detection -------------")
unset(MPIEXEC_EXECUTABLE CACHE)
find_package(MPI REQUIRED)

# Backward compatibility for older FindMPI.cmake
if(NOT MPIEXEC_EXECUTABLE)
  set(MPIEXEC_EXECUTABLE ${MPIEXEC} CACHE FILENAME "MPI Executable")
endif()

# Compatibility to Open-MPI 3.0.0: check whether MPI executable has option --oversubscribe and add it
execute_process(COMMAND ${MPIEXEC_EXECUTABLE} --oversubscribe ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ls ${MPIEXEC_POSTFLAGS} RESULT_VARIABLE HAS_NO_OVERSUBSCRIBE OUTPUT_QUIET ERROR_QUIET)
if(NOT HAS_NO_OVERSUBSCRIBE AND NOT "--oversubscribe" IN_LIST MPIEXEC_PREFLAGS)
 list(APPEND MPIEXEC_PREFLAGS --oversubscribe)
 set(MPIEXEC_PREFLAGS ${MPIEXEC_PREFLAGS} CACHE STRING "These flags will be directly before the executable that is being run by mpiexec." FORCE)
endif()
if(MPIEXEC_PREFLAGS)
  list(REMOVE_DUPLICATES MPIEXEC_PREFLAGS)
  set(MPIEXEC_PREFLAGS ${MPIEXEC_PREFLAGS} CACHE STRING "These flags will be directly before the executable that is being run by mpiexec." FORCE)
endif()

message(STATUS "MPI_FOUND: ${MPI_FOUND} VERSION: ${MPI_VERSION}")
message(STATUS "MPI_CXX_COMPILER: ${MPI_CXX_COMPILER}")
message(STATUS "MPIEXEC: ${MPIEXEC}")
message(STATUS "MPIEXEC_EXECUTABLE: ${MPIEXEC_EXECUTABLE}")
message(STATUS "MPIEXEC_NUMPROC_FLAG: ${MPIEXEC_NUMPROC_FLAG}")
message(STATUS "MPIEXEC_PREFLAGS: ${MPIEXEC_PREFLAGS}")
message(STATUS "MPIEXEC_POSTFLAGS: ${MPIEXEC_POSTFLAGS}")

# Create an interface target
add_library(mpi INTERFACE)
target_include_directories(mpi SYSTEM INTERFACE ${MPI_C_INCLUDE_PATH})
target_link_libraries(mpi INTERFACE ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})
separate_arguments(MPI_C_COMPILE_FLAGS) # Convert to list
target_compile_options(mpi INTERFACE ${MPI_C_COMPILE_FLAGS})
install(TARGETS mpi EXPORT nrgljubljana-targets)

# ############
# Dependencies

# HDF5
find_package(HDF5 MODULE REQUIRED COMPONENTS C HL)
add_library(hdf5 INTERFACE)
target_include_directories(hdf5 SYSTEM INTERFACE ${HDF5_INCLUDE_DIRS})
target_link_libraries(hdf5 INTERFACE "${HDF5_LIBRARIES}" ${HDF5_HL_LIBRARIES})
target_compile_options(hdf5 INTERFACE ${HDF5_DEFINITIONS})
if(NOT TARGET HDF5::HDF5)
  add_library(HDF5::HDF5 ALIAS hdf5)
endif()
install(TARGETS hdf5 EXPORT nrgljubljana-targets)

# OpenMP for application-level parallel regions. The default parallelism model is
# BLAS/LAPACK-internal threading; this target is intentionally empty unless the
# explicit opt-in below is enabled.
add_library(openmp INTERFACE)
target_compile_definitions(openmp INTERFACE NRG_ENABLE_APP_OPENMP=$<BOOL:${NRGLJUBLJANA_ENABLE_APP_OPENMP}>)
if(NRGLJUBLJANA_ENABLE_APP_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS CXX)
  target_link_libraries(openmp INTERFACE OpenMP::OpenMP_CXX)
endif()
install(TARGETS openmp EXPORT nrgljubljana-targets)

# Boost
# https://stackoverflow.com/questions/57415206/is-there-a-way-to-get-rid-of-the-new-boost-version-may-have-incorrect-or-missin
set(Boost_NO_WARN_NEW_VERSIONS 1)
set(Boost_NO_BOOST_CMAKE ON)
set(Boost_USE_STATIC_LIBS OFF)
find_package(Boost 1.53 REQUIRED COMPONENTS serialization mpi)
# note: header-only libraries such as multiprecision do not need to be listed
# as components in find_package()

# GSL
find_package(GSL REQUIRED)
add_library(gsl INTERFACE)
target_include_directories(gsl SYSTEM INTERFACE ${GSL_INCLUDE_DIRS})
target_link_libraries(gsl INTERFACE ${GSL_LIBRARIES})
install(TARGETS gsl EXPORT nrgljubljana-targets)

# GMP
find_package(GMP REQUIRED)
install(TARGETS gmp EXPORT nrgljubljana-targets)

if(NRGLJUBLJANA_USE_SYSTEM_DEPS)
  message(STATUS "Using dependencies supplied by the system or package manager")

  if(NRGLJUBLJANA_BUILD_TESTS_ENABLED)
    find_package(GTest CONFIG REQUIRED)
    if(TARGET GTest::gtest AND NOT TARGET gtest)
      add_library(gtest ALIAS GTest::gtest)
    endif()
    if(TARGET GTest::gtest_main AND NOT TARGET gtest_main)
      add_library(gtest_main ALIAS GTest::gtest_main)
    endif()
  endif()

  find_package(Eigen3 3.4 CONFIG REQUIRED)
  set(HIGHFIVE_FIND_HDF5 Off)
  find_package(HighFive CONFIG REQUIRED)
  find_package(fmt 11 CONFIG REQUIRED)
  find_package(range-v3 CONFIG REQUIRED)
else()
  # CPM for downloading dependencies
  include(CMake/CPM.cmake)

  if(NRGLJUBLJANA_BUILD_TESTS_ENABLED)
    CPMAddPackage("gh:google/googletest#v1.15.2")
  endif()

  # Eigen, tries finding in system first
  CPMFindPackage(
    NAME Eigen3
    GITLAB_REPOSITORY libeigen/eigen
    VERSION 3.4.0
    GIT_TAG 3.4.0
    OPTIONS
      "BUILD_TESTING OFF" # Disable tests for Eigen
      "EIGEN_BUILD_TESTING OFF"
      "EIGEN_BUILD_CMAKE_PACKAGE ON"
  )

  # High five library
  CPMAddPackage(
    GITHUB_REPOSITORY highfive-devs/highfive
    GIT_TAG v3.1.0
    OPTIONS
      "HIGHFIVE_BUILD_DOCS OFF" # Disable HighFive documentation
  )

  # fmt
  # Check if fmt is already a target; important when using as an outside library in nrgljubljana_interface
  if (NOT TARGET fmt::fmt)
    CPMAddPackage("gh:fmtlib/fmt#11.2.0")
  endif()

  # range-v3
  CPMAddPackage("gh:ericniebler/range-v3#0.12.0")
endif()

# MPFR
#find_package(MPFR REQUIRED)
#install(TARGETS mpfr EXPORT nrgljubljana-targets)
### Not yet. Long-term goal: promote thermodynamic quantities to multiple-precision variables.

# Mathematica (optional)
if(NRGLJUBLJANA_ENABLE_MATHEMATICA)
  find_package(Mathematica)
  if(Mathematica_FOUND)
    message(STATUS "Mathematica executable ${Mathematica_KERNEL_EXECUTABLE}")
  endif()
endif()

function(nrgljubljana_collect_target_link_items output target)
  set(items)
  if(TARGET ${target})
    foreach(prop IN ITEMS
        INTERFACE_LINK_LIBRARIES
        INTERFACE_LINK_OPTIONS
        INTERFACE_COMPILE_OPTIONS
        IMPORTED_LINK_INTERFACE_LIBRARIES
        IMPORTED_LOCATION
        IMPORTED_LOCATION_RELEASE
        IMPORTED_LOCATION_RELWITHDEBINFO
        IMPORTED_LOCATION_MINSIZEREL
        IMPORTED_LOCATION_DEBUG)
      get_target_property(values ${target} ${prop})
      if(values AND NOT values MATCHES "-NOTFOUND$")
        list(APPEND items ${values})
      endif()
    endforeach()
  endif()
  set(${output} ${items} PARENT_SCOPE)
endfunction()

function(nrgljubljana_detect_openmp_runtime_families output)
  set(families)
  set(sequential_mkl OFF)
  set(mkl_rt OFF)
  foreach(item IN LISTS ARGN)
    if(NOT item)
      continue()
    endif()
    string(TOLOWER "${item}" item_lower)
    if(item_lower MATCHES "mkl_sequential")
      set(sequential_mkl ON)
    endif()
    if(item_lower MATCHES "mkl_rt")
      set(mkl_rt ON)
    endif()
    if(item_lower MATCHES "libgomp|[-]lgomp|mkl_gnu_thread")
      list(APPEND families GNU)
    endif()
    if(item_lower MATCHES "libiomp5|[-]liomp5|libguide|[-]lguide|mkl_intel_thread")
      list(APPEND families Intel)
    endif()
    if(item_lower MATCHES "libomp\\.|[-]lomp|mkl_llvm_thread")
      list(APPEND families LLVM)
    endif()
  endforeach()
  if(families)
    list(REMOVE_DUPLICATES families)
  endif()
  set(${output} ${families} PARENT_SCOPE)
  set(NRGLJUBLJANA_DETECTED_SEQUENTIAL_MKL ${sequential_mkl} PARENT_SCOPE)
  set(NRGLJUBLJANA_DETECTED_MKL_RT ${mkl_rt} PARENT_SCOPE)
endfunction()

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake/BLAS_LAPACK)
# BLAS/Lapack
find_package(LAPACK REQUIRED) # also looks up BLAS
add_library(blas_lapack INTERFACE)
set(NRGLJUBLJANA_HAVE_MKL OFF)
set(NRGLJUBLJANA_HAVE_OPENBLAS OFF)
set(NRGLJUBLJANA_BLAS_VENDOR "unknown")
string(TOLOWER "${BLA_VENDOR};${BLAS_LIBRARIES};${LAPACK_LIBRARIES};${BLAS_LINKER_FLAGS};${LAPACK_LINKER_FLAGS}" NRGLJUBLJANA_LINALG_LINK_TEXT)
if(NRGLJUBLJANA_LINALG_LINK_TEXT MATCHES "mkl")
  set(NRGLJUBLJANA_HAVE_MKL ON)
  set(NRGLJUBLJANA_BLAS_VENDOR "Intel MKL")
  message(STATUS "Using Intel MKL for linear algebra")
  # https://software.intel.com/en-us/forums/intel-math-kernel-library/topic/759670
  # https://software.intel.com/en-us/articles/a-new-linking-model-single-dynamic-library-mkl_rt-since-intel-mkl-103
elseif(NRGLJUBLJANA_LINALG_LINK_TEXT MATCHES "openblas")
  set(NRGLJUBLJANA_HAVE_OPENBLAS ON)
  set(NRGLJUBLJANA_BLAS_VENDOR "OpenBLAS")
  message(STATUS "Using OpenBLAS for linear algebra")
else()
  message(STATUS "Using BLAS/LAPACK vendor: ${NRGLJUBLJANA_BLAS_VENDOR}")
endif()
string(REPLACE "\"" "\\\"" NRGLJUBLJANA_BLAS_VENDOR_ESCAPED "${NRGLJUBLJANA_BLAS_VENDOR}")
target_compile_definitions(blas_lapack INTERFACE
  "NRG_BLAS_VENDOR=\"${NRGLJUBLJANA_BLAS_VENDOR_ESCAPED}\""
  NRG_HAVE_MKL=$<BOOL:${NRGLJUBLJANA_HAVE_MKL}>
  NRG_HAVE_OPENBLAS=$<BOOL:${NRGLJUBLJANA_HAVE_OPENBLAS}>
)
message(STATUS "BLAS_LIBRARIES=${BLAS_LIBRARIES}")
message(STATUS "BLAS_LINKER_FLAGS=${BLAS_LINKER_FLAGS}")
message(STATUS "LAPACK_LIBRARIES=${LAPACK_LIBRARIES}")
message(STATUS "LAPACK_LINKER_FLAGS=${LAPACK_LINKER_FLAGS}")

set(NRGLJUBLJANA_OPENMP_RUNTIME_LINK_ITEMS
  ${LAPACK_LIBRARIES}
  ${BLAS_LIBRARIES}
  ${LAPACK_LINKER_FLAGS}
  ${BLAS_LINKER_FLAGS}
  ${MPI_C_LIBRARIES}
  ${MPI_CXX_LIBRARIES}
  ${MPI_C_LINK_FLAGS}
  ${MPI_CXX_LINK_FLAGS}
  ${HDF5_LIBRARIES}
  ${HDF5_HL_LIBRARIES}
  ${GSL_LIBRARIES}
)
nrgljubljana_collect_target_link_items(NRGLJUBLJANA_LAPACK_TARGET_ITEMS LAPACK::LAPACK)
nrgljubljana_collect_target_link_items(NRGLJUBLJANA_BLAS_TARGET_ITEMS BLAS::BLAS)
nrgljubljana_collect_target_link_items(NRGLJUBLJANA_MPI_TARGET_ITEMS mpi)
nrgljubljana_collect_target_link_items(NRGLJUBLJANA_HDF5_TARGET_ITEMS hdf5)
nrgljubljana_collect_target_link_items(NRGLJUBLJANA_GSL_TARGET_ITEMS gsl)
nrgljubljana_collect_target_link_items(NRGLJUBLJANA_BOOST_MPI_TARGET_ITEMS Boost::mpi)
nrgljubljana_collect_target_link_items(NRGLJUBLJANA_BOOST_SERIALIZATION_TARGET_ITEMS Boost::serialization)
list(APPEND NRGLJUBLJANA_OPENMP_RUNTIME_LINK_ITEMS
  ${NRGLJUBLJANA_LAPACK_TARGET_ITEMS}
  ${NRGLJUBLJANA_BLAS_TARGET_ITEMS}
  ${NRGLJUBLJANA_MPI_TARGET_ITEMS}
  ${NRGLJUBLJANA_HDF5_TARGET_ITEMS}
  ${NRGLJUBLJANA_GSL_TARGET_ITEMS}
  ${NRGLJUBLJANA_BOOST_MPI_TARGET_ITEMS}
  ${NRGLJUBLJANA_BOOST_SERIALIZATION_TARGET_ITEMS}
)
if(NRGLJUBLJANA_ENABLE_APP_OPENMP)
  nrgljubljana_collect_target_link_items(NRGLJUBLJANA_OPENMP_TARGET_ITEMS OpenMP::OpenMP_CXX)
  list(APPEND NRGLJUBLJANA_OPENMP_RUNTIME_LINK_ITEMS ${OpenMP_CXX_LIBRARIES} ${OpenMP_CXX_FLAGS} ${NRGLJUBLJANA_OPENMP_TARGET_ITEMS})
endif()
nrgljubljana_detect_openmp_runtime_families(NRGLJUBLJANA_OPENMP_RUNTIME_FAMILIES ${NRGLJUBLJANA_OPENMP_RUNTIME_LINK_ITEMS})
if(NRGLJUBLJANA_DETECTED_SEQUENTIAL_MKL)
  message(FATAL_ERROR "Sequential MKL was selected, but NRG Ljubljana expects threaded BLAS/LAPACK. Use a threaded MKL BLA_VENDOR, mkl_rt with a threaded MKL_THREADING_LAYER, or OpenBLAS.")
endif()
if(NRGLJUBLJANA_OPENMP_RUNTIME_FAMILIES)
  message(STATUS "OpenMP runtime families visible in BLAS/LAPACK/app link line: ${NRGLJUBLJANA_OPENMP_RUNTIME_FAMILIES}")
else()
  message(STATUS "No OpenMP runtime family visible in the link line (normal for pthread OpenBLAS or opaque mkl_rt); runtime startup reporting will inspect loaded libraries")
endif()
list(LENGTH NRGLJUBLJANA_OPENMP_RUNTIME_FAMILIES NRGLJUBLJANA_OPENMP_RUNTIME_FAMILY_COUNT)
if(NRGLJUBLJANA_OPENMP_RUNTIME_FAMILY_COUNT GREATER 1)
  message(FATAL_ERROR "More than one OpenMP runtime family appears in the link line: ${NRGLJUBLJANA_OPENMP_RUNTIME_FAMILIES}. Do not mix GNU libgomp, Intel libiomp5/libguide, and LLVM libomp. Disable application OpenMP or select a BLAS/LAPACK threading layer compatible with the compiler OpenMP runtime.")
endif()
if(NRGLJUBLJANA_ENABLE_APP_OPENMP AND NRGLJUBLJANA_DETECTED_MKL_RT)
  message(WARNING "Application OpenMP is enabled together with MKL's mkl_rt dispatcher. CMake cannot prove which MKL threading layer will be loaded at runtime; set MKL_THREADING_LAYER to the OpenMP runtime family used by the compiler, or disable NRGLJUBLJANA_ENABLE_APP_OPENMP.")
endif()
if(TARGET LAPACK::LAPACK)
  target_link_libraries(blas_lapack INTERFACE LAPACK::LAPACK)
else()
  message(WARNING "LAPACK::LAPACK target not provided by FindLAPACK; using LAPACK_LIBRARIES/BLAS_LIBRARIES directly")
  target_link_libraries(blas_lapack INTERFACE ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
endif()
install(TARGETS blas_lapack EXPORT nrgljubljana-targets)

include(CheckFunctionExists)
# Check for availability of cblas interface
set(CMAKE_REQUIRED_LIBRARIES ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
set(CMAKE_REQUIRED_FLAGS ${LAPACK_LINKER_FLAGS} ${BLAS_LINKER_FLAGS})
check_function_exists(cblas_dgemm HAVE_CBLAS)
unset(CMAKE_REQUIRED_LIBRARIES)
unset(CMAKE_REQUIRED_FLAGS)
if(HAVE_CBLAS)
  message(STATUS "cblas interface found in blas library")
  set(CBLAS_WORKAROUND OFF)
else()
  message(WARNING "cblas interface not found; will provide a workaround")
  set(CBLAS_WORKAROUND ON)
endif()

# Symmetry types to compile in
option(SYM_MORE "Compile in an extended set of symmetry types" ON)
option(SYM_ALL  "Compile in the full set of symmetry types (long compilation time)" ON)

if(SYM_MORE)
  message(STATUS "Building an extended set of symmetry types")
endif()
if(SYM_ALL)
  message(STATUS "Building the full set of symmetry types")
endif()

# Test suites
option(TEST_LONG "Enable tests with long runtime" OFF)

# Dynamic Analyzer Checks

option(ASAN  "Compile library and executables with LLVM Address Sanitizer" OFF)
option(UBSAN "Compile library and executables with LLVM Undefined Behavior Sanitizer" OFF)

if(ASAN)
  if(NOT TARGET asan)
    find_package(sanitizer REQUIRED "asan")
  endif()
  install(TARGETS asan EXPORT nrgljubljana-targets)
endif()
if(UBSAN)
  if(NOT TARGET ubsan)
    find_package(sanitizer REQUIRED "ubsan")
  endif()
  install(TARGETS ubsan EXPORT nrgljubljana-targets)
endif()

# Resolve Clang Linktime Problems
# CMake will adjust any linker flags from '-L path_to/mylib.so' to -lmylib
# if the proper mylib.so is automatically found by the linker, i.e.
# the directory comes first in LIBRARY_PATH.
# The clang linker however ignores LIBRARY_PATH.
# We thus explicitly add the content of LIBRARY_PATH to the LDFLAGS
# FIXME For future cmake versions we should populate the
# INTERFACE_LINK_DIRECTORIES of the triqs target
# ---------------------------------
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" AND DEFINED ENV{LIBRARY_PATH})
  string(REPLACE ":" ";" LINK_DIRS $ENV{LIBRARY_PATH})
  foreach(dir ${LINK_DIRS})
    string(APPEND CMAKE_SHARED_LINKER_FLAGS " -L${dir}")
    string(APPEND CMAKE_MODULE_LINKER_FLAGS " -L${dir}")
    string(APPEND CMAKE_EXE_LINKER_FLAGS " -L${dir}")
  endforeach()
endif()

# Workaround for ABI notes on aarch64 platforms
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64|arm64)$"
   AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
     add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-Wno-psabi>)
endif()

# #############
# Build Project

# Build and install the nrgljubljana library
add_subdirectory(c++)

# Tools
add_subdirectory(tools)

# Tests
if(NRGLJUBLJANA_BUILD_TESTS_ENABLED)
  enable_testing()
  add_subdirectory(test)
endif()

# Build the documentation
option(Build_Documentation "Build documentation" OFF)
if(Build_Documentation AND NOT IS_SUBPROJECT)
  message(STATUS "Documentation Build enabled")
  add_subdirectory(doc)
endif()

# Additional configuration files
add_subdirectory(share)

# Mathematica part of NRG Ljubljana (nrginit)
if(NRGLJUBLJANA_INSTALL_NRGINIT)
  add_subdirectory(nrginit)
else()
  message(STATUS "Skipping nrginit installation (NRGLJUBLJANA_INSTALL_NRGINIT=OFF)")
endif()

# #############
# Debian Package

option(BUILD_DEBIAN_PACKAGE "Build a deb package" OFF)
if(BUILD_DEBIAN_PACKAGE AND NOT IS_SUBPROJECT)
  if(NOT CMAKE_INSTALL_PREFIX STREQUAL "/usr")
    message(FATAL_ERROR "CMAKE_INSTALL_PREFIX must be /usr for packaging")
  endif()
  set(CPACK_PACKAGE_NAME nrgljubljana)
  set(CPACK_GENERATOR "DEB")
  set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
  set(CPACK_PACKAGE_CONTACT "https://github.com/rokzitko/nrgljubljana")
  execute_process(COMMAND dpkg --print-architecture OUTPUT_VARIABLE CMAKE_DEBIAN_PACKAGE_ARCHITECTURE OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
  set(CPACK_DEBIAN_PACKAGE_GENERATE_SHLIBS ON)
  include(CPack)
endif()
