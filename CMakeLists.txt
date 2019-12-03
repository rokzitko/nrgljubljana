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

cmake_minimum_required(VERSION 3.0.2 FATAL_ERROR)
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

# ############
# Define Project
project(NRGLJUBLJANA VERSION 2019.11 LANGUAGES C CXX)
get_directory_property(IS_SUBPROJECT PARENT_DIRECTORY)

# ############
# Load CPP2PY
# find_package(Cpp2Py 1.6 REQUIRED)

# Get the git hash & print status
execute_process(COMMAND git rev-parse HEAD WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_VARIABLE PROJECT_GIT_HASH OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "${PROJECT_NAME} version : ${PROJECT_VERSION}")
message(STATUS "${PROJECT_NAME} Git hash: ${PROJECT_GIT_HASH}")

# Assert that Install directory is given and invalid.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT OR (NOT IS_ABSOLUTE ${CMAKE_INSTALL_PREFIX}))
  message(FATAL_ERROR "No install prefix given (or invalid)")
endif()
message(STATUS "-------- CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX} --------")

# ############
# Options

# Make additional Find Modules available
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/share/cmake/Modules)

# FindMathematica
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake/Mathematica)

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
    -Wpedantic
    -Wno-sign-compare
    -Wno-extra-semi
    -Wno-unused-parameter
    -Wno-vla
    -Wno-vla-extension
    -Wno-unused-variable
#    $<$<CXX_COMPILER_ID:GNU>:-Wshadow=local>
#    $<$<CXX_COMPILER_ID:Clang>:-Wshadow>
    $<$<CXX_COMPILER_ID:Clang>:-Wno-gcc-compat>
)

# ############
# Dependencies

# OpenMP
find_package(OpenMP REQUIRED)
add_library(openmp INTERFACE)
target_compile_options(openmp INTERFACE ${OpenMP_CXX_FLAGS})
target_link_libraries(openmp INTERFACE ${OpenMP_CXX_FLAGS})
install(TARGETS openmp EXPORT nrgljubljana-targets)

# Boost
set(Boost_NO_BOOST_CMAKE ON)
set(Boost_USE_STATIC_LIBS OFF)
find_package(Boost 1.53 REQUIRED COMPONENTS serialization)
add_library(boost INTERFACE)
target_include_directories(boost SYSTEM INTERFACE ${Boost_INCLUDE_DIRS})
target_link_libraries(boost INTERFACE ${Boost_LIBRARIES})
install(TARGETS boost EXPORT nrgljubljana-targets)

# BLAS/Lapack
find_package(LAPACK REQUIRED)
add_library(blas_lapack INTERFACE)
target_compile_options(blas_lapack INTERFACE ${LAPACK_LINKER_FLAGS})
target_link_libraries(blas_lapack INTERFACE ${LAPACK_LIBRARIES})
install(TARGETS blas_lapack EXPORT nrgljubljana-targets)

# GMP
find_package(GMP REQUIRED)
install(TARGETS gmp EXPORT nrgljubljana-targets)

# Mathematica (optional)
find_package(Mathematica)
if(Mathematica_FOUND)
  message(STATUS "Mathematica executable ${Mathematica_KERNEL_EXECUTABLE}")
endif()

# Dynamic Analyzer Checks

option(ASAN OFF "Compile library and executables with LLVM Address Sanitizer")
option(UBSAN OFF "Compile library and executables with LLVM Undefined Behavior Sanitizer")

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

# #############
# Build Project

# Build and install the nrgljubljana library
add_subdirectory(c++)

# Tools
add_subdirectory(tools)

# Tests
option(Build_Tests "Build tests" ON)
if(Build_Tests AND NOT IS_SUBPROJECT)
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
add_subdirectory(nrginit)

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

