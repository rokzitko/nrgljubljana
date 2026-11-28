# The symmetry types are listed in symmetries.def. Each one compiled in gets a factory source generated from
# sym-instance.cc.in; SYM_MORE and SYM_ALL add the corresponding build sets to BASIC.
set(sources nrg-lib.cc)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS symmetries.def sym-instance.cc.in)
file(STRINGS symmetries.def symmetry_entries REGEX "^NRG_SYMMETRY\\(")
foreach(entry ${symmetry_entries})
  string(REGEX REPLACE "^NRG_SYMMETRY\\(([^)]*)\\).*$" "\\1" fields "${entry}")
  string(REPLACE " " "" fields "${fields}")
  string(REPLACE "," ";" fields "${fields}")
  list(GET fields 0 NRG_SYM_NAME)
  list(GET fields 1 NRG_SYM_IMPL)
  list(GET fields 2 set)
  list(GET fields 8 complex_only)
  if(set STREQUAL "BASIC" OR (set STREQUAL "MORE" AND SYM_MORE) OR (set STREQUAL "ALL" AND SYM_ALL))
    if(complex_only STREQUAL "true")
      set(NRG_SYM_COMPLEX_ONLY 1)
    else()
      set(NRG_SYM_COMPLEX_ONLY 0)
    endif()
    configure_file(sym-instance.cc.in ${CMAKE_CURRENT_BINARY_DIR}/sym-${NRG_SYM_NAME}.cc @ONLY)
    list(APPEND sources ${CMAKE_CURRENT_BINARY_DIR}/sym-${NRG_SYM_NAME}.cc)
  endif()
endforeach()

add_library(nrgljubljana_c SHARED ${sources})

# C++ Standard
target_compile_features(nrgljubljana_c PUBLIC cxx_std_20)

# Enable warnings
target_link_libraries(nrgljubljana_c     PRIVATE project_warnings)

# Set library version and PIC flag
set_target_properties(nrgljubljana_c PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  VERSION ${PROJECT_VERSION}
)

option(GPROF "Enable profiling using gprof" OFF)

# Configure compilation
target_compile_options(nrgljubljana_c     PRIVATE ${NRG_CXX_OPT} $<$<BOOL:${GPROF}>:-pg>)
nrgljubljana_apply_strict_floating_point(nrgljubljana_c)
target_include_directories(nrgljubljana_c     PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/c++>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/c++/symmetry>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/c++/symmetry/coefnew>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/c++/include>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/c++/include/range-v3>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/tools>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/nrgljubljana>)
target_compile_definitions(nrgljubljana_c PUBLIC
  NRG_SYM_BASIC
  $<$<BOOL:${SYM_MORE}>:NRG_SYM_MORE>
  $<$<BOOL:${SYM_ALL}>:NRG_SYM_ALL>
)

# Link dependencies
target_link_libraries(nrgljubljana_c
  PUBLIC
    Boost::headers
    Boost::mpi
    Boost::serialization
    mpi
    $<BUILD_INTERFACE:openmp>
    $<BUILD_INTERFACE:cuda_solver>
    $<BUILD_INTERFACE:blas_lapack>
    $<BUILD_INTERFACE:gmp>
    $<BUILD_INTERFACE:dl>
    $<BUILD_INTERFACE:hdf5>
    $<BUILD_INTERFACE:HighFive>
    $<BUILD_INTERFACE:fmt::fmt-header-only>
    $<BUILD_INTERFACE:range-v3>
    $<BUILD_INTERFACE:Eigen3::Eigen>
    $<BUILD_INTERFACE:$<$<BOOL:${ASAN}>:asan>>
    $<BUILD_INTERFACE:$<$<BOOL:${UBSAN}>:ubsan>>
)

target_link_options(nrgljubljana_c PRIVATE $<$<BOOL:${GPROF}>:-pg>)

# Compile nrg executable
add_executable(nrg nrg.cc)
target_link_libraries(nrg PRIVATE nrgljubljana_c)
target_include_directories(nrg PRIVATE $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/tools/generated>)
target_link_options(nrg PRIVATE $<$<BOOL:${GPROF}>:-pg>)
target_compile_options(nrg PRIVATE ${NRG_CXX_OPT} $<$<BOOL:${GPROF}>:-pg>)
nrgljubljana_apply_strict_floating_point(nrg)
string(REPLACE "\"" "\\\"" NRGLJUBLJANA_CXX_COMPILER_ID_ESCAPED "${CMAKE_CXX_COMPILER_ID}")
string(REPLACE "\"" "\\\"" NRGLJUBLJANA_CXX_COMPILER_VERSION_ESCAPED "${CMAKE_CXX_COMPILER_VERSION}")
string(REPLACE "\"" "\\\"" NRGLJUBLJANA_SYSTEM_PROCESSOR_ESCAPED "${CMAKE_SYSTEM_PROCESSOR}")
string(REPLACE "\"" "\\\"" NRGLJUBLJANA_SYSTEM_NAME_ESCAPED "${CMAKE_SYSTEM_NAME}")
string(REPLACE "\"" "\\\"" NRGLJUBLJANA_SYSTEM_VERSION_ESCAPED "${CMAKE_SYSTEM_VERSION}")
target_compile_definitions(nrg PRIVATE
  "NRG_BUILD_CXX_COMPILER_ID=\"${NRGLJUBLJANA_CXX_COMPILER_ID_ESCAPED}\""
  "NRG_BUILD_CXX_COMPILER_VERSION=\"${NRGLJUBLJANA_CXX_COMPILER_VERSION_ESCAPED}\""
  "NRG_BUILD_SYSTEM_PROCESSOR=\"${NRGLJUBLJANA_SYSTEM_PROCESSOR_ESCAPED}\""
  "NRG_BUILD_SYSTEM_NAME=\"${NRGLJUBLJANA_SYSTEM_NAME_ESCAPED}\""
  "NRG_BUILD_SYSTEM_VERSION=\"${NRGLJUBLJANA_SYSTEM_VERSION_ESCAPED}\""
  $<$<BOOL:${NRGLJUBLJANA_ENABLE_FP_TRAPS}>:NRGLJUBLJANA_ENABLE_FP_TRAPS=1>
)
set_target_properties(nrg PROPERTIES BUILD_RPATH "${PROJECT_BINARY_DIR}/c++")

# Install library and headers
install(TARGETS nrgljubljana_c EXPORT nrgljubljana-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(FILES nrg-lib.hpp workdir.hpp portabil.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nrgljubljana)
install(TARGETS nrg RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# ========= Static Analyzer Checks ==========

option(ANALYZE_SOURCES "Run static analyzer checks if found (clang-tidy, cppcheck)" OFF)
if(ANALYZE_SOURCES)

  # Locate static analyzer tools
  find_program(CPPCHECK_EXECUTABLE NAMES "cppcheck" PATHS ENV PATH)
  find_program(CLANG_TIDY_EXECUTABLE NAMES "clang-tidy" PATHS ENV PATH)

  # Run clang-tidy if found
  if(CLANG_TIDY_EXECUTABLE)
    message(STATUS "clang-tidy found: ${CLANG_TIDY_EXECUTABLE}")
    set_target_properties(nrgljubljana_c PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_EXECUTABLE}")
  else()
    message(STATUS "clang-tidy not found in $PATH. Please consider installing clang-tidy for additional checks!")
  endif()

  # Run cppcheck if found
  if(CPPCHECK_EXECUTABLE)
    message(STATUS "cppcheck found: ${CPPCHECK_EXECUTABLE}")
    add_custom_command(
	TARGET nrgljubljana_c POST_BUILD
	COMMAND ${CPPCHECK_EXECUTABLE}
	  --enable=warning,style,performance,portability
	  --std=c++20
	  --template=gcc
	  --verbose
	  --force
	  --quiet
    --suppressions-list=suppressions.txt
    --inline-suppr
    ${sources}
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
  else()
    message(STATUS "cppcheck not found in $PATH. Please consider installing cppcheck for additional checks!")
  endif()

endif()
