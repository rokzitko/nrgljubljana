# https://cmake.cmake.narkive.com/huyio5CB/using-m4-as-a-additional-preprocessor-for-c-c
# find_program(M4_EXECUTABLE m4 DOC "The M4 macro processor")
# file(GLOB m4_files "nrg-recalc-*.m4")
# add_custom_command(OUTPUT .. COMMAND m4 ARGS )

set(sources nrg-lib.cc sym-QS.cc sym-QSZ.cc)
if(SYM_MORE)
  list(APPEND sources sym-ISO.cc sym-ISOSZ.cc sym-SPSU2.cc sym-SPU1.cc)
endif()
if(SYM_ALL)
  list(APPEND sources sym-DBLQSZ.cc sym-DBLSU2.cc sym-DBLISOSZ.cc sym-ISOLR.cc sym-ISOSZLR.cc sym-NONE.cc sym-P.cc sym-PP.cc sym-SL.cc sym-SL3.cc sym-SPSU2LR.cc
  sym-SPSU2T.cc sym-SPU1LR.cc sym-SU2.cc sym-QSLR.cc sym-QST.cc sym-QSTZ.cc sym-QSZTZ.cc sym-QSZLR.cc sym-QJ.cc
  sym-U1.cc sym-SPSU2C3.cc sym-QSC3.cc)
endif()

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

option(ANALYZE_SOURCES OFF "Run static analyzer checks if found (clang-tidy, cppcheck)")
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
