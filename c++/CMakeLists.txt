# https://cmake.cmake.narkive.com/huyio5CB/using-m4-as-a-additional-preprocessor-for-c-c
# find_program(M4_EXECUTABLE m4 DOC "The M4 macro processor")
# file(GLOB m4_files "nrg-recalc-*.m4")
# add_custom_command(OUTPUT .. COMMAND m4 ARGS )

set(sources nrg-lib.cc sym-QS.cc sym-QSZ.cc)
if(SYM_MORE)
  list(APPEND sources sym-ISO.cc sym-ISOSZ.cc sym-SPSU2.cc sym-SPU1.cc)
endif()
if(SYM_ALL)
  list(APPEND sources sym-DBLSU2.cc sym-DBLISOSZ.cc sym-ISOLR.cc sym-ISOSZLR.cc sym-NONE.cc sym-P.cc sym-PP.cc sym-SL.cc sym-SL3.cc sym-SPSU2LR.cc
  sym-SPSU2T.cc sym-SPU1LR.cc sym-SU2.cc sym-QSLR.cc sym-QST.cc sym-QSTZ.cc sym-QSZTZ.cc sym-QSZLR.cc sym-QJ.cc
  sym-U1.cc sym-SPSU2C3.cc sym-QSC3.cc)
endif()

add_library(nrgljubljana_c  SHARED ${sources})

# Enable warnings
target_link_libraries(nrgljubljana_c     PRIVATE project_warnings)

option(GPROF "gprof" OFF)
# Configure compilation
target_compile_options(nrgljubljana_c     PRIVATE ${NRG_CXX_OPT} $<$<BOOL:${GPROF}>:-pg>)
target_include_directories(nrgljubljana_c     PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/c++>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/c++/symmetry>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/c++/symmetry/coefnew>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/c++/include>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/tools>)
target_compile_definitions(nrgljubljana_c PUBLIC
  $<$<CONFIG:Debug>:NRGLJUBLJANA_DEBUG>
  NRG_SYM_BASIC
  $<$<BOOL:${SYM_MORE}>:NRG_SYM_MORE>
  $<$<BOOL:${SYM_ALL}>:NRG_SYM_ALL>
  $<$<BOOL:${CBLAS_WORKAROUND}>:CBLAS_WORKAROUND>
)

# Link dependencies
target_link_libraries(nrgljubljana_c     PUBLIC openmp boost blas_lapack gmp dl mpi hdf5
  $<$<BOOL:${ASAN}>:asan>
  $<$<BOOL:${UBSAN}>:ubsan>
  PRIVATE
  $<$<BOOL:${GPROF}>:-pg>
)

# Compile nrg executable
add_executable(nrg nrg.cc)
target_link_libraries(nrg PRIVATE nrgljubljana_c $<$<BOOL:${GPROF}>:-pg>)
target_compile_options(nrg PRIVATE ${NRG_CXX_OPT} $<$<BOOL:${GPROF}>:-pg>)

# Install library and headers
install(TARGETS nrgljubljana_c EXPORT nrgljubljana-targets DESTINATION lib)
install(TARGETS nrg DESTINATION bin)

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

  # TODO: complex version?

  # Run cppcheck if found
  if(CPPCHECK_EXECUTABLE)
    message(STATUS "cppcheck found: ${CPPCHECK_EXECUTABLE}")
    add_custom_command(
	TARGET nrgljubljana_c
	COMMAND ${CPPCHECK_EXECUTABLE}
	  --enable=warning,style,performance,portability
	  --std=c++17
	  --template=gcc
	  --verbose
	  --force
	  --quiet
    --suppressions-list=suppressions.txt
    --inline-suppr
#    --suppress=*:*bindings*/*
#    --suppress=*:*cblas*
#    --suppress=*:*coewfnew*
    ${sources}
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
  else()
    message(STATUS "cppcheck not found in $PATH. Please consider installing cppcheck for additional checks!")
  endif()

endif()
