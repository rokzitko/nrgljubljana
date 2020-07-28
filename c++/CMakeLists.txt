set(sources nrg-lib.cc)
add_library(nrgljubljana_c     SHARED ${sources})
add_library(nrgljubljanacmpl_c SHARED ${sources})

# Enable warnings
target_link_libraries(nrgljubljana_c     PRIVATE project_warnings)
target_link_libraries(nrgljubljanacmpl_c PRIVATE project_warnings)

# Configure compilation
target_compile_options(nrgljubljana_c     PRIVATE -fPIC -std=c++14)
target_compile_options(nrgljubljanacmpl_c PRIVATE -fPIC -std=c++14)
target_include_directories(nrgljubljana_c     PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/c++>)
target_include_directories(nrgljubljanacmpl_c PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/c++>)
target_include_directories(nrgljubljana_c     PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/c++/coefnew>)
target_include_directories(nrgljubljanacmpl_c PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/c++/coefnew>)
target_compile_definitions(nrgljubljana_c PUBLIC
				NRG_REAL
				$<$<BOOL:${NRG_MPI}>:NRG_MPI>
				$<$<CONFIG:Debug>:NRGLJUBLJANA_DEBUG>
				NRG_SYM_BASIC
				$<$<BOOL:${SYM_MORE}>:NRG_SYM_MORE>
				$<$<BOOL:${SYM_ALL}>:NRG_SYM_ALL>
				$<$<BOOL:${CBLAS_WORKAROUND}>:CBLAS_WORKAROUND>
)
target_compile_definitions(nrgljubljanacmpl_c PUBLIC
				NRG_COMPLEX
				$<$<BOOL:${NRG_MPI}>:NRG_MPI>
				$<$<CONFIG:Debug>:NRGLJUBLJANA_DEBUG>
				NRG_SYM_BASIC
				$<$<BOOL:${SYM_MORE}>:NRG_SYM_MORE>
				$<$<BOOL:${SYM_ALL}>:NRG_SYM_ALL>
				$<$<BOOL:${CBLAS_WORKAROUND}>:CBLAS_WORKAROUND>
)

# Link dependencies
target_link_libraries(nrgljubljana_c     PUBLIC openmp boost blas_lapack gmp dl
  $<$<BOOL:${ASAN}>:asan>
  $<$<BOOL:${UBSAN}>:ubsan>
  $<$<BOOL:${NRG_MPI}>:mpi>
)
target_link_libraries(nrgljubljanacmpl_c PUBLIC openmp boost blas_lapack gmp dl
  $<$<BOOL:${ASAN}>:asan>
  $<$<BOOL:${UBSAN}>:ubsan>
  $<$<BOOL:${NRG_MPI}>:mpi>
)

# Compile nrg executable
add_executable(nrg nrg.cc)
target_link_libraries(nrg PRIVATE nrgljubljana_c
  $<$<BOOL:${NRG_MPI}>:mpi>
)
target_compile_definitions(nrg PUBLIC
  $<$<BOOL:${NRG_MPI}>:NRG_MPI>
  NRG_REAL
)

add_executable(nrgcmpl nrg.cc)
target_link_libraries(nrgcmpl PRIVATE nrgljubljanacmpl_c    # cmpl !!
  $<$<BOOL:${NRG_MPI}>:mpi>
)
target_compile_definitions(nrgcmpl PUBLIC
  $<$<BOOL:${NRG_MPI}>:NRG_MPI>
  NRG_COMPLEX
)

# Install library and headers
install(TARGETS nrgljubljana_c nrgljubljanacmpl_c EXPORT nrgljubljana-targets DESTINATION lib)
install(TARGETS nrg nrgcmpl DESTINATION bin)

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
	TARGET nrgljubljana_c
	COMMAND ${CPPCHECK_EXECUTABLE}
	  --enable=warning,style,performance,portability
	  --std=c++14
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
