set(sources nrg-lib.cc)
add_library(nrgljubljana_c STATIC ${sources})

# Enable warnings
target_link_libraries(nrgljubljana_c PRIVATE project_warnings)

# Configure compilation
target_compile_options(nrgljubljana_c PRIVATE -fPIC -std=c++14)
target_include_directories(nrgljubljana_c PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/c++>)
target_include_directories(nrgljubljana_c PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/c++/coefnew>)
target_compile_definitions(nrgljubljana_c PUBLIC
# (Not used)			NRGLJUBLJANA_GIT_HASH=${PROJECT_GIT_HASH}
				NRG_REAL
				$<$<CONFIG:Debug>:NRGLJUBLJANA_DEBUG>
				NRG_SYM_BASIC
				$<$<BOOL:${SYM_MORE}>:NRG_SYM_MORE>
				$<$<BOOL:${SYM_ALL}>:NRG_SYM_ALL>
				$<$<BOOL:${CBLAS_WORKAROUND}>:CBLAS_WORKAROUND>
)

# Link dependencies
target_link_libraries(nrgljubljana_c PUBLIC openmp boost blas_lapack gmp
  $<$<BOOL:${ASAN}>:asan>
  $<$<BOOL:${UBSAN}>:ubsan>
)

# Compile nrg executable
add_executable(nrg nrg.cc)
target_link_libraries(nrg PRIVATE nrgljubljana_c)

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
