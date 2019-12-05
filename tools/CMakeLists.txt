# GSL
find_package(GSL REQUIRED)
add_library(gsl INTERFACE)
target_include_directories(gsl SYSTEM INTERFACE ${GSL_INCLUDE_DIRS})
target_link_libraries(gsl INTERFACE ${GSL_LIBRARIES})
install(TARGETS gsl EXPORT nrgljubljana-targets)

macro(add_tool exec_name)
  add_executable(${exec_name} ${exec_name}/${exec_name}.cc)

  # Enable warnings
  target_link_libraries(${exec_name} PRIVATE project_warnings)

  # Configure compilation
  target_compile_options(${exec_name} PRIVATE -fPIC -std=c++17)
  target_include_directories(${exec_name} PRIVATE $<BUILD_INTERFACE:${CURRENT_SOURCE_DIR}/${exec_name}>)
  target_compile_definitions(${exec_name} PRIVATE
				NRGLJUBLJANA_GIT_HASH=${PROJECT_GIT_HASH}
				$<$<CONFIG:Debug>:NRGLJUBLJANA_DEBUG>
			  )

  # Link dependencies
  target_link_libraries(${exec_name} PRIVATE openmp boost blas_lapack gmp gsl
        $<$<BOOL:${ASAN}>:asan>
        $<$<BOOL:${UBSAN}>:ubsan>
  )

  # Install executables
  install(TARGETS ${exec_name} EXPORT nrgljubljana-targets DESTINATION bin)
endmacro()

set(all_executables adapt binavg broaden bw diag hilb intavg integ kk matrix mats nrgchain resample specmoments tdavg unitary)
foreach(exec ${all_executables})
  add_tool(${exec})
endforeach()

# Additional source file for 'matrix' tool
target_sources(matrix PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/matrix/parser.cc)
