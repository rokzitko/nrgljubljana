macro(add_tool exec_name)
  add_executable(${exec_name} ${exec_name}/${exec_name}.cc)

  # Enable warnings
  target_link_libraries(${exec_name} PRIVATE project_warnings)

  # C++ Standard
  target_compile_features(${exec_name} PRIVATE cxx_std_20)

  # Configure compilation
  target_compile_options(${exec_name} PRIVATE -fPIC ${NRG_CXX_OPT})
  target_include_directories(${exec_name} PRIVATE 
                     $<BUILD_INTERFACE:${CURRENT_SOURCE_DIR}/${exec_name}>
                     $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/c++>
                     $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/c++/include>
                     $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/c++/include/range-v3>)
  target_compile_definitions(${exec_name} PRIVATE
# (not used)				NRGLJUBLJANA_GIT_HASH=${PROJECT_GIT_HASH}
				$<$<CONFIG:Debug>:NRGLJUBLJANA_DEBUG>
			  )

  # Link dependencies
  target_link_libraries(${exec_name} PRIVATE openmp boost blas_lapack gmp gsl dl hdf5
        $<$<BOOL:${ASAN}>:asan>
        $<$<BOOL:${UBSAN}>:ubsan>
  )

  # Install executables
  install(TARGETS ${exec_name} EXPORT nrgljubljana-targets DESTINATION bin)
endmacro()

set(all_executables adapt binavg broaden bw diag h5write hilb intavg integ kk matrix mats nrgchain resample specmoments tdavg unitary)
foreach(exec ${all_executables})
  add_tool(${exec})
endforeach()

# Additional source file for 'matrix' tool
target_sources(matrix PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/matrix/parser.cc)
