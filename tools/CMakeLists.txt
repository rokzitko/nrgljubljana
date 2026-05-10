macro(add_tool exec_name)
  add_executable(${exec_name} ${exec_name}/${exec_name}.cc)

  # Enable warnings
  target_link_libraries(${exec_name} PRIVATE project_warnings)

  # C++ Standard
  target_compile_features(${exec_name} PRIVATE cxx_std_20)

  # Configure compilation
  target_compile_options(${exec_name} PRIVATE ${NRG_CXX_OPT})
  nrgljubljana_apply_strict_floating_point(${exec_name})
  target_include_directories(${exec_name} PRIVATE 
                     $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${exec_name}>
                     $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/c++>
                     $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/c++/include>
                     $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/c++/include/range-v3>)
  target_compile_definitions(${exec_name} PRIVATE)
  set_target_properties(${exec_name} PROPERTIES POSITION_INDEPENDENT_CODE ON BUILD_RPATH "${PROJECT_BINARY_DIR}/c++")

  # Link dependencies
  target_link_libraries(${exec_name} PRIVATE openmp Boost::boost Boost::mpi Boost::serialization blas_lapack gmp gsl dl hdf5 HighFive fmt::fmt-header-only range-v3 Eigen3::Eigen
        $<$<BOOL:${ASAN}>:asan>
        $<$<BOOL:${UBSAN}>:ubsan>
  )

  # Install executables
  install(TARGETS ${exec_name} EXPORT nrgljubljana-targets RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endmacro()

add_library(nrgchain_core OBJECT nrgchain/nrgchain.cc)
target_link_libraries(nrgchain_core PRIVATE project_warnings)
target_compile_features(nrgchain_core PRIVATE cxx_std_20)
target_compile_options(nrgchain_core PRIVATE ${NRG_CXX_OPT})
nrgljubljana_apply_strict_floating_point(nrgchain_core)
target_compile_definitions(nrgchain_core PRIVATE NRGCHAIN_NO_MAIN)
target_include_directories(nrgchain_core PRIVATE
                     $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/nrgchain>
                     $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/c++>
                     $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/c++/include>
                     $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/c++/include/range-v3>)

set(all_executables adapt binavg broaden bw diag h5write hilb instantiate intavg integ kk matrix mats nrgchain resample specmoments tdavg unitary)
foreach(exec ${all_executables})
  add_tool(${exec})
endforeach()

target_sources(instantiate PRIVATE $<TARGET_OBJECTS:nrgchain_core>)
target_include_directories(instantiate PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/nrgchain>)
target_link_libraries(instantiate PRIVATE nrgljubljana_c)

# Additional source file for 'matrix' tool
target_sources(matrix PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/matrix/parser.cc)
