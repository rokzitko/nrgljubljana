add_library(nrgljubljana_workflow SHARED workflow.cc discretize.cc chain.cc preparation.cc ../instantiate/instantiate.cc
  $<TARGET_OBJECTS:nrgchain_core>)
add_library(nrgljubljana::workflow ALIAS nrgljubljana_workflow)
target_compile_features(nrgljubljana_workflow PUBLIC cxx_std_20)
target_compile_definitions(nrgljubljana_workflow PRIVATE NRG_INSTANTIATE_NO_MAIN)
target_compile_options(nrgljubljana_workflow PRIVATE ${NRG_CXX_OPT})
nrgljubljana_apply_strict_floating_point(nrgljubljana_workflow)
target_include_directories(nrgljubljana_workflow
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE ${PROJECT_BINARY_DIR}/tools/generated ${CMAKE_CURRENT_SOURCE_DIR}/../nrgchain)
target_link_libraries(nrgljubljana_workflow PRIVATE $<BUILD_INTERFACE:nrgljubljana_c> project_warnings gsl Threads::Threads)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  # The matching private core must precede an older same-SONAME LD_LIBRARY_PATH entry.
  # Put this after link options: MPI may contribute --enable-new-dtags there.
  target_link_libraries(nrgljubljana_workflow PRIVATE "-Wl,--disable-new-dtags")
endif()
set_target_properties(nrgljubljana_workflow PROPERTIES
  EXPORT_NAME workflow VERSION ${PROJECT_VERSION}
  WORKFLOW_API_VERSION 3 EXPORT_PROPERTIES WORKFLOW_API_VERSION
  BUILD_RPATH "${PROJECT_BINARY_DIR}/c++"
  INSTALL_RPATH "$ORIGIN/nrgljubljana/workflow")

include(CMakePackageConfigHelpers)
configure_package_config_file(nrgljubljana_workflowConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/nrgljubljana_workflowConfig.configured.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/nrgljubljana_workflow)
file(GENERATE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/nrgljubljana_workflowConfig.cmake
  INPUT ${CMAKE_CURRENT_BINARY_DIR}/nrgljubljana_workflowConfig.configured.cmake)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/nrgljubljana_workflowConfigVersion.cmake
  VERSION ${PROJECT_VERSION} COMPATIBILITY SameMajorVersion)
install(TARGETS nrgljubljana_workflow EXPORT nrgljubljana-workflow-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT workflow
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT workflow)
# Keep the matching core private: no installed internal headers or old core package lookup.
install(FILES $<TARGET_FILE:nrgljubljana_c>
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/nrgljubljana/workflow COMPONENT workflow)
install(FILES include/nrgljubljana/solver_workflow.hpp include/nrgljubljana/solver_preparation.hpp include/nrgljubljana/solver_chain.hpp
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nrgljubljana COMPONENT workflow)
install(EXPORT nrgljubljana-workflow-targets NAMESPACE nrgljubljana::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/nrgljubljana_workflow COMPONENT workflow)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/nrgljubljana_workflowConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/nrgljubljana_workflowConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/nrgljubljana_workflow COMPONENT workflow)

if(NRGLJUBLJANA_BUILD_TESTS_ENABLED)
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
  enable_testing()
  # Exercise production routing/discretization without linking or running the solver.
  # The driver supplies small initializer/core stand-ins that write to all four streams.
  add_executable(workflow_diagnostics tests/diagnostics.cc workflow.cc discretize.cc)
  target_compile_features(workflow_diagnostics PRIVATE cxx_std_20)
  target_include_directories(workflow_diagnostics PRIVATE include ${PROJECT_SOURCE_DIR}/c++)
  target_link_libraries(workflow_diagnostics PRIVATE gmp gsl Boost::mpi fmt::fmt-header-only Threads::Threads project_warnings)
  target_link_libraries(workflow_diagnostics PRIVATE $<$<BOOL:${ASAN}>:asan> $<$<BOOL:${UBSAN}>:ubsan>)
  nrgljubljana_apply_strict_floating_point(workflow_diagnostics)
  add_test(NAME solver_workflow_diagnostics COMMAND workflow_diagnostics ${CMAKE_CURRENT_BINARY_DIR})
  add_test(NAME solver_workflow_diagnostics_unsynced COMMAND workflow_diagnostics ${CMAKE_CURRENT_BINARY_DIR} unsynced)
  set_tests_properties(solver_workflow_diagnostics solver_workflow_diagnostics_unsynced PROPERTIES TIMEOUT 30)

  add_executable(workflow_preparation tests/preparation.cpp)
  target_compile_features(workflow_preparation PRIVATE cxx_std_20)
  target_include_directories(workflow_preparation PRIVATE include)
  target_link_libraries(workflow_preparation PRIVATE nrgljubljana::workflow project_warnings
    $<$<BOOL:${ASAN}>:asan> $<$<BOOL:${UBSAN}>:ubsan>)
  nrgljubljana_apply_strict_floating_point(workflow_preparation)
  add_test(NAME solver_workflow_preparation COMMAND workflow_preparation ${PROJECT_SOURCE_DIR}/test/mixchain+nrgrun/templates/QS)
  set_tests_properties(solver_workflow_preparation PROPERTIES TIMEOUT 30)

  add_executable(workflow_contract tests/contract.cc)
  target_link_libraries(workflow_contract PRIVATE nrgljubljana::workflow gmp gsl project_warnings Boost::headers fmt::fmt-header-only)
  target_link_libraries(workflow_contract PRIVATE $<$<BOOL:${ASAN}>:asan> $<$<BOOL:${UBSAN}>:ubsan>)
  target_include_directories(workflow_contract PRIVATE ${PROJECT_SOURCE_DIR}/c++)
  add_test(NAME solver_workflow_contract COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/contract.py
    --driver $<TARGET_FILE:workflow_contract> --tools ${PROJECT_BINARY_DIR}/tools
    --template ${PROJECT_SOURCE_DIR}/test/mixchain+nrgrun/templates/QS
    --temporary-parent ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(solver_workflow_contract PROPERTIES TIMEOUT 180 ENVIRONMENT "OMP_NUM_THREADS=1;MKL_NUM_THREADS=1")
  add_test(NAME solver_workflow_comparator COMMAND ${Python3_EXECUTABLE} -B ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_contract.py)

  # discretize_hybridization against the chain.dat files of the mixchain executable
  add_executable(workflow_chain tests/chain.cc)
  target_compile_features(workflow_chain PRIVATE cxx_std_20)
  target_link_libraries(workflow_chain PRIVATE nrgljubljana::workflow project_warnings
    $<$<BOOL:${ASAN}>:asan> $<$<BOOL:${UBSAN}>:ubsan>)
  add_test(NAME solver_workflow_chain COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests/chain_test ${PROJECT_SOURCE_DIR}
    ${PROJECT_BINARY_DIR}/tools ${CMAKE_CURRENT_BINARY_DIR}/chain $<TARGET_FILE:workflow_chain>)
  set_tests_properties(solver_workflow_chain PROPERTIES TIMEOUT 300 ENVIRONMENT "OMP_NUM_THREADS=1;MKL_NUM_THREADS=1")

  # solve from a host that initializes MPI: serially, and one calculation per rank
  add_executable(workflow_solve_mpi tests/solve_mpi.cc)
  target_link_libraries(workflow_solve_mpi PRIVATE nrgljubljana::workflow MPI::MPI_CXX project_warnings
    $<$<BOOL:${ASAN}>:asan> $<$<BOOL:${UBSAN}>:ubsan>)
  add_test(NAME solver_workflow_solve COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests/solve_test ${PROJECT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}/solve $<TARGET_FILE:workflow_solve_mpi>)
  set_tests_properties(solver_workflow_solve PROPERTIES TIMEOUT 120 ENVIRONMENT "OMP_NUM_THREADS=1;MKL_NUM_THREADS=1")
  if(NRG_MPI)
    add_test(NAME solver_workflow_solve_np2 COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests/solve_test ${PROJECT_SOURCE_DIR}
      ${CMAKE_CURRENT_BINARY_DIR}/solve_np2 ${MPIEXEC_EXECUTABLE} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 2
      $<TARGET_FILE:workflow_solve_mpi> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(solver_workflow_solve_np2 PROPERTIES TIMEOUT 120 PROCESSORS 2
      ENVIRONMENT "OMP_NUM_THREADS=1;MKL_NUM_THREADS=1")
  endif()
endif()
