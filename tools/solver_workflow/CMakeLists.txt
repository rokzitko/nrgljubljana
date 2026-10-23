add_library(nrgljubljana_workflow SHARED workflow.cc discretize.cc ../instantiate/instantiate.cc
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
install(FILES include/nrgljubljana/solver_workflow.hpp
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nrgljubljana COMPONENT workflow)
install(EXPORT nrgljubljana-workflow-targets NAMESPACE nrgljubljana::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/nrgljubljana_workflow COMPONENT workflow)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/nrgljubljana_workflowConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/nrgljubljana_workflowConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/nrgljubljana_workflow COMPONENT workflow)

if(NRGLJUBLJANA_BUILD_TESTS_ENABLED)
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
  enable_testing()
  add_executable(workflow_contract tests/contract.cc)
  target_link_libraries(workflow_contract PRIVATE nrgljubljana::workflow gmp gsl project_warnings Boost::headers fmt::fmt-header-only)
  target_include_directories(workflow_contract PRIVATE ${PROJECT_SOURCE_DIR}/c++)
  add_test(NAME solver_workflow_contract COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/contract.py
    --driver $<TARGET_FILE:workflow_contract> --tools ${PROJECT_BINARY_DIR}/tools
    --template ${PROJECT_SOURCE_DIR}/test/mixchain+nrgrun/templates/QS
    --temporary-parent ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(solver_workflow_contract PROPERTIES TIMEOUT 180 ENVIRONMENT "OMP_NUM_THREADS=1;MKL_NUM_THREADS=1")
endif()
