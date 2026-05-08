include(CMakePackageConfigHelpers)

set(NRGLJUBLJANA_CMAKE_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/nrgljubljana)

configure_package_config_file(
  nrgljubljana-config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/nrgljubljana-config.cmake
  INSTALL_DESTINATION ${NRGLJUBLJANA_CMAKE_INSTALL_DIR}
)
configure_file(nrgljubljana-config-version.cmake.in nrgljubljana-config-version.cmake @ONLY)
install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/nrgljubljana-config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/nrgljubljana-config-version.cmake
  DESTINATION ${NRGLJUBLJANA_CMAKE_INSTALL_DIR}
)

install(EXPORT nrgljubljana-targets NAMESPACE nrgljubljana:: DESTINATION ${NRGLJUBLJANA_CMAKE_INSTALL_DIR})
