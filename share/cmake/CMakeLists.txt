configure_file(nrgljubljana-config.cmake.in nrgljubljana-config.cmake @ONLY)
configure_file(nrgljubljana-config-version.cmake.in nrgljubljana-config-version.cmake @ONLY)
install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/nrgljubljana-config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/nrgljubljana-config-version.cmake
  DESTINATION lib/cmake/nrgljubljana
)

install(EXPORT nrgljubljana-targets NAMESPACE nrgljubljana:: DESTINATION lib/cmake/nrgljubljana)
