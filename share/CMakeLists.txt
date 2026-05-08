add_subdirectory(cmake)

configure_file(nrgljubljana.modulefile.in nrgljubljana.modulefile @ONLY)
configure_file(nrgljubljanavars.sh.in nrgljubljanavars.sh @ONLY)

install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/nrgljubljana.modulefile
    ${CMAKE_CURRENT_BINARY_DIR}/nrgljubljanavars.sh
  DESTINATION ${CMAKE_INSTALL_DATADIR}
)

message(STATUS "***************************************************************")
message(STATUS "* Use:                                                         ")
message(STATUS "*                                                              ")
message(STATUS "*   source ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/nrgljubljanavars.sh      ")
message(STATUS "*                                                              ")
message(STATUS "* to set up the environment variables                          ")
message(STATUS "***************************************************************")
