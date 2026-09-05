install(FILES nrginit.m dmft.m initial.m initialparse.m sneg.m floquet.m models.m operators.m basis.m hamiltonian.m wilson.m sym.m misc.m DESTINATION ${NRGLJUBLJANA_INSTALL_NRGINIT_DIR})

set(NRGINIT_PROGRAM "${CMAKE_CURRENT_SOURCE_DIR}/nrginit")
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/nrginit.local")
  set(NRGINIT_PROGRAM "${CMAKE_CURRENT_SOURCE_DIR}/nrginit.local")
endif()
file(REAL_PATH "${NRGINIT_PROGRAM}" NRGINIT_PROGRAM)

install(PROGRAMS "${NRGINIT_PROGRAM}"
        DESTINATION ${NRGLJUBLJANA_INSTALL_NRGINIT_DIR}
        RENAME nrginit)
