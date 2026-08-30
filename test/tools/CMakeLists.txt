add_subdirectory(adapt)
add_subdirectory(binavg)
add_subdirectory(broaden)
add_subdirectory(bw)
add_subdirectory(diag)
add_subdirectory(h5write)
add_subdirectory(hilb)
add_subdirectory(intavg)
add_subdirectory(integ)
add_subdirectory(instantiate)
add_subdirectory(kk)
add_subdirectory(matrix)
add_subdirectory(mats)
add_subdirectory(nrgchain)
add_subdirectory(resample)
add_subdirectory(specmoments)
add_subdirectory(tdavg)
add_subdirectory(unitary)

set(NRG_VERSIONED_EXECUTABLES ${NRG_TOOL_EXECUTABLES} nrg)
foreach(tool IN LISTS NRG_VERSIONED_EXECUTABLES)
  add_test(
    NAME tool_version_${tool}
    COMMAND ${CMAKE_COMMAND}
      -DTOOL=$<TARGET_FILE:${tool}>
      -DTOOL_NAME=${tool}
      -DEXPECTED_VERSION=${PROJECT_VERSION}
      -P ${CMAKE_CURRENT_SOURCE_DIR}/check-version.cmake
  )
endforeach()
