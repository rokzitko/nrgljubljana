add_subdirectory(c++)

add_subdirectory(tools)

if(Mathematica_FOUND)
  add_subdirectory(nrginit)
endif()
