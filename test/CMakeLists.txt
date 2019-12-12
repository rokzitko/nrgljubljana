add_subdirectory(c++)
add_subdirectory(c++sym_basic)

if(SYM_MORE)
  add_subdirectory(c++sym_more)
endif()

if(SYM_ALL)
  add_subdirectory(c++sym_all)
endif()

add_subdirectory(tools)

if(Mathematica_FOUND)
  add_subdirectory(nrginit)
endif()
