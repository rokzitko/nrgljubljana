add_subdirectory(unit)

add_subdirectory(c++)
add_subdirectory(c++sym_basic)

if(SYM_MORE)
  add_subdirectory(c++sym_more)
endif()

if(SYM_ALL)
  add_subdirectory(c++sym_all)
  add_subdirectory(complex)
endif()

if(TEST_LONG AND SYM_ALL)
  add_subdirectory(test_long)
endif()

add_subdirectory(tools)

if(Mathematica_FOUND)
  add_subdirectory(nrginit)
endif()

