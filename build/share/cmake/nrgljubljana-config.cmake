# This file allows other CMake Projects to find us
# We provide general project information
# and reestablish the exported CMake Targets

# Multiple inclusion guard
if(NOT NRGLJUBLJANA_FOUND)
set(NRGLJUBLJANA_FOUND True)

# version
set(NRGLJUBLJANA_VERSION 2019.11)
set(NRGLJUBLJANA_GIT_HASH 7fe29c4aba3f7c23af5e524404b1d1552a7a134c)

# Root of the installation
set(NRGLJUBLJANA_ROOT  /Users/rokzitko/nrgljubljana/2019.11)

# Include the exported targets of this project
include(/Users/rokzitko/nrgljubljana/2019.11/lib/cmake/nrgljubljana/nrgljubljana-targets.cmake)

message(STATUS "Found nrgljubljana-config.cmake with version 2019.11, hash = 7fe29c4aba3f7c23af5e524404b1d1552a7a134c")

# Was the Project built with Documentation?
set(NRGLJUBLJANA_WITH_DOCUMENTATION )

endif()
