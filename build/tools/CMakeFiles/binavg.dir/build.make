# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.15.4/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.15.4/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/build

# Include any dependencies generated for this target.
include tools/CMakeFiles/binavg.dir/depend.make

# Include the progress variables for this target.
include tools/CMakeFiles/binavg.dir/progress.make

# Include the compile flags for this target's objects.
include tools/CMakeFiles/binavg.dir/flags.make

tools/CMakeFiles/binavg.dir/binavg/binavg.cc.o: tools/CMakeFiles/binavg.dir/flags.make
tools/CMakeFiles/binavg.dir/binavg/binavg.cc.o: ../tools/binavg/binavg.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tools/CMakeFiles/binavg.dir/binavg/binavg.cc.o"
	cd /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/build/tools && /Library/Developer/CommandLineTools/usr/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/binavg.dir/binavg/binavg.cc.o -c /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/tools/binavg/binavg.cc

tools/CMakeFiles/binavg.dir/binavg/binavg.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/binavg.dir/binavg/binavg.cc.i"
	cd /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/build/tools && /Library/Developer/CommandLineTools/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/tools/binavg/binavg.cc > CMakeFiles/binavg.dir/binavg/binavg.cc.i

tools/CMakeFiles/binavg.dir/binavg/binavg.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/binavg.dir/binavg/binavg.cc.s"
	cd /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/build/tools && /Library/Developer/CommandLineTools/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/tools/binavg/binavg.cc -o CMakeFiles/binavg.dir/binavg/binavg.cc.s

# Object files for target binavg
binavg_OBJECTS = \
"CMakeFiles/binavg.dir/binavg/binavg.cc.o"

# External object files for target binavg
binavg_EXTERNAL_OBJECTS =

tools/binavg: tools/CMakeFiles/binavg.dir/binavg/binavg.cc.o
tools/binavg: tools/CMakeFiles/binavg.dir/build.make
tools/binavg: /usr/local/lib/libboost_serialization.dylib
tools/binavg: /usr/local/lib/libgmp.dylib
tools/binavg: /usr/local/lib/libgmpxx.dylib
tools/binavg: /usr/local/Cellar/gsl/2.6/lib/libgsl.dylib
tools/binavg: /usr/local/Cellar/gsl/2.6/lib/libgslcblas.dylib
tools/binavg: tools/CMakeFiles/binavg.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable binavg"
	cd /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/build/tools && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/binavg.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tools/CMakeFiles/binavg.dir/build: tools/binavg

.PHONY : tools/CMakeFiles/binavg.dir/build

tools/CMakeFiles/binavg.dir/clean:
	cd /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/build/tools && $(CMAKE_COMMAND) -P CMakeFiles/binavg.dir/cmake_clean.cmake
.PHONY : tools/CMakeFiles/binavg.dir/clean

tools/CMakeFiles/binavg.dir/depend:
	cd /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/tools /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/build /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/build/tools /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/build/tools/CMakeFiles/binavg.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tools/CMakeFiles/binavg.dir/depend

