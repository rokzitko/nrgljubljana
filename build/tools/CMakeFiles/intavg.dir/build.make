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
include tools/CMakeFiles/intavg.dir/depend.make

# Include the progress variables for this target.
include tools/CMakeFiles/intavg.dir/progress.make

# Include the compile flags for this target's objects.
include tools/CMakeFiles/intavg.dir/flags.make

tools/CMakeFiles/intavg.dir/intavg/intavg.cc.o: tools/CMakeFiles/intavg.dir/flags.make
tools/CMakeFiles/intavg.dir/intavg/intavg.cc.o: ../tools/intavg/intavg.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tools/CMakeFiles/intavg.dir/intavg/intavg.cc.o"
	cd /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/build/tools && /Library/Developer/CommandLineTools/usr/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/intavg.dir/intavg/intavg.cc.o -c /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/tools/intavg/intavg.cc

tools/CMakeFiles/intavg.dir/intavg/intavg.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/intavg.dir/intavg/intavg.cc.i"
	cd /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/build/tools && /Library/Developer/CommandLineTools/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/tools/intavg/intavg.cc > CMakeFiles/intavg.dir/intavg/intavg.cc.i

tools/CMakeFiles/intavg.dir/intavg/intavg.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/intavg.dir/intavg/intavg.cc.s"
	cd /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/build/tools && /Library/Developer/CommandLineTools/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/tools/intavg/intavg.cc -o CMakeFiles/intavg.dir/intavg/intavg.cc.s

# Object files for target intavg
intavg_OBJECTS = \
"CMakeFiles/intavg.dir/intavg/intavg.cc.o"

# External object files for target intavg
intavg_EXTERNAL_OBJECTS =

tools/intavg: tools/CMakeFiles/intavg.dir/intavg/intavg.cc.o
tools/intavg: tools/CMakeFiles/intavg.dir/build.make
tools/intavg: /usr/local/lib/libboost_serialization.dylib
tools/intavg: /usr/local/lib/libgmp.dylib
tools/intavg: /usr/local/lib/libgmpxx.dylib
tools/intavg: /usr/local/Cellar/gsl/2.6/lib/libgsl.dylib
tools/intavg: /usr/local/Cellar/gsl/2.6/lib/libgslcblas.dylib
tools/intavg: tools/CMakeFiles/intavg.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable intavg"
	cd /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/build/tools && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/intavg.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tools/CMakeFiles/intavg.dir/build: tools/intavg

.PHONY : tools/CMakeFiles/intavg.dir/build

tools/CMakeFiles/intavg.dir/clean:
	cd /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/build/tools && $(CMAKE_COMMAND) -P CMakeFiles/intavg.dir/cmake_clean.cmake
.PHONY : tools/CMakeFiles/intavg.dir/clean

tools/CMakeFiles/intavg.dir/depend:
	cd /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/tools /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/build /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/build/tools /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/build/tools/CMakeFiles/intavg.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tools/CMakeFiles/intavg.dir/depend

