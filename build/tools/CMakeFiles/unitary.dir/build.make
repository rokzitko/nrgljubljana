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
include tools/CMakeFiles/unitary.dir/depend.make

# Include the progress variables for this target.
include tools/CMakeFiles/unitary.dir/progress.make

# Include the compile flags for this target's objects.
include tools/CMakeFiles/unitary.dir/flags.make

tools/CMakeFiles/unitary.dir/unitary/unitary.cc.o: tools/CMakeFiles/unitary.dir/flags.make
tools/CMakeFiles/unitary.dir/unitary/unitary.cc.o: ../tools/unitary/unitary.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tools/CMakeFiles/unitary.dir/unitary/unitary.cc.o"
	cd /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/build/tools && /Library/Developer/CommandLineTools/usr/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/unitary.dir/unitary/unitary.cc.o -c /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/tools/unitary/unitary.cc

tools/CMakeFiles/unitary.dir/unitary/unitary.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/unitary.dir/unitary/unitary.cc.i"
	cd /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/build/tools && /Library/Developer/CommandLineTools/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/tools/unitary/unitary.cc > CMakeFiles/unitary.dir/unitary/unitary.cc.i

tools/CMakeFiles/unitary.dir/unitary/unitary.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/unitary.dir/unitary/unitary.cc.s"
	cd /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/build/tools && /Library/Developer/CommandLineTools/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/tools/unitary/unitary.cc -o CMakeFiles/unitary.dir/unitary/unitary.cc.s

# Object files for target unitary
unitary_OBJECTS = \
"CMakeFiles/unitary.dir/unitary/unitary.cc.o"

# External object files for target unitary
unitary_EXTERNAL_OBJECTS =

tools/unitary: tools/CMakeFiles/unitary.dir/unitary/unitary.cc.o
tools/unitary: tools/CMakeFiles/unitary.dir/build.make
tools/unitary: /usr/local/lib/libboost_serialization.dylib
tools/unitary: /usr/local/lib/libgmp.dylib
tools/unitary: /usr/local/lib/libgmpxx.dylib
tools/unitary: /usr/local/Cellar/gsl/2.6/lib/libgsl.dylib
tools/unitary: /usr/local/Cellar/gsl/2.6/lib/libgslcblas.dylib
tools/unitary: tools/CMakeFiles/unitary.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable unitary"
	cd /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/build/tools && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/unitary.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tools/CMakeFiles/unitary.dir/build: tools/unitary

.PHONY : tools/CMakeFiles/unitary.dir/build

tools/CMakeFiles/unitary.dir/clean:
	cd /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/build/tools && $(CMAKE_COMMAND) -P CMakeFiles/unitary.dir/cmake_clean.cmake
.PHONY : tools/CMakeFiles/unitary.dir/clean

tools/CMakeFiles/unitary.dir/depend:
	cd /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/tools /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/build /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/build/tools /Users/rokzitko/common/delo/ljubljana/cmake/NRGLjubljana/build/tools/CMakeFiles/unitary.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tools/CMakeFiles/unitary.dir/depend

