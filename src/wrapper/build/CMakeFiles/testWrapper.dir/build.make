# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.10.2/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.10.2/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/derek/github/freestream-milne/src/wrapper

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/derek/github/freestream-milne/src/wrapper/build

# Include any dependencies generated for this target.
include CMakeFiles/testWrapper.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/testWrapper.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/testWrapper.dir/flags.make

CMakeFiles/testWrapper.dir/testWrapper.cpp.o: CMakeFiles/testWrapper.dir/flags.make
CMakeFiles/testWrapper.dir/testWrapper.cpp.o: ../testWrapper.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/derek/github/freestream-milne/src/wrapper/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/testWrapper.dir/testWrapper.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/testWrapper.dir/testWrapper.cpp.o -c /Users/derek/github/freestream-milne/src/wrapper/testWrapper.cpp

CMakeFiles/testWrapper.dir/testWrapper.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/testWrapper.dir/testWrapper.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/derek/github/freestream-milne/src/wrapper/testWrapper.cpp > CMakeFiles/testWrapper.dir/testWrapper.cpp.i

CMakeFiles/testWrapper.dir/testWrapper.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/testWrapper.dir/testWrapper.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/derek/github/freestream-milne/src/wrapper/testWrapper.cpp -o CMakeFiles/testWrapper.dir/testWrapper.cpp.s

CMakeFiles/testWrapper.dir/testWrapper.cpp.o.requires:

.PHONY : CMakeFiles/testWrapper.dir/testWrapper.cpp.o.requires

CMakeFiles/testWrapper.dir/testWrapper.cpp.o.provides: CMakeFiles/testWrapper.dir/testWrapper.cpp.o.requires
	$(MAKE) -f CMakeFiles/testWrapper.dir/build.make CMakeFiles/testWrapper.dir/testWrapper.cpp.o.provides.build
.PHONY : CMakeFiles/testWrapper.dir/testWrapper.cpp.o.provides

CMakeFiles/testWrapper.dir/testWrapper.cpp.o.provides.build: CMakeFiles/testWrapper.dir/testWrapper.cpp.o


# Object files for target testWrapper
testWrapper_OBJECTS = \
"CMakeFiles/testWrapper.dir/testWrapper.cpp.o"

# External object files for target testWrapper
testWrapper_EXTERNAL_OBJECTS =

testWrapper: CMakeFiles/testWrapper.dir/testWrapper.cpp.o
testWrapper: CMakeFiles/testWrapper.dir/build.make
testWrapper: /usr/local/Cellar/gsl/2.4/lib/libgsl.dylib
testWrapper: /usr/local/Cellar/gsl/2.4/lib/libgslcblas.dylib
testWrapper: CMakeFiles/testWrapper.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/derek/github/freestream-milne/src/wrapper/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable testWrapper"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/testWrapper.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/testWrapper.dir/build: testWrapper

.PHONY : CMakeFiles/testWrapper.dir/build

CMakeFiles/testWrapper.dir/requires: CMakeFiles/testWrapper.dir/testWrapper.cpp.o.requires

.PHONY : CMakeFiles/testWrapper.dir/requires

CMakeFiles/testWrapper.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/testWrapper.dir/cmake_clean.cmake
.PHONY : CMakeFiles/testWrapper.dir/clean

CMakeFiles/testWrapper.dir/depend:
	cd /Users/derek/github/freestream-milne/src/wrapper/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/derek/github/freestream-milne/src/wrapper /Users/derek/github/freestream-milne/src/wrapper /Users/derek/github/freestream-milne/src/wrapper/build /Users/derek/github/freestream-milne/src/wrapper/build /Users/derek/github/freestream-milne/src/wrapper/build/CMakeFiles/testWrapper.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/testWrapper.dir/depend
