# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /mnt/c/Users/e12209452/Documents/Uni/ScientificComputing/SEJ-ODE

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/c/Users/e12209452/Documents/Uni/ScientificComputing/SEJ-ODE/build

# Include any dependencies generated for this target.
include CMakeFiles/test_ode.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/test_ode.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/test_ode.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/test_ode.dir/flags.make

CMakeFiles/test_ode.dir/demos/test_ode.cc.o: CMakeFiles/test_ode.dir/flags.make
CMakeFiles/test_ode.dir/demos/test_ode.cc.o: ../demos/test_ode.cc
CMakeFiles/test_ode.dir/demos/test_ode.cc.o: CMakeFiles/test_ode.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/Users/e12209452/Documents/Uni/ScientificComputing/SEJ-ODE/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/test_ode.dir/demos/test_ode.cc.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/test_ode.dir/demos/test_ode.cc.o -MF CMakeFiles/test_ode.dir/demos/test_ode.cc.o.d -o CMakeFiles/test_ode.dir/demos/test_ode.cc.o -c /mnt/c/Users/e12209452/Documents/Uni/ScientificComputing/SEJ-ODE/demos/test_ode.cc

CMakeFiles/test_ode.dir/demos/test_ode.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_ode.dir/demos/test_ode.cc.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/Users/e12209452/Documents/Uni/ScientificComputing/SEJ-ODE/demos/test_ode.cc > CMakeFiles/test_ode.dir/demos/test_ode.cc.i

CMakeFiles/test_ode.dir/demos/test_ode.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_ode.dir/demos/test_ode.cc.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/Users/e12209452/Documents/Uni/ScientificComputing/SEJ-ODE/demos/test_ode.cc -o CMakeFiles/test_ode.dir/demos/test_ode.cc.s

# Object files for target test_ode
test_ode_OBJECTS = \
"CMakeFiles/test_ode.dir/demos/test_ode.cc.o"

# External object files for target test_ode
test_ode_EXTERNAL_OBJECTS =

test_ode: CMakeFiles/test_ode.dir/demos/test_ode.cc.o
test_ode: CMakeFiles/test_ode.dir/build.make
test_ode: /usr/lib/x86_64-linux-gnu/libopenblas.so
test_ode: CMakeFiles/test_ode.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/c/Users/e12209452/Documents/Uni/ScientificComputing/SEJ-ODE/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable test_ode"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_ode.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/test_ode.dir/build: test_ode
.PHONY : CMakeFiles/test_ode.dir/build

CMakeFiles/test_ode.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/test_ode.dir/cmake_clean.cmake
.PHONY : CMakeFiles/test_ode.dir/clean

CMakeFiles/test_ode.dir/depend:
	cd /mnt/c/Users/e12209452/Documents/Uni/ScientificComputing/SEJ-ODE/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/c/Users/e12209452/Documents/Uni/ScientificComputing/SEJ-ODE /mnt/c/Users/e12209452/Documents/Uni/ScientificComputing/SEJ-ODE /mnt/c/Users/e12209452/Documents/Uni/ScientificComputing/SEJ-ODE/build /mnt/c/Users/e12209452/Documents/Uni/ScientificComputing/SEJ-ODE/build /mnt/c/Users/e12209452/Documents/Uni/ScientificComputing/SEJ-ODE/build/CMakeFiles/test_ode.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/test_ode.dir/depend
