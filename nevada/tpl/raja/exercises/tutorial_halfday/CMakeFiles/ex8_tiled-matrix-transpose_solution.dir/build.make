# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.14

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
CMAKE_COMMAND = /usr/tce/backend/installations/linux-rhel8-x86_64/gcc-10.2.1/cmake-3.14.5-2o4l6ukye2v2vfcj3hrazdptk4vrejsn/bin/cmake

# The command to remove a file.
RM = /usr/tce/backend/installations/linux-rhel8-x86_64/gcc-10.2.1/cmake-3.14.5-2o4l6ukye2v2vfcj3hrazdptk4vrejsn/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /usr/workspace/rana2/Kripke

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /usr/workspace/rana2/Kripke/nevada

# Include any dependencies generated for this target.
include tpl/raja/exercises/tutorial_halfday/CMakeFiles/ex8_tiled-matrix-transpose_solution.dir/depend.make

# Include the progress variables for this target.
include tpl/raja/exercises/tutorial_halfday/CMakeFiles/ex8_tiled-matrix-transpose_solution.dir/progress.make

# Include the compile flags for this target's objects.
include tpl/raja/exercises/tutorial_halfday/CMakeFiles/ex8_tiled-matrix-transpose_solution.dir/flags.make

tpl/raja/exercises/tutorial_halfday/CMakeFiles/ex8_tiled-matrix-transpose_solution.dir/ex8_tiled-matrix-transpose_solution.cpp.o: tpl/raja/exercises/tutorial_halfday/CMakeFiles/ex8_tiled-matrix-transpose_solution.dir/flags.make
tpl/raja/exercises/tutorial_halfday/CMakeFiles/ex8_tiled-matrix-transpose_solution.dir/ex8_tiled-matrix-transpose_solution.cpp.o: ../tpl/raja/exercises/tutorial_halfday/ex8_tiled-matrix-transpose_solution.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/usr/workspace/rana2/Kripke/nevada/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tpl/raja/exercises/tutorial_halfday/CMakeFiles/ex8_tiled-matrix-transpose_solution.dir/ex8_tiled-matrix-transpose_solution.cpp.o"
	cd /usr/workspace/rana2/Kripke/nevada/tpl/raja/exercises/tutorial_halfday && /opt/rocm-4.2.0/llvm/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ex8_tiled-matrix-transpose_solution.dir/ex8_tiled-matrix-transpose_solution.cpp.o -c /usr/workspace/rana2/Kripke/tpl/raja/exercises/tutorial_halfday/ex8_tiled-matrix-transpose_solution.cpp

tpl/raja/exercises/tutorial_halfday/CMakeFiles/ex8_tiled-matrix-transpose_solution.dir/ex8_tiled-matrix-transpose_solution.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ex8_tiled-matrix-transpose_solution.dir/ex8_tiled-matrix-transpose_solution.cpp.i"
	cd /usr/workspace/rana2/Kripke/nevada/tpl/raja/exercises/tutorial_halfday && /opt/rocm-4.2.0/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /usr/workspace/rana2/Kripke/tpl/raja/exercises/tutorial_halfday/ex8_tiled-matrix-transpose_solution.cpp > CMakeFiles/ex8_tiled-matrix-transpose_solution.dir/ex8_tiled-matrix-transpose_solution.cpp.i

tpl/raja/exercises/tutorial_halfday/CMakeFiles/ex8_tiled-matrix-transpose_solution.dir/ex8_tiled-matrix-transpose_solution.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ex8_tiled-matrix-transpose_solution.dir/ex8_tiled-matrix-transpose_solution.cpp.s"
	cd /usr/workspace/rana2/Kripke/nevada/tpl/raja/exercises/tutorial_halfday && /opt/rocm-4.2.0/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /usr/workspace/rana2/Kripke/tpl/raja/exercises/tutorial_halfday/ex8_tiled-matrix-transpose_solution.cpp -o CMakeFiles/ex8_tiled-matrix-transpose_solution.dir/ex8_tiled-matrix-transpose_solution.cpp.s

# Object files for target ex8_tiled-matrix-transpose_solution
ex8_tiled__matrix__transpose_solution_OBJECTS = \
"CMakeFiles/ex8_tiled-matrix-transpose_solution.dir/ex8_tiled-matrix-transpose_solution.cpp.o"

# External object files for target ex8_tiled-matrix-transpose_solution
ex8_tiled__matrix__transpose_solution_EXTERNAL_OBJECTS =

bin/ex8_tiled-matrix-transpose_solution: tpl/raja/exercises/tutorial_halfday/CMakeFiles/ex8_tiled-matrix-transpose_solution.dir/ex8_tiled-matrix-transpose_solution.cpp.o
bin/ex8_tiled-matrix-transpose_solution: tpl/raja/exercises/tutorial_halfday/CMakeFiles/ex8_tiled-matrix-transpose_solution.dir/build.make
bin/ex8_tiled-matrix-transpose_solution: lib/libRAJA.a
bin/ex8_tiled-matrix-transpose_solution: tpl/raja/exercises/tutorial_halfday/CMakeFiles/ex8_tiled-matrix-transpose_solution.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/usr/workspace/rana2/Kripke/nevada/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../../../bin/ex8_tiled-matrix-transpose_solution"
	cd /usr/workspace/rana2/Kripke/nevada/tpl/raja/exercises/tutorial_halfday && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ex8_tiled-matrix-transpose_solution.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tpl/raja/exercises/tutorial_halfday/CMakeFiles/ex8_tiled-matrix-transpose_solution.dir/build: bin/ex8_tiled-matrix-transpose_solution

.PHONY : tpl/raja/exercises/tutorial_halfday/CMakeFiles/ex8_tiled-matrix-transpose_solution.dir/build

tpl/raja/exercises/tutorial_halfday/CMakeFiles/ex8_tiled-matrix-transpose_solution.dir/clean:
	cd /usr/workspace/rana2/Kripke/nevada/tpl/raja/exercises/tutorial_halfday && $(CMAKE_COMMAND) -P CMakeFiles/ex8_tiled-matrix-transpose_solution.dir/cmake_clean.cmake
.PHONY : tpl/raja/exercises/tutorial_halfday/CMakeFiles/ex8_tiled-matrix-transpose_solution.dir/clean

tpl/raja/exercises/tutorial_halfday/CMakeFiles/ex8_tiled-matrix-transpose_solution.dir/depend:
	cd /usr/workspace/rana2/Kripke/nevada && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /usr/workspace/rana2/Kripke /usr/workspace/rana2/Kripke/tpl/raja/exercises/tutorial_halfday /usr/workspace/rana2/Kripke/nevada /usr/workspace/rana2/Kripke/nevada/tpl/raja/exercises/tutorial_halfday /usr/workspace/rana2/Kripke/nevada/tpl/raja/exercises/tutorial_halfday/CMakeFiles/ex8_tiled-matrix-transpose_solution.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tpl/raja/exercises/tutorial_halfday/CMakeFiles/ex8_tiled-matrix-transpose_solution.dir/depend

