# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.28

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
CMAKE_SOURCE_DIR = /home/cvdg/Desktop/TSP_test/blossom_inequality_sep

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/build

# Include any dependencies generated for this target.
include CMakeFiles/blossom_separation.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/blossom_separation.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/blossom_separation.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/blossom_separation.dir/flags.make

CMakeFiles/blossom_separation.dir/blossom_detector.c.o: CMakeFiles/blossom_separation.dir/flags.make
CMakeFiles/blossom_separation.dir/blossom_detector.c.o: /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/blossom_detector.c
CMakeFiles/blossom_separation.dir/blossom_detector.c.o: CMakeFiles/blossom_separation.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/cvdg/Desktop/TSP_test/blossom_inequality_sep/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/blossom_separation.dir/blossom_detector.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/blossom_separation.dir/blossom_detector.c.o -MF CMakeFiles/blossom_separation.dir/blossom_detector.c.o.d -o CMakeFiles/blossom_separation.dir/blossom_detector.c.o -c /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/blossom_detector.c

CMakeFiles/blossom_separation.dir/blossom_detector.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/blossom_separation.dir/blossom_detector.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/blossom_detector.c > CMakeFiles/blossom_separation.dir/blossom_detector.c.i

CMakeFiles/blossom_separation.dir/blossom_detector.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/blossom_separation.dir/blossom_detector.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/blossom_detector.c -o CMakeFiles/blossom_separation.dir/blossom_detector.c.s

CMakeFiles/blossom_separation.dir/src/blossom.c.o: CMakeFiles/blossom_separation.dir/flags.make
CMakeFiles/blossom_separation.dir/src/blossom.c.o: /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/blossom.c
CMakeFiles/blossom_separation.dir/src/blossom.c.o: CMakeFiles/blossom_separation.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/cvdg/Desktop/TSP_test/blossom_inequality_sep/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/blossom_separation.dir/src/blossom.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/blossom_separation.dir/src/blossom.c.o -MF CMakeFiles/blossom_separation.dir/src/blossom.c.o.d -o CMakeFiles/blossom_separation.dir/src/blossom.c.o -c /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/blossom.c

CMakeFiles/blossom_separation.dir/src/blossom.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/blossom_separation.dir/src/blossom.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/blossom.c > CMakeFiles/blossom_separation.dir/src/blossom.c.i

CMakeFiles/blossom_separation.dir/src/blossom.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/blossom_separation.dir/src/blossom.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/blossom.c -o CMakeFiles/blossom_separation.dir/src/blossom.c.s

CMakeFiles/blossom_separation.dir/src/cliqwork.c.o: CMakeFiles/blossom_separation.dir/flags.make
CMakeFiles/blossom_separation.dir/src/cliqwork.c.o: /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/cliqwork.c
CMakeFiles/blossom_separation.dir/src/cliqwork.c.o: CMakeFiles/blossom_separation.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/cvdg/Desktop/TSP_test/blossom_inequality_sep/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object CMakeFiles/blossom_separation.dir/src/cliqwork.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/blossom_separation.dir/src/cliqwork.c.o -MF CMakeFiles/blossom_separation.dir/src/cliqwork.c.o.d -o CMakeFiles/blossom_separation.dir/src/cliqwork.c.o -c /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/cliqwork.c

CMakeFiles/blossom_separation.dir/src/cliqwork.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/blossom_separation.dir/src/cliqwork.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/cliqwork.c > CMakeFiles/blossom_separation.dir/src/cliqwork.c.i

CMakeFiles/blossom_separation.dir/src/cliqwork.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/blossom_separation.dir/src/cliqwork.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/cliqwork.c -o CMakeFiles/blossom_separation.dir/src/cliqwork.c.s

CMakeFiles/blossom_separation.dir/src/skeleton.c.o: CMakeFiles/blossom_separation.dir/flags.make
CMakeFiles/blossom_separation.dir/src/skeleton.c.o: /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/skeleton.c
CMakeFiles/blossom_separation.dir/src/skeleton.c.o: CMakeFiles/blossom_separation.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/cvdg/Desktop/TSP_test/blossom_inequality_sep/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object CMakeFiles/blossom_separation.dir/src/skeleton.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/blossom_separation.dir/src/skeleton.c.o -MF CMakeFiles/blossom_separation.dir/src/skeleton.c.o.d -o CMakeFiles/blossom_separation.dir/src/skeleton.c.o -c /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/skeleton.c

CMakeFiles/blossom_separation.dir/src/skeleton.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/blossom_separation.dir/src/skeleton.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/skeleton.c > CMakeFiles/blossom_separation.dir/src/skeleton.c.i

CMakeFiles/blossom_separation.dir/src/skeleton.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/blossom_separation.dir/src/skeleton.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/skeleton.c -o CMakeFiles/blossom_separation.dir/src/skeleton.c.s

CMakeFiles/blossom_separation.dir/src/cutpool.c.o: CMakeFiles/blossom_separation.dir/flags.make
CMakeFiles/blossom_separation.dir/src/cutpool.c.o: /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/cutpool.c
CMakeFiles/blossom_separation.dir/src/cutpool.c.o: CMakeFiles/blossom_separation.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/cvdg/Desktop/TSP_test/blossom_inequality_sep/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building C object CMakeFiles/blossom_separation.dir/src/cutpool.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/blossom_separation.dir/src/cutpool.c.o -MF CMakeFiles/blossom_separation.dir/src/cutpool.c.o.d -o CMakeFiles/blossom_separation.dir/src/cutpool.c.o -c /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/cutpool.c

CMakeFiles/blossom_separation.dir/src/cutpool.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/blossom_separation.dir/src/cutpool.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/cutpool.c > CMakeFiles/blossom_separation.dir/src/cutpool.c.i

CMakeFiles/blossom_separation.dir/src/cutpool.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/blossom_separation.dir/src/cutpool.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/cutpool.c -o CMakeFiles/blossom_separation.dir/src/cutpool.c.s

CMakeFiles/blossom_separation.dir/src/allocrus.c.o: CMakeFiles/blossom_separation.dir/flags.make
CMakeFiles/blossom_separation.dir/src/allocrus.c.o: /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/allocrus.c
CMakeFiles/blossom_separation.dir/src/allocrus.c.o: CMakeFiles/blossom_separation.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/cvdg/Desktop/TSP_test/blossom_inequality_sep/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building C object CMakeFiles/blossom_separation.dir/src/allocrus.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/blossom_separation.dir/src/allocrus.c.o -MF CMakeFiles/blossom_separation.dir/src/allocrus.c.o.d -o CMakeFiles/blossom_separation.dir/src/allocrus.c.o -c /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/allocrus.c

CMakeFiles/blossom_separation.dir/src/allocrus.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/blossom_separation.dir/src/allocrus.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/allocrus.c > CMakeFiles/blossom_separation.dir/src/allocrus.c.i

CMakeFiles/blossom_separation.dir/src/allocrus.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/blossom_separation.dir/src/allocrus.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/allocrus.c -o CMakeFiles/blossom_separation.dir/src/allocrus.c.s

CMakeFiles/blossom_separation.dir/src/urandom.c.o: CMakeFiles/blossom_separation.dir/flags.make
CMakeFiles/blossom_separation.dir/src/urandom.c.o: /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/urandom.c
CMakeFiles/blossom_separation.dir/src/urandom.c.o: CMakeFiles/blossom_separation.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/cvdg/Desktop/TSP_test/blossom_inequality_sep/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building C object CMakeFiles/blossom_separation.dir/src/urandom.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/blossom_separation.dir/src/urandom.c.o -MF CMakeFiles/blossom_separation.dir/src/urandom.c.o.d -o CMakeFiles/blossom_separation.dir/src/urandom.c.o -c /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/urandom.c

CMakeFiles/blossom_separation.dir/src/urandom.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/blossom_separation.dir/src/urandom.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/urandom.c > CMakeFiles/blossom_separation.dir/src/urandom.c.i

CMakeFiles/blossom_separation.dir/src/urandom.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/blossom_separation.dir/src/urandom.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/urandom.c -o CMakeFiles/blossom_separation.dir/src/urandom.c.s

CMakeFiles/blossom_separation.dir/src/gomoryhu.c.o: CMakeFiles/blossom_separation.dir/flags.make
CMakeFiles/blossom_separation.dir/src/gomoryhu.c.o: /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/gomoryhu.c
CMakeFiles/blossom_separation.dir/src/gomoryhu.c.o: CMakeFiles/blossom_separation.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/cvdg/Desktop/TSP_test/blossom_inequality_sep/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building C object CMakeFiles/blossom_separation.dir/src/gomoryhu.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/blossom_separation.dir/src/gomoryhu.c.o -MF CMakeFiles/blossom_separation.dir/src/gomoryhu.c.o.d -o CMakeFiles/blossom_separation.dir/src/gomoryhu.c.o -c /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/gomoryhu.c

CMakeFiles/blossom_separation.dir/src/gomoryhu.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/blossom_separation.dir/src/gomoryhu.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/gomoryhu.c > CMakeFiles/blossom_separation.dir/src/gomoryhu.c.i

CMakeFiles/blossom_separation.dir/src/gomoryhu.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/blossom_separation.dir/src/gomoryhu.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/gomoryhu.c -o CMakeFiles/blossom_separation.dir/src/gomoryhu.c.s

CMakeFiles/blossom_separation.dir/src/sortrus.c.o: CMakeFiles/blossom_separation.dir/flags.make
CMakeFiles/blossom_separation.dir/src/sortrus.c.o: /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/sortrus.c
CMakeFiles/blossom_separation.dir/src/sortrus.c.o: CMakeFiles/blossom_separation.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/cvdg/Desktop/TSP_test/blossom_inequality_sep/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building C object CMakeFiles/blossom_separation.dir/src/sortrus.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/blossom_separation.dir/src/sortrus.c.o -MF CMakeFiles/blossom_separation.dir/src/sortrus.c.o.d -o CMakeFiles/blossom_separation.dir/src/sortrus.c.o -c /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/sortrus.c

CMakeFiles/blossom_separation.dir/src/sortrus.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/blossom_separation.dir/src/sortrus.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/sortrus.c > CMakeFiles/blossom_separation.dir/src/sortrus.c.i

CMakeFiles/blossom_separation.dir/src/sortrus.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/blossom_separation.dir/src/sortrus.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/sortrus.c -o CMakeFiles/blossom_separation.dir/src/sortrus.c.s

CMakeFiles/blossom_separation.dir/src/safe_io.c.o: CMakeFiles/blossom_separation.dir/flags.make
CMakeFiles/blossom_separation.dir/src/safe_io.c.o: /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/safe_io.c
CMakeFiles/blossom_separation.dir/src/safe_io.c.o: CMakeFiles/blossom_separation.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/cvdg/Desktop/TSP_test/blossom_inequality_sep/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building C object CMakeFiles/blossom_separation.dir/src/safe_io.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/blossom_separation.dir/src/safe_io.c.o -MF CMakeFiles/blossom_separation.dir/src/safe_io.c.o.d -o CMakeFiles/blossom_separation.dir/src/safe_io.c.o -c /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/safe_io.c

CMakeFiles/blossom_separation.dir/src/safe_io.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/blossom_separation.dir/src/safe_io.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/safe_io.c > CMakeFiles/blossom_separation.dir/src/safe_io.c.i

CMakeFiles/blossom_separation.dir/src/safe_io.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/blossom_separation.dir/src/safe_io.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/safe_io.c -o CMakeFiles/blossom_separation.dir/src/safe_io.c.s

CMakeFiles/blossom_separation.dir/src/cliqhash.c.o: CMakeFiles/blossom_separation.dir/flags.make
CMakeFiles/blossom_separation.dir/src/cliqhash.c.o: /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/cliqhash.c
CMakeFiles/blossom_separation.dir/src/cliqhash.c.o: CMakeFiles/blossom_separation.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/cvdg/Desktop/TSP_test/blossom_inequality_sep/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building C object CMakeFiles/blossom_separation.dir/src/cliqhash.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/blossom_separation.dir/src/cliqhash.c.o -MF CMakeFiles/blossom_separation.dir/src/cliqhash.c.o.d -o CMakeFiles/blossom_separation.dir/src/cliqhash.c.o -c /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/cliqhash.c

CMakeFiles/blossom_separation.dir/src/cliqhash.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/blossom_separation.dir/src/cliqhash.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/cliqhash.c > CMakeFiles/blossom_separation.dir/src/cliqhash.c.i

CMakeFiles/blossom_separation.dir/src/cliqhash.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/blossom_separation.dir/src/cliqhash.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/cliqhash.c -o CMakeFiles/blossom_separation.dir/src/cliqhash.c.s

CMakeFiles/blossom_separation.dir/src/genhash.c.o: CMakeFiles/blossom_separation.dir/flags.make
CMakeFiles/blossom_separation.dir/src/genhash.c.o: /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/genhash.c
CMakeFiles/blossom_separation.dir/src/genhash.c.o: CMakeFiles/blossom_separation.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/cvdg/Desktop/TSP_test/blossom_inequality_sep/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building C object CMakeFiles/blossom_separation.dir/src/genhash.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/blossom_separation.dir/src/genhash.c.o -MF CMakeFiles/blossom_separation.dir/src/genhash.c.o.d -o CMakeFiles/blossom_separation.dir/src/genhash.c.o -c /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/genhash.c

CMakeFiles/blossom_separation.dir/src/genhash.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/blossom_separation.dir/src/genhash.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/genhash.c > CMakeFiles/blossom_separation.dir/src/genhash.c.i

CMakeFiles/blossom_separation.dir/src/genhash.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/blossom_separation.dir/src/genhash.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/genhash.c -o CMakeFiles/blossom_separation.dir/src/genhash.c.s

CMakeFiles/blossom_separation.dir/src/cut_st.c.o: CMakeFiles/blossom_separation.dir/flags.make
CMakeFiles/blossom_separation.dir/src/cut_st.c.o: /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/cut_st.c
CMakeFiles/blossom_separation.dir/src/cut_st.c.o: CMakeFiles/blossom_separation.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/cvdg/Desktop/TSP_test/blossom_inequality_sep/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building C object CMakeFiles/blossom_separation.dir/src/cut_st.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/blossom_separation.dir/src/cut_st.c.o -MF CMakeFiles/blossom_separation.dir/src/cut_st.c.o.d -o CMakeFiles/blossom_separation.dir/src/cut_st.c.o -c /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/cut_st.c

CMakeFiles/blossom_separation.dir/src/cut_st.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/blossom_separation.dir/src/cut_st.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/cut_st.c > CMakeFiles/blossom_separation.dir/src/cut_st.c.i

CMakeFiles/blossom_separation.dir/src/cut_st.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/blossom_separation.dir/src/cut_st.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/cut_st.c -o CMakeFiles/blossom_separation.dir/src/cut_st.c.s

CMakeFiles/blossom_separation.dir/src/util.c.o: CMakeFiles/blossom_separation.dir/flags.make
CMakeFiles/blossom_separation.dir/src/util.c.o: /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/util.c
CMakeFiles/blossom_separation.dir/src/util.c.o: CMakeFiles/blossom_separation.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/cvdg/Desktop/TSP_test/blossom_inequality_sep/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Building C object CMakeFiles/blossom_separation.dir/src/util.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/blossom_separation.dir/src/util.c.o -MF CMakeFiles/blossom_separation.dir/src/util.c.o.d -o CMakeFiles/blossom_separation.dir/src/util.c.o -c /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/util.c

CMakeFiles/blossom_separation.dir/src/util.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/blossom_separation.dir/src/util.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/util.c > CMakeFiles/blossom_separation.dir/src/util.c.i

CMakeFiles/blossom_separation.dir/src/util.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/blossom_separation.dir/src/util.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/src/util.c -o CMakeFiles/blossom_separation.dir/src/util.c.s

# Object files for target blossom_separation
blossom_separation_OBJECTS = \
"CMakeFiles/blossom_separation.dir/blossom_detector.c.o" \
"CMakeFiles/blossom_separation.dir/src/blossom.c.o" \
"CMakeFiles/blossom_separation.dir/src/cliqwork.c.o" \
"CMakeFiles/blossom_separation.dir/src/skeleton.c.o" \
"CMakeFiles/blossom_separation.dir/src/cutpool.c.o" \
"CMakeFiles/blossom_separation.dir/src/allocrus.c.o" \
"CMakeFiles/blossom_separation.dir/src/urandom.c.o" \
"CMakeFiles/blossom_separation.dir/src/gomoryhu.c.o" \
"CMakeFiles/blossom_separation.dir/src/sortrus.c.o" \
"CMakeFiles/blossom_separation.dir/src/safe_io.c.o" \
"CMakeFiles/blossom_separation.dir/src/cliqhash.c.o" \
"CMakeFiles/blossom_separation.dir/src/genhash.c.o" \
"CMakeFiles/blossom_separation.dir/src/cut_st.c.o" \
"CMakeFiles/blossom_separation.dir/src/util.c.o"

# External object files for target blossom_separation
blossom_separation_EXTERNAL_OBJECTS =

lib/libblossom_separation.so: CMakeFiles/blossom_separation.dir/blossom_detector.c.o
lib/libblossom_separation.so: CMakeFiles/blossom_separation.dir/src/blossom.c.o
lib/libblossom_separation.so: CMakeFiles/blossom_separation.dir/src/cliqwork.c.o
lib/libblossom_separation.so: CMakeFiles/blossom_separation.dir/src/skeleton.c.o
lib/libblossom_separation.so: CMakeFiles/blossom_separation.dir/src/cutpool.c.o
lib/libblossom_separation.so: CMakeFiles/blossom_separation.dir/src/allocrus.c.o
lib/libblossom_separation.so: CMakeFiles/blossom_separation.dir/src/urandom.c.o
lib/libblossom_separation.so: CMakeFiles/blossom_separation.dir/src/gomoryhu.c.o
lib/libblossom_separation.so: CMakeFiles/blossom_separation.dir/src/sortrus.c.o
lib/libblossom_separation.so: CMakeFiles/blossom_separation.dir/src/safe_io.c.o
lib/libblossom_separation.so: CMakeFiles/blossom_separation.dir/src/cliqhash.c.o
lib/libblossom_separation.so: CMakeFiles/blossom_separation.dir/src/genhash.c.o
lib/libblossom_separation.so: CMakeFiles/blossom_separation.dir/src/cut_st.c.o
lib/libblossom_separation.so: CMakeFiles/blossom_separation.dir/src/util.c.o
lib/libblossom_separation.so: CMakeFiles/blossom_separation.dir/build.make
lib/libblossom_separation.so: CMakeFiles/blossom_separation.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/cvdg/Desktop/TSP_test/blossom_inequality_sep/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_15) "Linking C shared library lib/libblossom_separation.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/blossom_separation.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/blossom_separation.dir/build: lib/libblossom_separation.so
.PHONY : CMakeFiles/blossom_separation.dir/build

CMakeFiles/blossom_separation.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/blossom_separation.dir/cmake_clean.cmake
.PHONY : CMakeFiles/blossom_separation.dir/clean

CMakeFiles/blossom_separation.dir/depend:
	cd /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/cvdg/Desktop/TSP_test/blossom_inequality_sep /home/cvdg/Desktop/TSP_test/blossom_inequality_sep /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/build /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/build /home/cvdg/Desktop/TSP_test/blossom_inequality_sep/build/CMakeFiles/blossom_separation.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/blossom_separation.dir/depend

