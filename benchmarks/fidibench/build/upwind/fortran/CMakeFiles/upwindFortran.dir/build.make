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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/tolnaia/COMBS/benchmarks/fidibench

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/tolnaia/COMBS/benchmarks/fidibench/build

# Include any dependencies generated for this target.
include upwind/fortran/CMakeFiles/upwindFortran.dir/depend.make

# Include the progress variables for this target.
include upwind/fortran/CMakeFiles/upwindFortran.dir/progress.make

# Include the compile flags for this target's objects.
include upwind/fortran/CMakeFiles/upwindFortran.dir/flags.make

upwind/fortran/CMakeFiles/upwindFortran.dir/upwind.F90.o: upwind/fortran/CMakeFiles/upwindFortran.dir/flags.make
upwind/fortran/CMakeFiles/upwindFortran.dir/upwind.F90.o: ../upwind/fortran/upwind.F90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/tolnaia/COMBS/benchmarks/fidibench/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object upwind/fortran/CMakeFiles/upwindFortran.dir/upwind.F90.o"
	cd /home/tolnaia/COMBS/benchmarks/fidibench/build/upwind/fortran && /opt/moose/mpich-3.3/gcc-9.2.0/bin/mpif90 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/tolnaia/COMBS/benchmarks/fidibench/upwind/fortran/upwind.F90 -o CMakeFiles/upwindFortran.dir/upwind.F90.o

upwind/fortran/CMakeFiles/upwindFortran.dir/upwind.F90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/upwindFortran.dir/upwind.F90.i"
	cd /home/tolnaia/COMBS/benchmarks/fidibench/build/upwind/fortran && /opt/moose/mpich-3.3/gcc-9.2.0/bin/mpif90 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/tolnaia/COMBS/benchmarks/fidibench/upwind/fortran/upwind.F90 > CMakeFiles/upwindFortran.dir/upwind.F90.i

upwind/fortran/CMakeFiles/upwindFortran.dir/upwind.F90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/upwindFortran.dir/upwind.F90.s"
	cd /home/tolnaia/COMBS/benchmarks/fidibench/build/upwind/fortran && /opt/moose/mpich-3.3/gcc-9.2.0/bin/mpif90 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/tolnaia/COMBS/benchmarks/fidibench/upwind/fortran/upwind.F90 -o CMakeFiles/upwindFortran.dir/upwind.F90.s

upwind/fortran/CMakeFiles/upwindFortran.dir/upwind.F90.o.requires:

.PHONY : upwind/fortran/CMakeFiles/upwindFortran.dir/upwind.F90.o.requires

upwind/fortran/CMakeFiles/upwindFortran.dir/upwind.F90.o.provides: upwind/fortran/CMakeFiles/upwindFortran.dir/upwind.F90.o.requires
	$(MAKE) -f upwind/fortran/CMakeFiles/upwindFortran.dir/build.make upwind/fortran/CMakeFiles/upwindFortran.dir/upwind.F90.o.provides.build
.PHONY : upwind/fortran/CMakeFiles/upwindFortran.dir/upwind.F90.o.provides

upwind/fortran/CMakeFiles/upwindFortran.dir/upwind.F90.o.provides.build: upwind/fortran/CMakeFiles/upwindFortran.dir/upwind.F90.o


# Object files for target upwindFortran
upwindFortran_OBJECTS = \
"CMakeFiles/upwindFortran.dir/upwind.F90.o"

# External object files for target upwindFortran
upwindFortran_EXTERNAL_OBJECTS =

upwind/fortran/upwindFortran: upwind/fortran/CMakeFiles/upwindFortran.dir/upwind.F90.o
upwind/fortran/upwindFortran: upwind/fortran/CMakeFiles/upwindFortran.dir/build.make
upwind/fortran/upwindFortran: upwind/fortran/CMakeFiles/upwindFortran.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/tolnaia/COMBS/benchmarks/fidibench/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking Fortran executable upwindFortran"
	cd /home/tolnaia/COMBS/benchmarks/fidibench/build/upwind/fortran && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/upwindFortran.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
upwind/fortran/CMakeFiles/upwindFortran.dir/build: upwind/fortran/upwindFortran

.PHONY : upwind/fortran/CMakeFiles/upwindFortran.dir/build

upwind/fortran/CMakeFiles/upwindFortran.dir/requires: upwind/fortran/CMakeFiles/upwindFortran.dir/upwind.F90.o.requires

.PHONY : upwind/fortran/CMakeFiles/upwindFortran.dir/requires

upwind/fortran/CMakeFiles/upwindFortran.dir/clean:
	cd /home/tolnaia/COMBS/benchmarks/fidibench/build/upwind/fortran && $(CMAKE_COMMAND) -P CMakeFiles/upwindFortran.dir/cmake_clean.cmake
.PHONY : upwind/fortran/CMakeFiles/upwindFortran.dir/clean

upwind/fortran/CMakeFiles/upwindFortran.dir/depend:
	cd /home/tolnaia/COMBS/benchmarks/fidibench/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/tolnaia/COMBS/benchmarks/fidibench /home/tolnaia/COMBS/benchmarks/fidibench/upwind/fortran /home/tolnaia/COMBS/benchmarks/fidibench/build /home/tolnaia/COMBS/benchmarks/fidibench/build/upwind/fortran /home/tolnaia/COMBS/benchmarks/fidibench/build/upwind/fortran/CMakeFiles/upwindFortran.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : upwind/fortran/CMakeFiles/upwindFortran.dir/depend
