# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.7

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


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
CMAKE_COMMAND = /opt/software/cmake/3.7.2--GCC-4.4.5/bin/cmake

# The command to remove a file.
RM = /opt/software/cmake/3.7.2--GCC-4.4.5/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /mnt/home/dunan/Job/2018/GroupK/GroupK/yass

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/home/dunan/Job/2018/GroupK/GroupK/yass

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/opt/software/cmake/3.7.2--GCC-4.4.5/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/opt/software/cmake/3.7.2--GCC-4.4.5/bin/ccmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /mnt/home/dunan/Job/2018/GroupK/GroupK/yass/CMakeFiles /mnt/home/dunan/Job/2018/GroupK/GroupK/yass/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /mnt/home/dunan/Job/2018/GroupK/GroupK/yass/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named yass

# Build rule for target.
yass: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 yass
.PHONY : yass

# fast build rule for target.
yass/fast:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/build
.PHONY : yass/fast

src/PyYASS.o: src/PyYASS.cpp.o

.PHONY : src/PyYASS.o

# target to build an object file
src/PyYASS.cpp.o:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/PyYASS.cpp.o
.PHONY : src/PyYASS.cpp.o

src/PyYASS.i: src/PyYASS.cpp.i

.PHONY : src/PyYASS.i

# target to preprocess a source file
src/PyYASS.cpp.i:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/PyYASS.cpp.i
.PHONY : src/PyYASS.cpp.i

src/PyYASS.s: src/PyYASS.cpp.s

.PHONY : src/PyYASS.s

# target to generate assembly for a file
src/PyYASS.cpp.s:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/PyYASS.cpp.s
.PHONY : src/PyYASS.cpp.s

src/align.o: src/align.c.o

.PHONY : src/align.o

# target to build an object file
src/align.c.o:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/align.c.o
.PHONY : src/align.c.o

src/align.i: src/align.c.i

.PHONY : src/align.i

# target to preprocess a source file
src/align.c.i:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/align.c.i
.PHONY : src/align.c.i

src/align.s: src/align.c.s

.PHONY : src/align.s

# target to generate assembly for a file
src/align.c.s:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/align.c.s
.PHONY : src/align.c.s

src/assemble.o: src/assemble.c.o

.PHONY : src/assemble.o

# target to build an object file
src/assemble.c.o:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/assemble.c.o
.PHONY : src/assemble.c.o

src/assemble.i: src/assemble.c.i

.PHONY : src/assemble.i

# target to preprocess a source file
src/assemble.c.i:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/assemble.c.i
.PHONY : src/assemble.c.i

src/assemble.s: src/assemble.c.s

.PHONY : src/assemble.s

# target to generate assembly for a file
src/assemble.c.s:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/assemble.c.s
.PHONY : src/assemble.c.s

src/avl.o: src/avl.c.o

.PHONY : src/avl.o

# target to build an object file
src/avl.c.o:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/avl.c.o
.PHONY : src/avl.c.o

src/avl.i: src/avl.c.i

.PHONY : src/avl.i

# target to preprocess a source file
src/avl.c.i:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/avl.c.i
.PHONY : src/avl.c.i

src/avl.s: src/avl.c.s

.PHONY : src/avl.s

# target to generate assembly for a file
src/avl.c.s:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/avl.c.s
.PHONY : src/avl.c.s

src/display.o: src/display.c.o

.PHONY : src/display.o

# target to build an object file
src/display.c.o:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/display.c.o
.PHONY : src/display.c.o

src/display.i: src/display.c.i

.PHONY : src/display.i

# target to preprocess a source file
src/display.c.i:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/display.c.i
.PHONY : src/display.c.i

src/display.s: src/display.c.s

.PHONY : src/display.s

# target to generate assembly for a file
src/display.c.s:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/display.c.s
.PHONY : src/display.c.s

src/global_var.o: src/global_var.c.o

.PHONY : src/global_var.o

# target to build an object file
src/global_var.c.o:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/global_var.c.o
.PHONY : src/global_var.c.o

src/global_var.i: src/global_var.c.i

.PHONY : src/global_var.i

# target to preprocess a source file
src/global_var.c.i:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/global_var.c.i
.PHONY : src/global_var.c.i

src/global_var.s: src/global_var.c.s

.PHONY : src/global_var.s

# target to generate assembly for a file
src/global_var.c.s:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/global_var.c.s
.PHONY : src/global_var.c.s

src/kword.o: src/kword.c.o

.PHONY : src/kword.o

# target to build an object file
src/kword.c.o:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/kword.c.o
.PHONY : src/kword.c.o

src/kword.i: src/kword.c.i

.PHONY : src/kword.i

# target to preprocess a source file
src/kword.c.i:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/kword.c.i
.PHONY : src/kword.c.i

src/kword.s: src/kword.c.s

.PHONY : src/kword.s

# target to generate assembly for a file
src/kword.c.s:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/kword.c.s
.PHONY : src/kword.c.s

src/list.o: src/list.c.o

.PHONY : src/list.o

# target to build an object file
src/list.c.o:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/list.c.o
.PHONY : src/list.c.o

src/list.i: src/list.c.i

.PHONY : src/list.i

# target to preprocess a source file
src/list.c.i:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/list.c.i
.PHONY : src/list.c.i

src/list.s: src/list.c.s

.PHONY : src/list.s

# target to generate assembly for a file
src/list.c.s:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/list.c.s
.PHONY : src/list.c.s

src/main.o: src/main.c.o

.PHONY : src/main.o

# target to build an object file
src/main.c.o:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/main.c.o
.PHONY : src/main.c.o

src/main.i: src/main.c.i

.PHONY : src/main.i

# target to preprocess a source file
src/main.c.i:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/main.c.i
.PHONY : src/main.c.i

src/main.s: src/main.c.s

.PHONY : src/main.s

# target to generate assembly for a file
src/main.c.s:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/main.c.s
.PHONY : src/main.c.s

src/prdyn.o: src/prdyn.c.o

.PHONY : src/prdyn.o

# target to build an object file
src/prdyn.c.o:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/prdyn.c.o
.PHONY : src/prdyn.c.o

src/prdyn.i: src/prdyn.c.i

.PHONY : src/prdyn.i

# target to preprocess a source file
src/prdyn.c.i:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/prdyn.c.i
.PHONY : src/prdyn.c.i

src/prdyn.s: src/prdyn.c.s

.PHONY : src/prdyn.s

# target to generate assembly for a file
src/prdyn.c.s:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/prdyn.c.s
.PHONY : src/prdyn.c.s

src/proba.o: src/proba.c.o

.PHONY : src/proba.o

# target to build an object file
src/proba.c.o:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/proba.c.o
.PHONY : src/proba.c.o

src/proba.i: src/proba.c.i

.PHONY : src/proba.i

# target to preprocess a source file
src/proba.c.i:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/proba.c.i
.PHONY : src/proba.c.i

src/proba.s: src/proba.c.s

.PHONY : src/proba.s

# target to generate assembly for a file
src/proba.c.s:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/proba.c.s
.PHONY : src/proba.c.s

src/red_black.o: src/red_black.c.o

.PHONY : src/red_black.o

# target to build an object file
src/red_black.c.o:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/red_black.c.o
.PHONY : src/red_black.c.o

src/red_black.i: src/red_black.c.i

.PHONY : src/red_black.i

# target to preprocess a source file
src/red_black.c.i:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/red_black.c.i
.PHONY : src/red_black.c.i

src/red_black.s: src/red_black.c.s

.PHONY : src/red_black.s

# target to generate assembly for a file
src/red_black.c.s:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/red_black.c.s
.PHONY : src/red_black.c.s

src/regroup.o: src/regroup.c.o

.PHONY : src/regroup.o

# target to build an object file
src/regroup.c.o:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/regroup.c.o
.PHONY : src/regroup.c.o

src/regroup.i: src/regroup.c.i

.PHONY : src/regroup.i

# target to preprocess a source file
src/regroup.c.i:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/regroup.c.i
.PHONY : src/regroup.c.i

src/regroup.s: src/regroup.c.s

.PHONY : src/regroup.s

# target to generate assembly for a file
src/regroup.c.s:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/regroup.c.s
.PHONY : src/regroup.c.s

src/threads.o: src/threads.c.o

.PHONY : src/threads.o

# target to build an object file
src/threads.c.o:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/threads.c.o
.PHONY : src/threads.c.o

src/threads.i: src/threads.c.i

.PHONY : src/threads.i

# target to preprocess a source file
src/threads.c.i:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/threads.c.i
.PHONY : src/threads.c.i

src/threads.s: src/threads.c.s

.PHONY : src/threads.s

# target to generate assembly for a file
src/threads.c.s:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/threads.c.s
.PHONY : src/threads.c.s

src/tuple.o: src/tuple.c.o

.PHONY : src/tuple.o

# target to build an object file
src/tuple.c.o:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/tuple.c.o
.PHONY : src/tuple.c.o

src/tuple.i: src/tuple.c.i

.PHONY : src/tuple.i

# target to preprocess a source file
src/tuple.c.i:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/tuple.c.i
.PHONY : src/tuple.c.i

src/tuple.s: src/tuple.c.s

.PHONY : src/tuple.s

# target to generate assembly for a file
src/tuple.c.s:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/tuple.c.s
.PHONY : src/tuple.c.s

src/util.o: src/util.c.o

.PHONY : src/util.o

# target to build an object file
src/util.c.o:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/util.c.o
.PHONY : src/util.c.o

src/util.i: src/util.c.i

.PHONY : src/util.i

# target to preprocess a source file
src/util.c.i:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/util.c.i
.PHONY : src/util.c.i

src/util.s: src/util.c.s

.PHONY : src/util.s

# target to generate assembly for a file
src/util.c.s:
	$(MAKE) -f CMakeFiles/yass.dir/build.make CMakeFiles/yass.dir/src/util.c.s
.PHONY : src/util.c.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... yass"
	@echo "... rebuild_cache"
	@echo "... edit_cache"
	@echo "... src/PyYASS.o"
	@echo "... src/PyYASS.i"
	@echo "... src/PyYASS.s"
	@echo "... src/align.o"
	@echo "... src/align.i"
	@echo "... src/align.s"
	@echo "... src/assemble.o"
	@echo "... src/assemble.i"
	@echo "... src/assemble.s"
	@echo "... src/avl.o"
	@echo "... src/avl.i"
	@echo "... src/avl.s"
	@echo "... src/display.o"
	@echo "... src/display.i"
	@echo "... src/display.s"
	@echo "... src/global_var.o"
	@echo "... src/global_var.i"
	@echo "... src/global_var.s"
	@echo "... src/kword.o"
	@echo "... src/kword.i"
	@echo "... src/kword.s"
	@echo "... src/list.o"
	@echo "... src/list.i"
	@echo "... src/list.s"
	@echo "... src/main.o"
	@echo "... src/main.i"
	@echo "... src/main.s"
	@echo "... src/prdyn.o"
	@echo "... src/prdyn.i"
	@echo "... src/prdyn.s"
	@echo "... src/proba.o"
	@echo "... src/proba.i"
	@echo "... src/proba.s"
	@echo "... src/red_black.o"
	@echo "... src/red_black.i"
	@echo "... src/red_black.s"
	@echo "... src/regroup.o"
	@echo "... src/regroup.i"
	@echo "... src/regroup.s"
	@echo "... src/threads.o"
	@echo "... src/threads.i"
	@echo "... src/threads.s"
	@echo "... src/tuple.o"
	@echo "... src/tuple.i"
	@echo "... src/tuple.s"
	@echo "... src/util.o"
	@echo "... src/util.i"
	@echo "... src/util.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

