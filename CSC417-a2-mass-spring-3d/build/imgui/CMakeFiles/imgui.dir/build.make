# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.27

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
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.27.5/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.27.5/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/build

# Include any dependencies generated for this target.
include imgui/CMakeFiles/imgui.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include imgui/CMakeFiles/imgui.dir/compiler_depend.make

# Include the progress variables for this target.
include imgui/CMakeFiles/imgui.dir/progress.make

# Include the compile flags for this target's objects.
include imgui/CMakeFiles/imgui.dir/flags.make

imgui/CMakeFiles/imgui.dir/__/imgui/imgui.cpp.o: imgui/CMakeFiles/imgui.dir/flags.make
imgui/CMakeFiles/imgui.dir/__/imgui/imgui.cpp.o: /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/extern/libigl/external/imgui/imgui.cpp
imgui/CMakeFiles/imgui.dir/__/imgui/imgui.cpp.o: imgui/CMakeFiles/imgui.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object imgui/CMakeFiles/imgui.dir/__/imgui/imgui.cpp.o"
	cd /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/build/imgui && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT imgui/CMakeFiles/imgui.dir/__/imgui/imgui.cpp.o -MF CMakeFiles/imgui.dir/__/imgui/imgui.cpp.o.d -o CMakeFiles/imgui.dir/__/imgui/imgui.cpp.o -c /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/extern/libigl/external/imgui/imgui.cpp

imgui/CMakeFiles/imgui.dir/__/imgui/imgui.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/imgui.dir/__/imgui/imgui.cpp.i"
	cd /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/build/imgui && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/extern/libigl/external/imgui/imgui.cpp > CMakeFiles/imgui.dir/__/imgui/imgui.cpp.i

imgui/CMakeFiles/imgui.dir/__/imgui/imgui.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/imgui.dir/__/imgui/imgui.cpp.s"
	cd /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/build/imgui && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/extern/libigl/external/imgui/imgui.cpp -o CMakeFiles/imgui.dir/__/imgui/imgui.cpp.s

imgui/CMakeFiles/imgui.dir/__/imgui/imgui_demo.cpp.o: imgui/CMakeFiles/imgui.dir/flags.make
imgui/CMakeFiles/imgui.dir/__/imgui/imgui_demo.cpp.o: /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/extern/libigl/external/imgui/imgui_demo.cpp
imgui/CMakeFiles/imgui.dir/__/imgui/imgui_demo.cpp.o: imgui/CMakeFiles/imgui.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object imgui/CMakeFiles/imgui.dir/__/imgui/imgui_demo.cpp.o"
	cd /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/build/imgui && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT imgui/CMakeFiles/imgui.dir/__/imgui/imgui_demo.cpp.o -MF CMakeFiles/imgui.dir/__/imgui/imgui_demo.cpp.o.d -o CMakeFiles/imgui.dir/__/imgui/imgui_demo.cpp.o -c /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/extern/libigl/external/imgui/imgui_demo.cpp

imgui/CMakeFiles/imgui.dir/__/imgui/imgui_demo.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/imgui.dir/__/imgui/imgui_demo.cpp.i"
	cd /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/build/imgui && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/extern/libigl/external/imgui/imgui_demo.cpp > CMakeFiles/imgui.dir/__/imgui/imgui_demo.cpp.i

imgui/CMakeFiles/imgui.dir/__/imgui/imgui_demo.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/imgui.dir/__/imgui/imgui_demo.cpp.s"
	cd /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/build/imgui && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/extern/libigl/external/imgui/imgui_demo.cpp -o CMakeFiles/imgui.dir/__/imgui/imgui_demo.cpp.s

imgui/CMakeFiles/imgui.dir/__/imgui/imgui_draw.cpp.o: imgui/CMakeFiles/imgui.dir/flags.make
imgui/CMakeFiles/imgui.dir/__/imgui/imgui_draw.cpp.o: /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/extern/libigl/external/imgui/imgui_draw.cpp
imgui/CMakeFiles/imgui.dir/__/imgui/imgui_draw.cpp.o: imgui/CMakeFiles/imgui.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object imgui/CMakeFiles/imgui.dir/__/imgui/imgui_draw.cpp.o"
	cd /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/build/imgui && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT imgui/CMakeFiles/imgui.dir/__/imgui/imgui_draw.cpp.o -MF CMakeFiles/imgui.dir/__/imgui/imgui_draw.cpp.o.d -o CMakeFiles/imgui.dir/__/imgui/imgui_draw.cpp.o -c /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/extern/libigl/external/imgui/imgui_draw.cpp

imgui/CMakeFiles/imgui.dir/__/imgui/imgui_draw.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/imgui.dir/__/imgui/imgui_draw.cpp.i"
	cd /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/build/imgui && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/extern/libigl/external/imgui/imgui_draw.cpp > CMakeFiles/imgui.dir/__/imgui/imgui_draw.cpp.i

imgui/CMakeFiles/imgui.dir/__/imgui/imgui_draw.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/imgui.dir/__/imgui/imgui_draw.cpp.s"
	cd /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/build/imgui && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/extern/libigl/external/imgui/imgui_draw.cpp -o CMakeFiles/imgui.dir/__/imgui/imgui_draw.cpp.s

imgui/CMakeFiles/imgui.dir/__/imgui/imgui_widgets.cpp.o: imgui/CMakeFiles/imgui.dir/flags.make
imgui/CMakeFiles/imgui.dir/__/imgui/imgui_widgets.cpp.o: /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/extern/libigl/external/imgui/imgui_widgets.cpp
imgui/CMakeFiles/imgui.dir/__/imgui/imgui_widgets.cpp.o: imgui/CMakeFiles/imgui.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object imgui/CMakeFiles/imgui.dir/__/imgui/imgui_widgets.cpp.o"
	cd /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/build/imgui && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT imgui/CMakeFiles/imgui.dir/__/imgui/imgui_widgets.cpp.o -MF CMakeFiles/imgui.dir/__/imgui/imgui_widgets.cpp.o.d -o CMakeFiles/imgui.dir/__/imgui/imgui_widgets.cpp.o -c /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/extern/libigl/external/imgui/imgui_widgets.cpp

imgui/CMakeFiles/imgui.dir/__/imgui/imgui_widgets.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/imgui.dir/__/imgui/imgui_widgets.cpp.i"
	cd /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/build/imgui && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/extern/libigl/external/imgui/imgui_widgets.cpp > CMakeFiles/imgui.dir/__/imgui/imgui_widgets.cpp.i

imgui/CMakeFiles/imgui.dir/__/imgui/imgui_widgets.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/imgui.dir/__/imgui/imgui_widgets.cpp.s"
	cd /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/build/imgui && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/extern/libigl/external/imgui/imgui_widgets.cpp -o CMakeFiles/imgui.dir/__/imgui/imgui_widgets.cpp.s

imgui/CMakeFiles/imgui.dir/__/imgui/examples/imgui_impl_glfw.cpp.o: imgui/CMakeFiles/imgui.dir/flags.make
imgui/CMakeFiles/imgui.dir/__/imgui/examples/imgui_impl_glfw.cpp.o: /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/extern/libigl/external/imgui/examples/imgui_impl_glfw.cpp
imgui/CMakeFiles/imgui.dir/__/imgui/examples/imgui_impl_glfw.cpp.o: imgui/CMakeFiles/imgui.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object imgui/CMakeFiles/imgui.dir/__/imgui/examples/imgui_impl_glfw.cpp.o"
	cd /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/build/imgui && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT imgui/CMakeFiles/imgui.dir/__/imgui/examples/imgui_impl_glfw.cpp.o -MF CMakeFiles/imgui.dir/__/imgui/examples/imgui_impl_glfw.cpp.o.d -o CMakeFiles/imgui.dir/__/imgui/examples/imgui_impl_glfw.cpp.o -c /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/extern/libigl/external/imgui/examples/imgui_impl_glfw.cpp

imgui/CMakeFiles/imgui.dir/__/imgui/examples/imgui_impl_glfw.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/imgui.dir/__/imgui/examples/imgui_impl_glfw.cpp.i"
	cd /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/build/imgui && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/extern/libigl/external/imgui/examples/imgui_impl_glfw.cpp > CMakeFiles/imgui.dir/__/imgui/examples/imgui_impl_glfw.cpp.i

imgui/CMakeFiles/imgui.dir/__/imgui/examples/imgui_impl_glfw.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/imgui.dir/__/imgui/examples/imgui_impl_glfw.cpp.s"
	cd /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/build/imgui && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/extern/libigl/external/imgui/examples/imgui_impl_glfw.cpp -o CMakeFiles/imgui.dir/__/imgui/examples/imgui_impl_glfw.cpp.s

imgui/CMakeFiles/imgui.dir/__/imgui/examples/imgui_impl_opengl3.cpp.o: imgui/CMakeFiles/imgui.dir/flags.make
imgui/CMakeFiles/imgui.dir/__/imgui/examples/imgui_impl_opengl3.cpp.o: /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/extern/libigl/external/imgui/examples/imgui_impl_opengl3.cpp
imgui/CMakeFiles/imgui.dir/__/imgui/examples/imgui_impl_opengl3.cpp.o: imgui/CMakeFiles/imgui.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object imgui/CMakeFiles/imgui.dir/__/imgui/examples/imgui_impl_opengl3.cpp.o"
	cd /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/build/imgui && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT imgui/CMakeFiles/imgui.dir/__/imgui/examples/imgui_impl_opengl3.cpp.o -MF CMakeFiles/imgui.dir/__/imgui/examples/imgui_impl_opengl3.cpp.o.d -o CMakeFiles/imgui.dir/__/imgui/examples/imgui_impl_opengl3.cpp.o -c /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/extern/libigl/external/imgui/examples/imgui_impl_opengl3.cpp

imgui/CMakeFiles/imgui.dir/__/imgui/examples/imgui_impl_opengl3.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/imgui.dir/__/imgui/examples/imgui_impl_opengl3.cpp.i"
	cd /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/build/imgui && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/extern/libigl/external/imgui/examples/imgui_impl_opengl3.cpp > CMakeFiles/imgui.dir/__/imgui/examples/imgui_impl_opengl3.cpp.i

imgui/CMakeFiles/imgui.dir/__/imgui/examples/imgui_impl_opengl3.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/imgui.dir/__/imgui/examples/imgui_impl_opengl3.cpp.s"
	cd /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/build/imgui && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/extern/libigl/external/imgui/examples/imgui_impl_opengl3.cpp -o CMakeFiles/imgui.dir/__/imgui/examples/imgui_impl_opengl3.cpp.s

# Object files for target imgui
imgui_OBJECTS = \
"CMakeFiles/imgui.dir/__/imgui/imgui.cpp.o" \
"CMakeFiles/imgui.dir/__/imgui/imgui_demo.cpp.o" \
"CMakeFiles/imgui.dir/__/imgui/imgui_draw.cpp.o" \
"CMakeFiles/imgui.dir/__/imgui/imgui_widgets.cpp.o" \
"CMakeFiles/imgui.dir/__/imgui/examples/imgui_impl_glfw.cpp.o" \
"CMakeFiles/imgui.dir/__/imgui/examples/imgui_impl_opengl3.cpp.o"

# External object files for target imgui
imgui_EXTERNAL_OBJECTS =

imgui/libimgui.a: imgui/CMakeFiles/imgui.dir/__/imgui/imgui.cpp.o
imgui/libimgui.a: imgui/CMakeFiles/imgui.dir/__/imgui/imgui_demo.cpp.o
imgui/libimgui.a: imgui/CMakeFiles/imgui.dir/__/imgui/imgui_draw.cpp.o
imgui/libimgui.a: imgui/CMakeFiles/imgui.dir/__/imgui/imgui_widgets.cpp.o
imgui/libimgui.a: imgui/CMakeFiles/imgui.dir/__/imgui/examples/imgui_impl_glfw.cpp.o
imgui/libimgui.a: imgui/CMakeFiles/imgui.dir/__/imgui/examples/imgui_impl_opengl3.cpp.o
imgui/libimgui.a: imgui/CMakeFiles/imgui.dir/build.make
imgui/libimgui.a: imgui/CMakeFiles/imgui.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Linking CXX static library libimgui.a"
	cd /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/build/imgui && $(CMAKE_COMMAND) -P CMakeFiles/imgui.dir/cmake_clean_target.cmake
	cd /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/build/imgui && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/imgui.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
imgui/CMakeFiles/imgui.dir/build: imgui/libimgui.a
.PHONY : imgui/CMakeFiles/imgui.dir/build

imgui/CMakeFiles/imgui.dir/clean:
	cd /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/build/imgui && $(CMAKE_COMMAND) -P CMakeFiles/imgui.dir/cmake_clean.cmake
.PHONY : imgui/CMakeFiles/imgui.dir/clean

imgui/CMakeFiles/imgui.dir/depend:
	cd /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/extern/libigl/external/libigl-imgui /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/build /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/build/imgui /Users/natalia/Dropbox/CSQUEEN/CSC2549-Phys-based-animation/CSC417-a2-mass-spring-3d/build/imgui/CMakeFiles/imgui.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : imgui/CMakeFiles/imgui.dir/depend

