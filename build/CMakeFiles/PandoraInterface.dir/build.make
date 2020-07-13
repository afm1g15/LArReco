# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /var/clus/usera/afm67/2020/May/CosmicJorisMay2020/LArReco

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /var/clus/usera/afm67/2020/May/CosmicJorisMay2020/LArReco/build

# Include any dependencies generated for this target.
include CMakeFiles/PandoraInterface.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/PandoraInterface.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/PandoraInterface.dir/flags.make

CMakeFiles/PandoraInterface.dir/test/PandoraInterface.cxx.o: CMakeFiles/PandoraInterface.dir/flags.make
CMakeFiles/PandoraInterface.dir/test/PandoraInterface.cxx.o: ../test/PandoraInterface.cxx
	$(CMAKE_COMMAND) -E cmake_progress_report /var/clus/usera/afm67/2020/May/CosmicJorisMay2020/LArReco/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/PandoraInterface.dir/test/PandoraInterface.cxx.o"
	/cvmfs/larsoft.opensciencegrid.org/products/gcc/v8_2_0/Linux64bit+3.10-2.17/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/PandoraInterface.dir/test/PandoraInterface.cxx.o -c /var/clus/usera/afm67/2020/May/CosmicJorisMay2020/LArReco/test/PandoraInterface.cxx

CMakeFiles/PandoraInterface.dir/test/PandoraInterface.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PandoraInterface.dir/test/PandoraInterface.cxx.i"
	/cvmfs/larsoft.opensciencegrid.org/products/gcc/v8_2_0/Linux64bit+3.10-2.17/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /var/clus/usera/afm67/2020/May/CosmicJorisMay2020/LArReco/test/PandoraInterface.cxx > CMakeFiles/PandoraInterface.dir/test/PandoraInterface.cxx.i

CMakeFiles/PandoraInterface.dir/test/PandoraInterface.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PandoraInterface.dir/test/PandoraInterface.cxx.s"
	/cvmfs/larsoft.opensciencegrid.org/products/gcc/v8_2_0/Linux64bit+3.10-2.17/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /var/clus/usera/afm67/2020/May/CosmicJorisMay2020/LArReco/test/PandoraInterface.cxx -o CMakeFiles/PandoraInterface.dir/test/PandoraInterface.cxx.s

CMakeFiles/PandoraInterface.dir/test/PandoraInterface.cxx.o.requires:
.PHONY : CMakeFiles/PandoraInterface.dir/test/PandoraInterface.cxx.o.requires

CMakeFiles/PandoraInterface.dir/test/PandoraInterface.cxx.o.provides: CMakeFiles/PandoraInterface.dir/test/PandoraInterface.cxx.o.requires
	$(MAKE) -f CMakeFiles/PandoraInterface.dir/build.make CMakeFiles/PandoraInterface.dir/test/PandoraInterface.cxx.o.provides.build
.PHONY : CMakeFiles/PandoraInterface.dir/test/PandoraInterface.cxx.o.provides

CMakeFiles/PandoraInterface.dir/test/PandoraInterface.cxx.o.provides.build: CMakeFiles/PandoraInterface.dir/test/PandoraInterface.cxx.o

# Object files for target PandoraInterface
PandoraInterface_OBJECTS = \
"CMakeFiles/PandoraInterface.dir/test/PandoraInterface.cxx.o"

# External object files for target PandoraInterface
PandoraInterface_EXTERNAL_OBJECTS =

bin/PandoraInterface: CMakeFiles/PandoraInterface.dir/test/PandoraInterface.cxx.o
bin/PandoraInterface: CMakeFiles/PandoraInterface.dir/build.make
bin/PandoraInterface: /var/clus/usera/afm67/2020/May/CosmicJorisMay2020/PandoraSDK/lib/libPandoraSDK.so
bin/PandoraInterface: /var/clus/usera/afm67/2020/May/CosmicJorisMay2020/LArContent/lib/libLArContent.so
bin/PandoraInterface: /var/clus/usera/afm67/2020/May/CosmicJorisMay2020/PandoraMonitoring/lib/libPandoraMonitoring.so
bin/PandoraInterface: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_18_04b/Linux64bit+3.10-2.17-e19-prof/lib/libCore.so
bin/PandoraInterface: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_18_04b/Linux64bit+3.10-2.17-e19-prof/lib/libImt.so
bin/PandoraInterface: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_18_04b/Linux64bit+3.10-2.17-e19-prof/lib/libRIO.so
bin/PandoraInterface: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_18_04b/Linux64bit+3.10-2.17-e19-prof/lib/libNet.so
bin/PandoraInterface: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_18_04b/Linux64bit+3.10-2.17-e19-prof/lib/libHist.so
bin/PandoraInterface: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_18_04b/Linux64bit+3.10-2.17-e19-prof/lib/libGraf.so
bin/PandoraInterface: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_18_04b/Linux64bit+3.10-2.17-e19-prof/lib/libGraf3d.so
bin/PandoraInterface: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_18_04b/Linux64bit+3.10-2.17-e19-prof/lib/libGpad.so
bin/PandoraInterface: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_18_04b/Linux64bit+3.10-2.17-e19-prof/lib/libROOTDataFrame.so
bin/PandoraInterface: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_18_04b/Linux64bit+3.10-2.17-e19-prof/lib/libTree.so
bin/PandoraInterface: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_18_04b/Linux64bit+3.10-2.17-e19-prof/lib/libTreePlayer.so
bin/PandoraInterface: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_18_04b/Linux64bit+3.10-2.17-e19-prof/lib/libRint.so
bin/PandoraInterface: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_18_04b/Linux64bit+3.10-2.17-e19-prof/lib/libPostscript.so
bin/PandoraInterface: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_18_04b/Linux64bit+3.10-2.17-e19-prof/lib/libMatrix.so
bin/PandoraInterface: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_18_04b/Linux64bit+3.10-2.17-e19-prof/lib/libPhysics.so
bin/PandoraInterface: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_18_04b/Linux64bit+3.10-2.17-e19-prof/lib/libMathCore.so
bin/PandoraInterface: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_18_04b/Linux64bit+3.10-2.17-e19-prof/lib/libThread.so
bin/PandoraInterface: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_18_04b/Linux64bit+3.10-2.17-e19-prof/lib/libMultiProc.so
bin/PandoraInterface: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_18_04b/Linux64bit+3.10-2.17-e19-prof/lib/libEve.so
bin/PandoraInterface: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_18_04b/Linux64bit+3.10-2.17-e19-prof/lib/libGeom.so
bin/PandoraInterface: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_18_04b/Linux64bit+3.10-2.17-e19-prof/lib/libRGL.so
bin/PandoraInterface: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_18_04b/Linux64bit+3.10-2.17-e19-prof/lib/libEG.so
bin/PandoraInterface: CMakeFiles/PandoraInterface.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable bin/PandoraInterface"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/PandoraInterface.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/PandoraInterface.dir/build: bin/PandoraInterface
.PHONY : CMakeFiles/PandoraInterface.dir/build

CMakeFiles/PandoraInterface.dir/requires: CMakeFiles/PandoraInterface.dir/test/PandoraInterface.cxx.o.requires
.PHONY : CMakeFiles/PandoraInterface.dir/requires

CMakeFiles/PandoraInterface.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/PandoraInterface.dir/cmake_clean.cmake
.PHONY : CMakeFiles/PandoraInterface.dir/clean

CMakeFiles/PandoraInterface.dir/depend:
	cd /var/clus/usera/afm67/2020/May/CosmicJorisMay2020/LArReco/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /var/clus/usera/afm67/2020/May/CosmicJorisMay2020/LArReco /var/clus/usera/afm67/2020/May/CosmicJorisMay2020/LArReco /var/clus/usera/afm67/2020/May/CosmicJorisMay2020/LArReco/build /var/clus/usera/afm67/2020/May/CosmicJorisMay2020/LArReco/build /var/clus/usera/afm67/2020/May/CosmicJorisMay2020/LArReco/build/CMakeFiles/PandoraInterface.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/PandoraInterface.dir/depend

