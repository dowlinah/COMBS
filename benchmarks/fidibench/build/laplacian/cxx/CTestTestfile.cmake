# CMake generated Testfile for 
# Source directory: /home/tolnaia/COMBS/benchmarks/fidibench/laplacian/cxx
# Build directory: /home/tolnaia/COMBS/benchmarks/fidibench/build/laplacian/cxx
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(testStencil2d "/opt/moose/mpich-3.3/gcc-9.2.0/bin/mpiexec" "-n" "8" "./testStencil2d" "-numCells" "32")
add_test(laplacian2D "/opt/moose/mpich-3.3/gcc-9.2.0/bin/mpiexec" "-n" "8" "./laplacian" "-numDims" "2" "-numCells" "32")
set_tests_properties(laplacian2D PROPERTIES  PASS_REGULAR_EXPRESSION "[Cc]heck sums:[ ]*")
add_test(laplacian3D "/opt/moose/mpich-3.3/gcc-9.2.0/bin/mpiexec" "-n" "8" "./laplacian" "-numDims" "3" "-numCells" "32")
set_tests_properties(laplacian3D PROPERTIES  PASS_REGULAR_EXPRESSION "[Cc]heck sums:[ ]*")
