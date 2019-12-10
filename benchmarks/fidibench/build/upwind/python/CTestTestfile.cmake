# CMake generated Testfile for 
# Source directory: /home/tolnaia/COMBS/benchmarks/fidibench/upwind/python
# Build directory: /home/tolnaia/COMBS/benchmarks/fidibench/build/upwind/python
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(upwindPython "/opt/moose/miniconda/bin/python" "/home/tolnaia/COMBS/benchmarks/fidibench/upwind/python/upwind.py" "32" "10")
set_tests_properties(upwindPython PROPERTIES  PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1\\.0|0\\.999]")
add_test(upwindPython3 "/opt/moose/miniconda/bin/python" "/home/tolnaia/COMBS/benchmarks/fidibench/upwind/python/upwind3.py" "32" "10")
set_tests_properties(upwindPython3 PROPERTIES  PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1\\.0|0\\.999]")
add_test(upwindPythonMpi1 "/opt/moose/mpich-3.3/gcc-9.2.0/bin/mpiexec" "-n" "1" "/opt/moose/miniconda/bin/python" "/home/tolnaia/COMBS/benchmarks/fidibench/upwind/python/upwindMPI.py" "32" "10")
set_tests_properties(upwindPythonMpi1 PROPERTIES  PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1\\.0|0\\.999]")
add_test(upwindPythonMpiN "/opt/moose/mpich-3.3/gcc-9.2.0/bin/mpiexec" "-n" "16" "/opt/moose/miniconda/bin/python" "/home/tolnaia/COMBS/benchmarks/fidibench/upwind/python/upwindMPI.py" "32" "10")
set_tests_properties(upwindPythonMpiN PROPERTIES  PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1\\.0|0\\.999]")
