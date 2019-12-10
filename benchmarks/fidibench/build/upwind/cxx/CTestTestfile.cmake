# CMake generated Testfile for 
# Source directory: /home/tolnaia/COMBS/benchmarks/fidibench/upwind/cxx
# Build directory: /home/tolnaia/COMBS/benchmarks/fidibench/build/upwind/cxx
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(upwindCxx1 "/home/tolnaia/COMBS/benchmarks/fidibench/build/upwind/cxx/upwindCxx" "32" "10")
set_tests_properties(upwindCxx1 PROPERTIES  ENVIRONMENT "OMP_NUM_THREADS=1" PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1|0\\.999]")
add_test(upwindCxx2 "/home/tolnaia/COMBS/benchmarks/fidibench/build/upwind/cxx/upwindCxx" "32" "10")
set_tests_properties(upwindCxx2 PROPERTIES  ENVIRONMENT "OMP_NUM_THREADS=2" PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1|0\\.999]")
add_test(upwindCxx4 "/home/tolnaia/COMBS/benchmarks/fidibench/build/upwind/cxx/upwindCxx" "32" "10")
set_tests_properties(upwindCxx4 PROPERTIES  ENVIRONMENT "OMP_NUM_THREADS=4" PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1|0\\.999]")
add_test(upwindCxx8 "/home/tolnaia/COMBS/benchmarks/fidibench/build/upwind/cxx/upwindCxx" "32" "10")
set_tests_properties(upwindCxx8 PROPERTIES  ENVIRONMENT "OMP_NUM_THREADS=8" PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1|0\\.999]")
add_test(upwindCxx "/home/tolnaia/COMBS/benchmarks/fidibench/build/upwind/cxx/upwindCxx" "32" "10")
set_tests_properties(upwindCxx PROPERTIES  PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1|0\\.999]")
add_test(upwindAccCxx "/home/tolnaia/COMBS/benchmarks/fidibench/build/upwind/cxx/upwindAccCxx" "32" "10")
set_tests_properties(upwindAccCxx PROPERTIES  PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1|0\\.999]")
add_test(upwindAcc2Cxx "/home/tolnaia/COMBS/benchmarks/fidibench/build/upwind/cxx/upwindAcc2Cxx" "32" "10")
set_tests_properties(upwindAcc2Cxx PROPERTIES  PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1|0\\.999]")
