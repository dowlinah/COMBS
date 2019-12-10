# CMake generated Testfile for 
# Source directory: /home/tolnaia/COMBS/benchmarks/fidibench/upwind/fortran
# Build directory: /home/tolnaia/COMBS/benchmarks/fidibench/build/upwind/fortran
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(upwindFortran1 "/home/tolnaia/COMBS/benchmarks/fidibench/build/upwind/fortran/upwindFortran" "32" "10")
set_tests_properties(upwindFortran1 PROPERTIES  ENVIRONMENT "OMP_NUM_THREADS=1" PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1|0\\.999]")
add_test(upwindFortran2 "/home/tolnaia/COMBS/benchmarks/fidibench/build/upwind/fortran/upwindFortran" "32" "10")
set_tests_properties(upwindFortran2 PROPERTIES  ENVIRONMENT "OMP_NUM_THREADS=2" PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1|0\\.999]")
add_test(upwindFortran4 "/home/tolnaia/COMBS/benchmarks/fidibench/build/upwind/fortran/upwindFortran" "32" "10")
set_tests_properties(upwindFortran4 PROPERTIES  ENVIRONMENT "OMP_NUM_THREADS=4" PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1|0\\.999]")
add_test(upwindFortran8 "/home/tolnaia/COMBS/benchmarks/fidibench/build/upwind/fortran/upwindFortran" "32" "10")
set_tests_properties(upwindFortran8 PROPERTIES  ENVIRONMENT "OMP_NUM_THREADS=8" PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1|0\\.999]")
add_test(upwindF03 "/home/tolnaia/COMBS/benchmarks/fidibench/build/upwind/fortran/upwindF03" "32" "10")
set_tests_properties(upwindF03 PROPERTIES  PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1|0\\.999]")
add_test(upwindF08 "/home/tolnaia/COMBS/benchmarks/fidibench/build/upwind/fortran/upwindF08" "32" "10")
set_tests_properties(upwindF08 PROPERTIES  PASS_REGULAR_EXPRESSION "[Cc]heck sum:[ ]*[1|0\\.999]")
