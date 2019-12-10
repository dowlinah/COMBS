MOOSE
=====

[![Build status](https://www.moosebuild.org/idaholab/moose/master/branch_status.svg)](https://www.moosebuild.org/repo/idaholab/moose/)

The Multiphysics Object-Oriented Simulation Environment (MOOSE) is a finite-element, multiphysics framework primarily developed by [Idaho National Laboratory](http://www.inl.gov). It provides a high-level interface to some of the most sophisticated [nonlinear solver technology](http://www.mcs.anl.gov/petsc/) on the planet. MOOSE presents a straightforward API that aligns well with the real-world problems scientists and engineers need to tackle. Every detail about how an engineer interacts with MOOSE has been thought through, from the installation process through running your simulation on state of the art supercomputers, the MOOSE system will accelerate your research.

Some of the capability at your fingertips:

* Fully-coupled, fully-implicit multiphysics solver
* Dimension independent physics
* Automatically parallel (largest runs >100,000 CPU cores!)
* Modular development simplifies code reuse
* Built-in mesh adaptivity
* Continuous and Discontinuous Galerkin (DG) (at the same time!)
* Intuitive parallel multiscale solves (see videos below)
* Dimension agnostic, parallel geometric search (for contact related applications)
* Flexible, plugable graphical user interface
* ~30 plugable interfaces allow specialization of every part of the solve

SETUP libmesh
=============
Build Instructions For libmesh
The default is to build libmesh "out of tree," i.e. within a separate build directory, rather than in the source tree itself. This simplifies the process of having multiple, independently-configured builds.

cd to location of libmesh clone
(Only if using a git clone) git submodule update --init
mkdir build
./configure
make
make check (optional, runs the example programs and unit tests when possible)
make install

SETUP MOOSE
=============

Prerequisites
--------------

sudo apt-get install 
  build-essential \
  gfortran \
  tcl \
  git \
  m4 \
  freeglut3 \
  doxygen \
  libblas-dev \
  liblapack-dev \
  libx11-dev \
  libnuma-dev \
  libcurl4-gnutls-dev \
  zlib1g-dev \
  libhwloc-dev \
  libxml2-dev \
  libpng-dev \
  pkg-config \
  liblzma-dev
  
  Environment
  --------------
  In instructions change ~/projects to pathway to this benchmark
  https://mooseframework.inl.gov/old/wiki/BasicManualInstallation/Linux/
  
  
  Modify your Bash Profile
  ------------------------
  
  echo "module load moose-dev-gcc" >> ~/.bash_profile
  
  Compile libMesh
  --------------
  
  ./scripts/update_and_rebuild_libmesh.sh
  
  Compile and Test MOOSE
  ----------------------
  
  cd test
  
  make -j Num_Threads
  
  ./run_tests - Num_Threads
  

More Information
================

**For more information, including installation instructions, please see the official website: [http://mooseframework.org](http://mooseframework.org)**

Contributing
============

For information on how to contribute code changes to MOOSE please [see the mooseframework.org wiki](http://mooseframework.org/wiki/Contributing/).
