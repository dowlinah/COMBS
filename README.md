# COMBS: Clarkson Open-Source Multi-Physics Benchmark Suite

COMBS is an open-source multi-physics benchmark suite. This
suite draws on an assortment of open-source high performance
computing benchmarks, and aggregates them into one unified
benchmark suite for testing an assortment of features of
a system's performance.

This system has been tested and is designed for use on an
Ubuntu 18.04 system. It may be possible to use it on other
operating systems, but that is currently untested and
unsupported.

## Usage

There are three main scripts for using the suite. These scripts do
what their names imply:

- deps.sh
  - Installs the dependencies for the benchmark suite.
- build.sh
  - Builds each benchmark from their source code to prepare for running 
  it and gathering the data from it
- run.sh
  - Runs each benchmark and stores the output data in each benchmark's 
  directory before aggregating the data in the top level directory as a 
  CSV file and a TeX table.
- clean.sh
  - Cleans each benchmark directory, removing build files.

So, to run the benchmarks, first deps.sh should be used to
install the dependencies for the suite. Then, build.sh mush
be ran successfully, before run.sh can be ran. This script
can take a fairly long time, as it has to run every
benchmark and gather the data from running it.

Once run.sh completes, the data will be in the top directory as a csv file and a tex file as well. After running the benchmarks, the clean.sh
script may be used to remove the build files to save room on the system.
