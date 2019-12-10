#!/bin/bash

#SBATCH -J upwindMpiCxx_16
#SBATCH -A nesi99999
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntask=20

exe="/home/tolnaia/COMBS/benchmarks/fidibench/build/upwind/cxx/upwindMpiCxx"
srun $exe -numCells 200 -numSteps 4 


