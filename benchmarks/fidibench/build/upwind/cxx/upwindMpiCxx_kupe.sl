#!/bin/bash
#SBATCH --job-name=upwindMpiCxx
#SBATCH --time=00:10:00
#SBATCH --partition=NeSI
#SBATCH --account=nesi99999

exe="/home/tolnaia/COMBS/benchmarks/fidibench/build/upwind/cxx/upwindMpiCxx"
time srun $exe -numCells 800 -numSteps 2

