#!/bin/bash
#SBATCH --job-name=upwindAcc2
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=4096
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=gpu
#SBATCH --hint=nomultithread    # don't use hyperthreading

exe="/home/tolnaia/COMBS/benchmarks/fidibench/build/upwind/cxx/upwindAcc2Cxx"
time srun $exe -numCells 512 -numSteps 10

