#!/bin/bash

#SBATCH --job-name=test       ###Job Name
#SBATCH --partition=kamiak,cas,clark    ###Partition on which to run
#SBATCH --nodes=1             ###Number of nodes to use
#SBATCH --ntasks-per-node=20   ###Number of tasks per node (aka MPI processes)
#SBATCH --cpus-per-task=1    ###Number of cpus per task (aka OpenMP threads)
#SBATCH --time=7-00:00:00
module load gaussian          ###Load gaussian module on Kamiak

finit=$1                      ###The base name of your input file  
fend='.com'                     
foutend='.out'

export GAUSS_SCRDIR="$(mkworkspace -q)"  ###Creates a workspace

g09 < ${finit}${fend} > ${finit}${foutend}                                     
