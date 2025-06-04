#!/bin/bash
#SBATCH -N 2
#SBATCH --exclusive
#SBATCH --partition=Blade
#SBATCH --tasks-per-node=1
#SBATCH -o /nethome/sdyp10/MPI-Pthreads/output.txt
#SBATCH -e /nethome/sdyp10/MPI-Pthreads/errors.txt
mpirun --bind-to none mpi_g10 $1 $2 $3 $4