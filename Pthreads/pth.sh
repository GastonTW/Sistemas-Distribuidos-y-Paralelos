#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --partition=Blade
#SBATCH -o /nethome/sdyp10/Pthreads/output.txt
#SBATCH -e /nethome/sdyp10/Pthreads/errors.txt
./pthread_g10 $1 $2 $3 $4