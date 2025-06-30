#!/bin/bash

# Request an hour of runtime:
#SBATCH --time=23:59:00
#SBATCH -p batch
# Controls the minimum/maximum number of nodes allocated to the job
#SBATCH -N 1
#SBATCH -c 16

# Default resources are 1 core with 2.8GB of memory.

# Use more memory (4GB):
#SBATCH --mem=32G

# Specify a job name:
#SBATCH -J matlabjob

# Specify an output file
#SBATCH -o job%j.out
#SBATCH -e job%j.err

module load matlab/R2019a

# Run a matlab function called 'foo.m' in the same directory as this batch script.
matlab-threaded -nodisplay -nojvm -r "s_gen_Nt_resolution, exit"
