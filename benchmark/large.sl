#!/bin/bash -l
#SBATCH -p regular
#SBATCH -N 768
#SBATCH -t 00:01:00
#SBATCH -J my_job
#SBATCH -o my_job.o%j
#SBATCH -A mpccc

export OMP_NUM_THREADS=6
export KMP_AFFINITY=compact,granularity=core,1

srun -n 3072 -c 6 ../src/mini_dft -nbgrp 96 -in large.in > large.out
