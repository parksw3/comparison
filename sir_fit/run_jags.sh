#!/bin/bash
#SBATCH --account=def-bolker # replace this with your own account
#SBATCH --ntasks=1               # number of MPI processes
#SBATCH --mem-per-cpu=2048M      # memory; default unit is megabytes
#SBATCH --time=1-00:00           # time (DD-HH:MM)
fn=jags_fit_$I_FOR
mpirun -np 1 R CMD BATCH --vanilla "--args $I_FOR" jags_fit.R $fn.Rout