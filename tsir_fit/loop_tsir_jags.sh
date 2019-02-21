module load r/3.4.0
module load jags/4.3.0
module load openmpi/1.10.7
for ((i=0; i<=9; i++))
 do
 # Using environment variable I_FOR to communicate the case number to individual jobs:
 export I_FOR=$i
 sbatch run_tsir_jags.sh
 done