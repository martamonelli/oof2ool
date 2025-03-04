#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=18
#SBATCH --cpus-per-task=1
#SBATCH --clusters=inter
#SBATCH --partition=cm4_inter
#SBATCH -t 01:30:00
#SBATCH -J test
#SBATCH -o logs/submit_inter.log
#SBATCH --export=ALL
#SBATCH --mail-type=ALL
#SBATCH --mail-user=marta.monelli@ipmu.jp

eval "$(conda shell.bash hook)"
conda activate lbs

cd $SLURM_SUBMIT_DIR
export OMP_NUM_THREADS=1

### ARGUMENTS                                   spin_rate_rpm  realization    net_ukrts  fknee_mhz  alpha  days  ptg_only  tod_only
for realization in {0..00}; do
  srun python simulating_TOD_noise-only.py      1.           $realization   50.        100.       2.     20.  
done
