#!/bin/bash -l
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_std
#SBATCH --qos=cm4_std
#SBATCH -t 06:00:00
#SBATCH -J test
#SBATCH -o logs/submit_cm4.log
#SBATCH --export=ALL
#SBATCH --mail-type=ALL
#SBATCH --mail-user=marta.monelli@ipmu.jp

eval "$(conda shell.bash hook)"
conda activate lbs

cd $SLURM_SUBMIT_DIR
export OMP_NUM_THREADS=3

### ARGUMENTS                              spin_rate_rpm  realization    net_ukrts  fknee_mhz  alpha  days  ptg_only  tod_only
for realization in {0..49}; do
  srun python simulating_TOD_noise-only.py 0.05           $realization   50.        50.        2.     1.    
done
