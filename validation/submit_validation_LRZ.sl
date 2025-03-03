#!/bin/bash -l
#SBATCH -D . # work dir == submit dir
#SBATCH --ntasks=1 # you need 1 CPU, right?
#SBATCH --clusters=serial
#SBATCH --partition=serial_std # or serial_long ?
#SBATCH -t 24:00:00 # no! serial_std is ok.
#SBATCH -J test
#SBATCH -o LOG/validation.log
#SBATCH --export=none # really better ... so we oblige that to avoid side effects
#SBATCH --mail-type=ALL
#SBATCH --mail-user=marta.monelli@ipmu.jp

module load slurm_setup

# please setup here your environment; load modules; activate conda environments; ...
python SCRIPTS/validation.py 
