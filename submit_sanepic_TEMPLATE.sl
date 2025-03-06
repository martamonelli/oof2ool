#!/bin/bash
#SBATCH -o ./myjob.%j.%N.out
#SBATCH -D ./
#SBATCH --clusters=
#SBATCH --partition=
#SBATCH --qos=
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=
### --- Maximum number of CPU cores is 112 for cm4 - Use CM4 resources carefully and efficiently ! ---
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --mail-type=ALL
#SBATCH --mail-user=marta.monelli@ipmu.jp
#SBATCH --time=8:00:00

eval "$(conda shell.bash hook)"
conda activate lbs

cd $SLURM_SUBMIT_DIR
export OMP_NUM_THREADS=1

module load intel/2021.4.0 intel-mkl/2021.4.0 cfitsio/4.3.0-intel23 intel-mpi/2021.9.0 fftw/3.3.10-intel23-impi-openmp

nproc=
nnn=

days_str=
spin_str=

fknee_str=
alpha_str=

det_str=all_together
seglength=

for simnum in {00, 01, 02, 03, 04, 05, 06, 07, 08, 09, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49}
do
mpirun -np $nproc /dss/dsshome1/05/di38wul/sanepic/Sanepic -F /dss/lxclscratch/05/di38wul/oof2ool/${days_str}/${fknee_str}_${alpha_str} -B ${simnum} -f /dss/dssfs02/lwp-dss-0001/pn36hu/pn36hu-dss-0000/oof2ool/PROCESSED/for_sanepic/samplist_${days_str}_${det_str}.bi -d 1 -C ${det_str} -e _${simnum}_${det_str} -k /dss/dssfs02/lwp-dss-0001/pn36hu/pn36hu-dss-0000/oof2ool/PROCESSED/for_sanepic/iSpf_${days_str}_${fknee_str}_${alpha_str}_ -u $seglength -O /dss/dssfs02/lwp-dss-0001/pn36hu/pn36hu-dss-0000/oof2ool/PROCESSED/${days_str}_${fknee_str}_${alpha_str}_${spin_str}/ -p 1 -l $seglength -n $nnn -N 256 -X /dss/dssfs02/lwp-dss-0001/pn36hu/pn36hu-dss-0000/oof2ool/PROCESSED/for_sanepic/BoloFile_00.txt -Z /dss/lxclscratch/05/di38wul/oof2ool/pointings/${spin_str}/${det_str}_
done
