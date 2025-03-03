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
#SBATCH --time=0:19:00

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

for simnum in 00
do
mpirun -np $nproc /dss/dsshome1/05/di38wul/sanepic/Sanepic \
### dirfile:
-F /dss/lxclscratch/05/di38wul/oof2ool/${days_str}/${fknee_str}mhz_${alpha_str}/${simnum}_ \
### bextension:
-B ${simnum} \
### filesegnum:
-f /dss/dssfs02/lwp-dss-0001/pn36hu/pn36hu-dss-0000/oof2ool/PROCESSED/for_sanepic/samplist_${days_str}_${det_str}.bi \
### ndet:
-d 1 \
### bolos:
-C ${det_str} \
### termin_out:
-e _${simnum}_${det_str} \
### noiseSppreffile:
-k /dss/dssfs02/lwp-dss-0001/pn36hu/pn36hu-dss-0000/oof2ool/PROCESSED/for_sanepic/iSpf_${days_str}_${fknee_str}_${alpha_str}\
### nssp:
-u $seglength \
### outdir:
-O /dss/dssfs02/lwp-dss-0001/pn36hu/pn36hu-dss-0000/oof2ool/PROCESSED/${days_str}_${fknee_str}_${alpha_str}_${spin_str} \
### polar:
-p 1 \
### ns:
-l $seglength \
### nsegtot:
-n $nnn \
### nside:
-N 256 \
### filedetinfo:
-X /dss/dssfs02/lwp-dss-0001/pn36hu/pn36hu-dss-0000/oof2ool/PROCESSED/for_sanepic/BoloFile_00.txt \
### dirfilep:
-Z /dss/lxclscratch/05/di38wul/oof2ool/pointings/${spin_str}/${det_str}_
done
