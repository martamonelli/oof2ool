'''

This scripts checks 

'''

import numpy as np 

import sys
import os
import re

spin_rate_rpm = float(sys.argv[1])
fknee_mhz = float(sys.argv[2])
alpha = float(sys.argv[3])
days = float(sys.argv[4])

if spin_rate_rpm == 0.05:
    spin_sun_angle_deg = 45.
elif spin_rate_rpm == 0.3:
    spin_sun_angle_deg = 37.5
elif spin_rate_rpm == 1.:
    spin_sun_angle_deg = 26.
else:
    print('spin_rate_rpm should be 0.05, 0.3 or 1.')
    quit()

if days == 20:
    clusters_str = '#SBATCH --clusters=inter'
    partition_str = '#SBATCH --partition=cm4_inter'
    qos_str = '###'
elif days == 1:
    clusters_str = '#SBATCH --clusters=cm4'
    partition_str = '#SBATCH --partition=cm4_standard'
    qos_str = '#SBATCH --qos=cm4_standard'
else:
    print('days should be 20 or 1')
    quit()

#hard-coded params
sampling_rate_hz = 19.1
ndet = 4

samples_in_a_day = sampling_rate_hz * 3600 * 24
samples_in_a_year = samples_in_a_day * 365
samples_in_a_chunk = int(2 ** np.ceil(np.log2(days * samples_in_a_day))) #FIXME: used to force circularity, to be relaxed
nchunks = int(np.ceil(samples_in_a_year/samples_in_a_chunk))
duration_d = nchunks*samples_in_a_chunk/samples_in_a_day

scratch_path = '/dss/lxclscratch/05/di38wul/oof2ool'
container_path = '/dss/dssfs02/lwp-dss-0001/pn36hu/pn36hu-dss-0000/oof2ool/'

### checking if the stacked pointings exist

### checking if the stacked TOD exist

### checking if the [...] exists

### checking if the output directory exists
out_path = container_path+'PROCESSED/'+str(int(days)).zfill(2)+'days_'+str(int(fknee_mhz)).zfill(3)+"mHz_"+str(int(alpha*10)).zfill(2)+'_'+str(int(spin_rate_rpm*100)).zfill(3)+'rpm_'+str(int(10*spin_sun_angle_deg))

os.system('cp submit_sanepic_TEMPLATE.sl '+out_path+'/submit_sanepic.sl')
file_path = out_path+'/submit_sanepic.sl'

with open(file_path, "r") as file:
    lines = file.readlines()

# Look for the line containing "base_first" and replace the number
for i, line in enumerate(lines):
    if re.search(r"#SBATCH --clusters=.*", line):
        lines[i] = (clusters_str+"\n")
    if re.search(r"#SBATCH --partition=.*", line):
        lines[i] = (partition_str+"\n")
    if re.search(r"#SBATCH --qos=.*", line):
        lines[i] = (qos_str+"\n")
    if re.search(r"#SBATCH --ntasks-per-node=.*", line):
        lines[i] = ("#SBATCH --ntasks-per-node="+str(nchunks*ndet+1)+"\n")
    if re.search(r"nproc=.*", line):
        lines[i] = ("nproc="+str(nchunks*ndet+1)+"\n")
    if re.search(r"nnn=.*", line):
        lines[i] = ("nnn="+str(nchunks*ndet)+"\n")
    if re.search(r"seglength=.*", line):
        lines[i] = ("seglength="+str(samples_in_a_chunk)+"\n")
    if re.search(r"fknee_str=.*", line):
        lines[i] = ("fknee_str="+str(int(fknee_mhz)).zfill(3)+'mHz\n')
    if re.search(r"alpha_str=.*", line):
        lines[i] = ("alpha_str="+str(int(alpha*10))+"\n")
    if re.search(r"spin_str=.*", line):
        lines[i] = ("spin_str="+str(int(spin_rate_rpm*100)).zfill(3)+'rpm_'+str(int(10*spin_sun_angle_deg))+"\n")
    if re.search(r"days_str=.*", line):
        lines[i] = ("days_str="+str(int(days)).zfill(2)+"days\n")

# Write the updated content back to the file
with open(file_path, "w") as file:
    file.writelines(lines)

