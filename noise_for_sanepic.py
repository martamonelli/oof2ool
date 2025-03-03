import matplotlib
matplotlib.use('Agg')
import numpy as np
import healpy as hp
import scipy as sc
import matplotlib.pyplot as plt
import math as mt
import array
import csv
import getopt, sys
import argparse

fknee_mhz = float(sys.argv[1])
alpha = float(sys.argv[2])
days = float(sys.argv[3])

sampling_rate_hz = 19.1

fknee_hz = fknee_mhz*1e-3

samples_in_a_day = sampling_rate_hz * 3600 * 24
samples_in_a_year = samples_in_a_day * 365
samples_in_a_chunk = int(2 ** np.ceil(np.log2(days * samples_in_a_day))) #FIXME: used to force circularity, to be relaxed
nchunks = int(np.ceil(samples_in_a_year/samples_in_a_chunk))
duration_d = nchunks*samples_in_a_chunk/samples_in_a_day

container_path = '/dss/dssfs02/lwp-dss-0001/pn36hu/pn36hu-dss-0000/oof2ool/'

bolonames = ['all_together']

ff = np.fft.fftfreq(samples_in_a_chunk)*sampling_rate_hz
ff[0] = ff[1]
sp = (fknee_hz/abs(ff))**alpha + 1
isp = 1.0/abs(sp)

for detname in bolonames:
    isp.tofile(container_path+'PROCESSED/for_sanepic/iSpf_'+str(int(days)).zfill(2)+'days_'+str(int(fknee_mhz)).zfill(3)+"mHz_"+str(int(alpha*10)).zfill(2)+'_'+detname)  
