'''

script to reformat the pointings to be read from sanepic
it takes spin_rate_rpm and days in input (#TODO: allow arbitrary values for spin_sun_angle)

'''

import numpy as np
import scipy as sc
import math as mt
from astropy.io import fits
import array
import csv
import sys
import struct
import os

spin_rate_rpm = float(sys.argv[1])
days = float(sys.argv[2])

if spin_rate_rpm == 0.05:
    spin_sun_angle_deg = 45.
elif spin_rate_rpm == 0.3:
    spin_sun_angle_deg = 37.5
elif spin_rate_rpm == 1.:
    spin_sun_angle_deg = 26.
else:
    print('spin_rate_rpm should be 0.05, 0.3 or 1.')
    quit()

sampling_rate_hz = 19.1

samples_in_a_day = sampling_rate_hz * 3600 * 24
samples_in_a_year = samples_in_a_day * 365
samples_in_a_chunk = int(2 ** np.ceil(np.log2(days * samples_in_a_day))) #FIXME: used to force circularity, to be relaxed
nchunks = int(np.ceil(samples_in_a_year/samples_in_a_chunk))
duration_d = nchunks*samples_in_a_chunk/samples_in_a_day

### redefining Guillaume's variables
#nn = 18
#interv = 8388608*4
nn = nchunks
interv = samples_in_a_chunk 

scratch_path = '/dss/lxclscratch/05/di38wul/oof2ool'
container_path = '/dss/dssfs02/lwp-dss-0001/pn36hu/pn36hu-dss-0000/oof2ool'

dird = container_path+'/pointings/'+str(int(spin_rate_rpm*100)).zfill(3)+'rpm_'+str(int(10*spin_sun_angle_deg))+'/'
dirout = scratch_path+'/pointings/'+str(int(spin_rate_rpm*100)).zfill(3)+'rpm_'+str(int(10*spin_sun_angle_deg))+'/'
termin = 'all_together'
fnames = ["pnt_0A", "pnt_0B", "pnt_1A", "pnt_1B"]
ndet = len(fnames)

### Guillaume's code
nnn=nn*np.size(fnames)
iilist = np.arange(0,nn*interv,interv)
dd1 = np.arange(0,nn*interv*np.size(fnames),interv)
dd2 = dd1 + interv - 1

dd = np.zeros((np.size(dd1)*2))
dd[0:2*nnn:2] = dd1
dd[1:2*nnn+11:2] = dd2
ddd = dd.astype(np.int64)

### my version
i_first = np.arange(nchunks*ndet)*samples_in_a_chunk
i_last = np.arange(1,nchunks*ndet+1)*samples_in_a_chunk - 1

test = np.empty(2*nchunks*ndet)
test[0::2] = i_first
test[1::2] = i_last
test = test.astype(np.int64)

test.tofile(container_path+'/PROCESSED/for_sanepic/samplist_'+str(int(days)).zfill(2)+'days_'+termin+'.bi')

filenameraout = dirout+termin+'_ra.bi'
filenamedecout = dirout+termin+'_dec.bi'
filenamepsiout = dirout+termin+'_psi.bi'
fra = open(filenameraout,'wb')
fdec = open(filenamedecout,'wb')
fpsi = open(filenamepsiout,'wb')

isampall=0
iiall=0
for fname in fnames:
    iii=0
    for ii in iilist:
        filenamedata = dird+fname+'_'+str(iii).zfill(5)+'.fits'
        print('going through '+filenamedata)
        hdul = fits.open(filenamedata)
        ptdata = hdul[1].data
        
        thetap = ptdata['THETA']
        #print(np.std(thetap))
        rap = ptdata['PHI']
        #print(np.std(rap))
        psip = ptdata['PSI']
        
        decp = np.pi/2.0 - thetap
        
        sra = struct.pack('d'*len(rap), *rap)
        sdec = struct.pack('d'*len(decp), *decp)
        spsi = struct.pack('d'*len(psip), *psip)
        fra.seek(isampall*8)
        fdec.seek(isampall*8)
        fpsi.seek(isampall*8)
        fra.write(sra)
        fdec.write(sdec)
        fpsi.write(spsi)
        
        #print(isampall,thetap[0])
        
        isampall += np.size(thetap)
        iiall+=1
        iii+=1
        #print(ii)

fra.close()
fdec.close()
fpsi.close()

psiall = np.fromfile(filenamepsiout, dtype=np.double)
#print(np.std(psiall[0:8388608]))
#print(np.std(psiall[8388608:8388608*2]))
#print(psiall[0:10])
#print(psiall[8388608-10:8388608+10])
