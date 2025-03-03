import numpy as np
import scipy as sc
import math as mt
from astropy.io import fits
import array
import csv
import sys
import struct

spin_rate_rpm = float(sys.argv[1])
realizations = int(sys.argv[2])
net_ukrts = float(sys.argv[3])
fknee_mhz = float(sys.argv[4])
alpha = float(sys.argv[5])
days = float(sys.argv[6])

sampling_rate_hz = 19.1

samples_in_a_day = sampling_rate_hz * 3600 * 24
samples_in_a_year = samples_in_a_day * 365
samples_in_a_chunk = int(2 ** np.ceil(np.log2(days * samples_in_a_day))) #FIXME: used to force circularity, to be relaxed
nchunks = int(np.ceil(samples_in_a_year/samples_in_a_chunk))
duration_d = nchunks*samples_in_a_chunk/samples_in_a_day

### redefining Guillaume's variables
nn = nchunks
interv = samples_in_a_chunk
fk = fknee_mhz
simnum = str(realizations).zfill(2)

scratch_path = '/dss/lxclscratch/05/di38wul/oof2ool'
container_path = '/dss/dssfs02/lwp-dss-0001/pn36hu/pn36hu-dss-0000/oof2ool'
termin = 'all_together'

bolos = ['0A','0B','1A','1B']

print('realization '+simnum)
    
dirdn = container_path+'/'+str(int(days)).zfill(2)+'days/'+str(int(fknee_mhz)).zfill(3)+"mHz_"+str(int(alpha*10)).zfill(2)+'/'+simnum+'/'
dirout = scratch_path+'/'+str(int(days)).zfill(2)+'days/'+str(int(fknee_mhz)).zfill(3)+"mHz_"+str(int(alpha*10)).zfill(2)

#dirdn = '/dss/dsshome1/05/di38wul/SCRATCH/out/cmb/' #FIXME: remove 
#dirout = container_path                             #FIXME: remove

filenameoutall_n = dirout+'/'+simnum+'_'+termin+'.bi'
fdataall_n = open(filenameoutall_n,'wb')
    
iiall = 0
isampall = 0
for bolo in bolos:
    fname = "tod_"+bolo
        
    iilist = np.arange(0,nn*interv,interv)
    #print(iilist)
        
    isamp=0
    iii=0
    for ii in iilist:
        filenamedata_n = dirdn+fname+'_'+str(iii).zfill(5)+'.fits'
        print(filenamedata_n)
        hduln = fits.open(filenamedata_n)
        data_n = hduln[1].data['TOD']

        fdataall_n.seek(isampall*8)
        sdat_n = struct.pack('d'*len(data_n),*data_n)
        fdataall_n.write(sdat_n)
            
        isamp += np.size(data_n)
        isampall += np.size(data_n)
        iii+=1
        iiall+=1
        #print(ii)
        
fdataall_n.close()
