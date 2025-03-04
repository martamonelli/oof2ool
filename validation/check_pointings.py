import numpy as np
import healpy as hp
from scipy.optimize import curve_fit

from astropy.io import fits
from matplotlib import pyplot as plt

'''
plt.rcParams.update({
    'text.usetex': True,
    'font.family': 'serif'
})
'''

base_path = '/dss/dssfs02/lwp-dss-0001/pn36hu/pn36hu-dss-0000/oof2ool/pointings'

file_path = base_path + '/005rpm_450/pnt_0A_00000.fits'
hdul = fits.open(file_path)
pnt = hdul[1].data

theta_005_450 = pnt['THETA']


file_path = base_path + '/030rpm_450/pnt_0A_00000.fits'
hdul = fits.open(file_path)
pnt = hdul[1].data

theta_030_450 = pnt['THETA']


file_path = base_path + '/030rpm_375/pnt_0A_00000.fits'
hdul = fits.open(file_path)
pnt = hdul[1].data

theta_030_375 = pnt['THETA']

nsamp_030 = 500
ratio = int(30/5)

samps_030 = np.arange(nsamp_030)
samps_005 = np.arange(nsamp_030*ratio)/ratio

plt.plot(samps_005, theta_005_450[:nsamp_030*ratio], label='theta for 0.05rpm and 45deg, stretched')
plt.plot(samps_030, theta_030_450[:nsamp_030], label='theta for 0.30rpm and 45deg')
plt.legend()
plt.savefig('thetas.png')

