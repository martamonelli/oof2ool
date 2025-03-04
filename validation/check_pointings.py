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


file_path = base_path + '/030rpm_375/pnt_0A_00000.fits'
hdul = fits.open(file_path)
pnt = hdul[1].data

theta_030_375 = pnt['THETA']


file_path = base_path + '/100rpm_260/pnt_0A_00000.fits'
hdul = fits.open(file_path)
pnt = hdul[1].data

theta_100_260 = pnt['THETA']

nsamp = 2000

plt.plot(theta_005_450[:nsamp], label='0.05rpm and 45.0deg')
plt.plot(theta_030_375[:nsamp], label='0.30rpm and 37.5deg')
plt.plot(theta_100_260[:nsamp], label='1.00rpm and 26.0deg')
plt.legend()
plt.xlabel('samples')
plt.ylabel('theta')
plt.savefig('thetas.png')

