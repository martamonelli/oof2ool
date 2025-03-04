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

plt.plot(theta_005_450[:500], label='theta for 0.05rpm and 45.0deg')
plt.plot(theta_030_450[:500], label='theta for 0.30rpm and 45.0deg')
plt.plot(theta_030_375[:500], label='theta for 0.30rpm and 37.0deg')
plt.legend()
plt.savefig('thetas.png')

