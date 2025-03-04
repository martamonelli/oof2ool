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

###

file_path = base_path + '/005rpm_450/pnt_0A_00000.fits'
hdul = fits.open(file_path)
pnt = hdul[1].data

theta_005_450 = pnt['THETA']
phi_005_450 = pnt['PHI']
psi_005_450 = pnt['PSI']

###

file_path = base_path + '/030rpm_375/pnt_0A_00000.fits'
hdul = fits.open(file_path)
pnt = hdul[1].data

theta_030_375 = pnt['THETA']
phi_030_375 = pnt['PHI']
psi_030_375 = pnt['PSI']

###

file_path = base_path + '/100rpm_260/pnt_0A_00000.fits'
hdul = fits.open(file_path)
pnt = hdul[1].data

theta_100_260 = pnt['THETA']
phi_100_260 = pnt['PHI']
psi_100_260 = pnt['PSI']

###

nsamp = 50000

plt.plot(theta_005_450[:nsamp], label='0.05rpm and 45.0deg')
plt.plot(theta_030_375[:nsamp], label='0.30rpm and 37.5deg')
plt.plot(theta_100_260[:nsamp], label='1.00rpm and 26.0deg')
plt.legend()
plt.xlabel('samples')
plt.ylabel('theta')
plt.savefig('thetas.png')
plt.clf()

plt.plot(phi_005_450[:nsamp], label='0.05rpm and 45.0deg')
plt.plot(phi_030_375[:nsamp], label='0.30rpm and 37.5deg')
plt.plot(phi_100_260[:nsamp], label='1.00rpm and 26.0deg')
plt.legend()
plt.xlabel('samples')
plt.ylabel('phi')
plt.savefig('phis.png')
plt.clf()

plt.plot(psi_005_450[:nsamp], label='0.05rpm and 45.0deg')
plt.plot(psi_030_375[:nsamp], label='0.30rpm and 37.5deg')
plt.plot(psi_100_260[:nsamp], label='1.00rpm and 26.0deg')
plt.legend()
plt.xlabel('samples')
plt.ylabel('psi')
plt.savefig('psis.png')
plt.clf()

###

pix_005_450 = hp.ang2pix(256,theta_005_450[:nsamp],phi_005_450[:nsamp])
pix_030_375 = hp.ang2pix(256,theta_030_375[:nsamp],phi_030_375[:nsamp])
pix_100_260 = hp.ang2pix(256,theta_100_260[:nsamp],phi_100_260[:nsamp])

plt.plot(pix_005_450[:nsamp], label='0.05rpm and 45.0deg')
plt.plot(pix_030_375[:nsamp], label='0.30rpm and 37.5deg')
plt.plot(pix_100_260[:nsamp], label='1.00rpm and 26.0deg')
plt.legend()
plt.xlabel('samples')
plt.ylabel('pixel (nside=256)')
plt.savefig('pixs.png')
plt.clf()

