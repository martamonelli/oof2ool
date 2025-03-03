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

nreal = 50
nchunks = 18

base_path = '/dss/dssfs02/lwp-dss-0001/pn36hu/pn36hu-dss-0000/oof2ool/20days'

fsamp  = 19.1
sigma  = 50.*np.sqrt(fsamp)
f_knee = 0.01
slope  = 1.

################################
### reading Guillaume's iSpf
################################

iSpf_long = np.fromfile('/dss/dssfs02/lwp-dss-0001/pn36hu/pn36hu-dss-0000/oof2ool/SanepicMaps/iSpf0.01_19.1_index1.0_long__quadr1', np.float64)
iSpf_short = np.fromfile('/dss/dssfs02/lwp-dss-0001/pn36hu/pn36hu-dss-0000/oof2ool/SanepicMaps/iSpf0.01_19.1_index1.0_short__quadr1', np.float64)

################################
### pure white noise (hist)
################################

white_tod = np.array([])

for i in np.arange(nchunks):
    filepath = base_path + '/010mHz_10_SPLIT/00/tod_0A_' + str(i).zfill(5) + '.fits'
    hdul = fits.open(filepath)
    tod = hdul[1].data
    white_tod = np.concatenate((white_tod, tod['TOD']), axis=None)

nbins = 100

white_hist, bin_edges = np.histogram(white_tod, bins=nbins)
bin_centers = (bin_edges[1:] + bin_edges[:-1])/2

def zero_mean_gaussian(x,a,sigma):
    return a*np.exp(-x**2/(2*sigma**2))

popt, pcov = curve_fit(zero_mean_gaussian, bin_centers, white_hist)#, p0=[len(white_tod),200])
print('fitting the binned TOD with a Gaussian:')
print('sigma=%5.2f \pm %5.2f\n' % tuple([popt[1],pcov[1,1]**0.5])) 

plt.figure(figsize=(5,3.5))
plt.hist(white_tod, bins=nbins, color='lightgray')
plt.plot(bin_centers, (bin_edges[1]-bin_edges[0])*zero_mean_gaussian(bin_centers,len(white_tod)/(sigma*(2*np.pi)**0.5),sigma), color='black', label='expectation') 
plt.plot(bin_centers, zero_mean_gaussian(bin_centers, *popt), color='orangered', label='best fit')#label=r'$\sigma=%5.2f\pm%5.2f$' % tuple([popt[1],pcov[1,1]**0.5]))
plt.legend(loc='upper right')
plt.xlabel(r'binned TOD [$\mu$K]')
plt.ylabel(r'occurrencies')
plt.tight_layout()
plt.savefig('hist.pdf')
plt.clf()

################################
### pure white noise (PS)
################################

white_ps = np.empty((nreal,len(white_tod)))
oof_ps = np.empty((nreal,len(white_tod)))

for i in np.arange(nreal):    
    white_tod = np.array([])
    oof_tod = np.array([])

    for j in np.arange(nchunks):
        filepath = base_path + '/010mHz_10_SPLIT/' + str(i).zfill(2) + '/tod_0A_' + str(j).zfill(5) + '.fits'
        hdul = fits.open(filepath)
        
        tod = hdul[1].data
        white_tod = np.concatenate((white_tod, tod['TOD']), axis=None)
        
        tod = hdul[2].data
        oof_tod = np.concatenate((oof_tod, tod['TOD']), axis=None)

    white_ps[i] = np.abs(np.fft.fft(white_tod))**2/len(white_tod)
    oof_ps[i] = np.abs(np.fft.fft(white_tod+oof_tod))**2/len(oof_tod)

###preparing freqs array(s)
time_step = 1. / fsamp
freqs = np.fft.fftfreq(tod.size, time_step)
freqs_pos = freqs[:freqs.size//2]             #keeping only positive freqs

fmin = fsamp / len(tod)
fidx = np.argwhere(freqs_pos>=fmin)
freqs_sub = freqs_pos[fidx].flatten()         #keeping only sub-chunk information

###preparing power spectrum array(s)
white_ps_pos = white_ps[:,:freqs.size//2]     #keeping only positive freqs

oof_ps_pos = oof_ps[:,:freqs.size//2]         #keeping only positive freqs
oof_ps_sub = oof_ps_pos[:,fidx]               #keeping only sub-chunk information

iSpf_long_pos = iSpf_long[:freqs.size//2]     #keeping only positive freqs
iSpf_long_sub = iSpf_long_pos[fidx].flatten() #keeping only sub-chunk information
Spf_long = 1/iSpf_long_sub

def white_ps_theo(freqs, sigma):
    return sigma**2*np.ones(len(freqs))

ps = white_ps_pos
mean_ps = np.mean(ps,axis=0)
popt, pcov = curve_fit(white_ps_theo, freqs_pos, mean_ps, sigma=np.std(ps,axis=0), absolute_sigma=True, p0=200)

print('')
print('fitting avg white noise power spectrum with a constant:')
print('sigma=%5.5f \pm %5.5f\n' % tuple([popt[0],pcov[0,0]**0.5]))

plt.figure(figsize=(5,3.5))
plt.loglog(freqs_pos,mean_ps,'lightgray',label='50 realizations, avg')
plt.loglog(freqs_pos,white_ps_theo(freqs_pos,sigma),color='black',label='expectation')
plt.loglog(freqs_pos,white_ps_theo(freqs_pos,popt[0]),color='orangered',label='best fit')#,label=r'$\sigma=%5.1f$' % popt[0])
plt.legend(loc='lower left')
plt.xlabel(r'frequency [Hz]')
plt.ylabel(r'$P(f)$ [$\mu$K$^2$]')
plt.tight_layout()
plt.savefig('white_ps.pdf')
plt.show()

################################
### pure oof noise (PS)
################################

def oof_ps_theo(freqs, sigma, f_knee, slope):
    return sigma**2*(1+(f_knee/freqs)**(slope))

ps = oof_ps_sub
mean_ps = np.mean(ps,axis=0).flatten()

popt, pcov = curve_fit(oof_ps_theo, freqs_sub, mean_ps, sigma=np.std(ps,axis=0).flatten(), absolute_sigma=True, p0=tuple([200,0.01,1]))

print('')
print('fitting avg 1/f noise power spectrum:')
print('sigma=%5.3f \pm %5.3f, f_knee=%5.3f \pm %5.3f, slope=%5.3f \pm %5.3f\n' % tuple([popt[0],pcov[0,0]**0.5,popt[1],pcov[1,1]**0.5,popt[2],pcov[2,2]**0.5]))

plt.figure(figsize=(5,3.5))
plt.loglog(freqs_sub,mean_ps,'lightgray',label='50 realizations, avg')
plt.loglog(freqs_sub,oof_ps_theo(freqs_sub,sigma,f_knee,slope),color='black',label='expectation')
plt.loglog(freqs_sub,oof_ps_theo(freqs_sub,*popt),color='orangered',label='best fit')#,label=r'$\sigma=%5.1f$, $f_{knee}=%5.3f$, $\alpha=%5.1f$' % tuple(popt))
plt.loglog(freqs_sub,sigma**2*Spf_long,'orange',linestyle='dashed',label='SANEPIC weights, rescaled')
plt.ylim(bottom=2e4)
plt.legend(loc='upper right')
plt.xlabel(r'frequency [Hz]')
plt.ylabel(r'$P(f)$ [$\mu$K$^2$]')
plt.tight_layout()
plt.savefig('oof_ps.pdf')
plt.show()

################################
### pure oof noise (circularity)
################################

filepath = base_path + '/999mHz_20/00/tod_0A_00000.fits'
hdul = fits.open(filepath)

tod = hdul[1].data
oof_tod = tod['TOD']

oof_ps = np.abs(np.fft.fft(oof_tod))**2 / len(oof_tod)

pos_idxs = np.arange(500)
neg_idxs = -np.flip(pos_idxs) - 1

plt.figure(figsize=(5,3.5))
plt.plot(neg_idxs, tod[-500:],color='teal',label='last 500 samples')
plt.plot(pos_idxs, tod[:500], color='orange',label='first 500 samples')
plt.legend(loc='upper right')
plt.xlabel(r'sample')
plt.ylabel(r'TOD [$\mu$K]')
plt.tight_layout()
plt.savefig('circularity_ps.pdf')
plt.plot
