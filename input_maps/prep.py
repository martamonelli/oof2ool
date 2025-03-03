import numpy as np
import healpy as hp
import camb
from matplotlib import pyplot as plt

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif"
})

nside = 256
lmax = nside*2
ell = np.arange(lmax+1)

C2D = ell*(ell+1)/(2*np.pi)
D2C = np.append([1,1],1/C2D[2:])

pars = camb.read_ini('inifiles/planck_2018.ini')
results = camb.get_results(pars)
powers = results.get_cmb_power_spectra(pars, CMB_unit='muK')

# Returns the angular power spectra induced by scalar perturbations from 0 to lmax.
CAMB = powers['lensed_scalar']
Dls_in = np.array([CAMB[:,0], CAMB[:,1], CAMB[:,2], CAMB[:,3]])[:,:lmax+1]
#Dls_in = Cls_in*C2D
Dls_in[:,:2] *= 0
Cls_in = Dls_in*D2C

Dls_avg = np.empty((6,lmax+1))

realizations = 1 # = 100

for realization in np.arange(realizations):
    np.random.seed(realization)
    map_real = hp.synfast(Cls_in, nside, lmax=lmax, new=True)
    hp.write_map('maps/in_map_'+str(realization).zfill(2)+'.fits', map_real, overwrite=True)

    Cls_real = hp.anafast(map_real, lmax=lmax)
    Dls_real = Cls_real*C2D

    Dls_avg += Dls_real

    for comp in [3,0,1,2]:
        color = ['firebrick', 'darkorange', 'gold', 'lightsteelblue'][comp]
        plt.loglog(ell[2:], Dls_real[comp,2:], color=color, alpha=0.05, linewidth=0.5)
        plt.loglog(ell[2:], Dls_real[comp,2:], color='lightgray', alpha=0.05, linewidth=0.5)

Dls_avg = Dls_avg/realizations

for comp in [3,0,1,2]:
    color = ['firebrick', 'darkorange', 'gold', 'lightsteelblue'][comp]
    label = ['$TT$', '$EE$', '$BB$', '$TE$'][comp]
    plt.loglog(ell[2:], Dls_in[comp,2:], color='black', alpha=0.5)
    plt.loglog(ell[2:], Dls_avg[comp,2:], color=color, label=label)

plt.xlabel('$\ell$')
plt.ylabel('$\ell(\ell+1)C_\ell^{XY}/2\pi$ [$\mu$K$^2$]')
plt.legend(loc='lower right')
plt.savefig('in_Dls.pdf')
plt.clf()
