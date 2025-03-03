import numpy as np
import healpy as hp
import astropy.units as u
from numpy.random import MT19937, RandomState, SeedSequence

import litebird_sim as lbs

import sys
import os

#reading parameters initialized in the .sl script
spin_rate_rpm = float(sys.argv[1])
realization = int(sys.argv[2])
net_ukrts = float(sys.argv[3])
fknee_mhz = float(sys.argv[4])
alpha = float(sys.argv[5])
days = float(sys.argv[6])
ptg_only = bool(sys.argv[7])
tod_only = bool(sys.argv[8])

if ptg_only == tod_only:
    print('one (and only one) between ptg_only and tod_only should be true!')
    quit()

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

'''
#FIXME: needed to take care of weird LBsim convention, to be removed!
if tod_only:
    sampling_rate_hz /= 2
'''

samples_in_a_day = sampling_rate_hz * 3600 * 24
samples_in_a_year = samples_in_a_day * 365
samples_in_a_chunk = int(2 ** np.ceil(np.log2(days * samples_in_a_day))) #FIXME: used to force circularity, to be relaxed
nchunks = int(np.ceil(samples_in_a_year/samples_in_a_chunk))
duration_d = nchunks*samples_in_a_chunk/samples_in_a_day

print(f'samples in a chunk = {samples_in_a_chunk}')
print(f'{nchunks} chunks')
print(f'duration in days = {duration_d}')      

base_path = '/dss/dssfs02/lwp-dss-0001/pn36hu/pn36hu-dss-0000/oof2ool/'+str(int(days)).zfill(2)+'days'

sim = lbs.Simulation(
    base_path=base_path,
    start_time=0,
    duration_s=(duration_d * u.d).to("s").value,
    random_seed=realization, #from in line arguments
    imo = lbs.Imo(flatfile_location=lbs.PTEP_IMO_LOCATION),
    )

sim.set_scanning_strategy(
    scanning_strategy=lbs.SpinningScanningStrategy(
        spin_sun_angle_rad=np.deg2rad(spin_sun_angle_deg),           #defined above
        spin_rate_hz=1.0 / (1/spin_rate_rpm * u.min).to("s").value,  #defined above
        precession_rate_hz=1.0 / (192.348 * u.min).to("s").value,    #nominal LB
    )
)

instr = lbs.InstrumentInfo(
    name="LB",
    spin_boresight_angle_rad=np.deg2rad(95-spin_sun_angle_deg),      #beta = 95 - alpha (Yusuke)
)

sim.set_instrument(instr)

# We create four detectors oriented along the 0, π/2, π/4 and 3π/4 directions
dets=[
     lbs.DetectorInfo(name="0A", sampling_rate_hz=sampling_rate_hz, net_ukrts=net_ukrts, alpha=alpha, fknee_mhz=fknee_mhz, fwhm_arcmin=0.),
     lbs.DetectorInfo(name="0B", sampling_rate_hz=sampling_rate_hz, net_ukrts=net_ukrts, alpha=alpha, fknee_mhz=fknee_mhz, quat=lbs.quat_rotation_z(np.pi / 2), fwhm_arcmin=0.),
     lbs.DetectorInfo(name="1A", sampling_rate_hz=sampling_rate_hz, net_ukrts=net_ukrts, alpha=alpha, fknee_mhz=fknee_mhz, quat=lbs.quat_rotation_z(np.pi / 4), fwhm_arcmin=0.),
     lbs.DetectorInfo(name="1B", sampling_rate_hz=sampling_rate_hz, net_ukrts=net_ukrts, alpha=alpha, fknee_mhz=fknee_mhz, quat=lbs.quat_rotation_z(3 * np.pi / 4), fwhm_arcmin=0.),
     ]

print(f'running a {duration_d} days simulation for {len(dets)} detectors with spin rate = {spin_rate_rpm} rpm, net = {net_ukrts} ukrts, alpha = {alpha} and fknee = {fknee_mhz} mhz')

# We create the observations
obs = sim.create_observations(
    detectors=dets,
    tod_dtype=np.float64,  # Needed if you use the TOAST destriper
    #n_blocks_det=lbs.MPI_COMM_WORLD.size,
    n_blocks_time = nchunks,
    split_list_over_processes=False,
)

sim.prepare_pointings()
sim.precompute_pointings(pointings_dtype=np.float64) #Avinash's suggestion - to compute pointings in double prec

'''
in_map = hp.fitsfunc.read_map(base_path+'/../in/in_map_'+str(realization).zfill(2)+'.fits', field=[0,1,2])
nside = hp.npix2nside(len(in_map[0]))

lbs.scan_map_in_observations(
    obs,
    in_map,
    input_map_in_galactic = False, #to avoid rotations
    interpolation = "" #FIXME: fine as long as we only simulate noise
)

for cur_obs in sim.observations:
    # We include a CMB-only TOD in the signal
    cur_obs.cmb_only_tod = np.copy(cur_obs.tod)      # already in uK because in_map in uK
    cur_obs.tod *= 0.                                # setting tod to zero
'''

# Here we add 1/f noise using the detector noise parameters from the detector object
lbs.noise.add_noise_to_observations(obs, 'one_over_f', random=sim.random)

for cur_obs in sim.observations:
    # We include a oof-only TOD in the signal
    cur_obs.noise_tod = np.copy(cur_obs.tod)*1e6     # litebird_sim returns noise TOD in K, I'm saving them in uK

params = lbs.ExternalDestriperParameters(
    nside=64,                   #useless
    return_hit_map=False,       #useless
    return_binned_map=False,    #useless
    return_destriped_map=False, #useless
)

if ptg_only:
       lbs.save_simulation_for_madam(
           sim=sim,
           params=params,
           madam_subfolder_name=base_path+'/../pointings/'+str(int(spin_rate_rpm*100)).zfill(3)+'rpm_'+str(int(10*spin_sun_angle_deg)),
           components_to_bin="",
           save_pointings=True,
           save_tods=False,
       )

elif tod_only:
    tod_path = base_path+'/'+str(int(fknee_mhz)).zfill(3)+"mHz_"+str(int(alpha*10)).zfill(2)+'/'+str(realization).zfill(2) 

    # Iterate over all the maps we want to produce. For each of
    # them we specify the name of the subfolder where the Madam
    # files will be saved and the list of components to include
    for (subfolder, components_to_bin) in [(tod_path, [])]:
        lbs.save_simulation_for_madam(
            sim=sim,
            params=params,
            madam_subfolder_name=subfolder,
            components=["noise_tod"], #list all the components
            components_to_bin=components_to_bin,
            save_pointings=False,
            save_tods=True,
        )

    save_tods = False
