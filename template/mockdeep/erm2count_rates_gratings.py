#!/usr/bin/env python
#########################################################################
# IDENT			p3_process_scene_20936pzs_01.py
# LANGUAGE		Python
# AUTHOR		P.FERRUIT
# PURPOSE		Process the exposures generated using the IPS and the
#               family of input scenes G20936pzs_01x simulating a
#               4-point-nod dither sequence with the galaxy 20936pzs
#               cube prepared by Camilla.
#               #01 = standard synthetic FPA106 model from
#               nrspydet.misc.fpa106_toolbox.
#
# VERSION
# 1.0.0 03.03.2018 PF Creation
# 1.0.1 07.03.2018 PF Bug correction
#   1) The IPC was not included in the computation. Corrected the calls
#   to the erm_to_sci.f_standard() function
# 1.0.2 20.03.2018 PF Bug correction.
#   1) I had forgotten the daily folders in the structure of the
#   output folders. Corrected that.
# 1.0.3 27.04.2018 PF Bug correction & update
#   1) The order of the JLAB and pipeline IDs in the folder names was
#   swapped compared to the expected one. Corrected that.
#
#########################################################################
# ===============================================================
# Imports
# ===============================================================
import argparse
import datetime
import os.path
import sys
import numpy
from nrspydet.model import configurationfpa as c_configurationfpa
from nrspydet.model import multiaccum as c_multiaccum
from nrspydet.misc import fpa106_toolbox as fpa106_toolbox
from nrspysim.products import electronratemap as c_electronratemap
from nips.archives import toolbox as archive_toolbox
from nrspysim.misc import erm_to_sci as erm_to_sci
from astropy.io import fits
import numpy as np

# ===============================================================
# Module-wide variables
# ===============================================================
_name = "erm2count_rates_gratings.py"
_version = "1.0.3"

# =======================================================================
# Script specific variables
# =======================================================================
base_path = '/Users/ggiardin/JWST/IPSWork/00MockDEEP/'
input_path = base_path+'erm/'
print("# Input path: {:s}".format(input_path))
output_path = base_path+'crates/'
#output_path = base_path+'coadd_crates/'

print("# Output path: {:s}".format(output_path))

# ===============================================================
# Spectrograph config
# ===============================================================
#Change as applicable
FWA = 'F290LP'
GWA = 'G395M'
daily_folder = 'Day2018139'

target_erm =[[FWA+'_'+GWA+'_MOS_dither_00_n0_000.erm', FWA+'_'+GWA+'_MOS_dither_00_n1_000.erm', FWA+'_'+GWA+'_MOS_dither_00_n2_000.erm'],
             [FWA+'_'+GWA+'_MOS_dither_1a_n0_000.erm', FWA+'_'+GWA+'_MOS_dither_1a_n1_000.erm',	FWA+'_'+GWA+'_MOS_dither_1a_n2_000.erm'],
             [FWA+'_'+GWA+'_MOS_dither_2a_n0_000.erm', FWA+'_'+GWA+'_MOS_dither_2a_n1_000.erm',	FWA+'_'+GWA+'_MOS_dither_2a_n2_000.erm']]

bkg_mos_erm = [FWA+'_'+GWA+'_MOS_background_000.erm',
               FWA+'_'+GWA+'_MOS_background_001.erm',
               FWA+'_'+GWA+'_MOS_background_002.erm']

bkg_slit_erm = [FWA+'_'+GWA+'_SLIT_background_000.erm',
                FWA+'_'+GWA+'_SLIT_background_001.erm',
                FWA+'_'+GWA+'_SLIT_background_002.erm']


# Filenames of real darks used to set the detectors quality flags
dark491fname = base_path+'data/NRSDET-DARK-IRS2-7245080538_3_491_SE_2017-09-02T10h49m27.cts.fits'
dark492fname = base_path+'data/NRSDET-DARK-IRS2-7245080538_3_492_SE_2017-09-02T10h49m27.cts.fits'

for d in range(3):
        for n in range(3):
                print(target_erm[d][n])


#background_erm = []




# ===============================================================
# Expousre parameters / straibgs
# ===============================================================
print("# Spectral configuration: {:s}/{:s}".format(FWA, GWA))
base_nid = 1000
print("# Base NID number: {:d}".format(base_nid))
env = 'IPS'
pipeline_id = archive_toolbox._dic_env[env][0]
jlab_id = archive_toolbox._dic_env[env][1]
print("# Environment: {:s} ({:s} , {:s})".format(env, pipeline_id, jlab_id))
base_obs_id = 'DEEP-'+GWA
print("# Base observation ID: {:s}".format(base_obs_id))
datetime_start = '20180518T120000.000'
datetime_end = '20180518T124800.000'
print("# Datetime start and end: {} : {}".format(datetime_start, datetime_end))
datetime_file = '2018-03-03T12h50m00'
print("# Datetime file: {}".format(datetime_file))
print("# Daily folder: {:s}".format(daily_folder))
print("# =====================================")
print("# Detector configuration:")
nexp = 1
print("# Number of exposures per scene: {:d}".format(nexp))
ng = 19
nint = 2
# nint=6
nf = 5 #IRS2
print(" Number of groups {:d} - number of intgrations {:d}".format(ng, nint))
print("")

configfpa = c_configurationfpa.ConfigurationFPA()
configfpa.m_set_mode('IRS2', 'full-frame')
configfpa.m_set_array_parameters()
configfpa.m_set_header('none', 'none', 'none')

multiaccum = c_multiaccum.Multiaccum()
multiaccum.m_set_configuration(configfpa)
multiaccum.m_set_header('none', 'none', 'none')
multiaccum.m_set(nint, ng, nf)
multiaccum.m_info()

seed = None
if (seed is None):
    print("# No random generator seed provided.")
else:
    print("# Random generator seed: {:d}".format(seed))
    numpy.random.seed(seed=seed)

# to Pierre, What is this?
# factor = args.factor
# print("# Normalisation factor: {:8.4e}".format(factor))
# print("# =====================================")

# =======================================================================
# Generating the instances of the reference files (nrspydet.references)
# =======================================================================
print("# Generating the instances of the reference files.")
dark = fpa106_toolbox.f_generate_dark(uniform=False, mode='IRS2', seed=seed)
gain = fpa106_toolbox.f_generate_gain('FULL-FRAME')
ctm = fpa106_toolbox.f_generate_ctm()
readout = fpa106_toolbox.f_generate_readout_noise('IRS2')
noise_model = fpa106_toolbox.f_generate_noise_model('IRS2')
print("# =====================================")


# =======================================================================
# Generating the count-rate maps
# =======================================================================
work_folder = os.path.join(output_path, daily_folder)
try:
    os.makedirs(work_folder)
except Exception as error:
    if (error.args[0] != 17):
        print('ERROR - Encountered an error when trying to create the daily folder:')
        print('ERROR - {:s}'.format(work_folder))
        print('ERROR - {}'.format(error.args))
        raise ValueError

print("# Generating the count-rate maps.")
b_ctms = None

# Looping over dither pointings
for idither in range(3):

    # Reading in background maps
    b_mos = bkg_mos_erm[idither]
    b_slit = bkg_slit_erm[idither]

    # -------------------------------------------------------------------
    # Background ERM
    # -------------------------------------------------------------------
    b_erm = c_electronratemap.ElectronRateMap()
    b_erm.m_read_from_fits(os.path.join(input_path, b_mos))
    b_erm1 = c_electronratemap.ElectronRateMap()
    b_erm1.m_read_from_fits(os.path.join(input_path, b_slit))
    b_erm.m_add(b_erm1)

    # Looping over nodding positions
    for inod in range(3):
        # -------------------------------------------------------------------
        # Name of the various files
        # -------------------------------------------------------------------
        g_mos = target_erm[idither][inod]

        # Reading in targets e-/rate maps
        g_erm = c_electronratemap.ElectronRateMap()
        g_erm.m_read_from_fits(os.path.join(input_path, g_mos))
        print('# Processing '+g_mos )

        # Adding target and background together
        g_erm.m_add(b_erm)

        # -------------------------------------------------------------------
        # Generating the exposure folders and the associated exposures - 4 exposures for pointing
        # -------------------------------------------------------------------
        for iexp in range(nexp):
            # Note nexp in this function is only used to scale noise -> nexp=1
            gb_ctms = erm_to_sci.f_standard(g_erm, multiaccum, noise_category='extended', noise_model=noise_model,
                                            scale=1.0,
                                            nexp=1, planes=None, bias=None, use_superbias=False, dark=dark, dark_sub=True,
                                            readout=readout, gain=gain, ctm=ctm, verbose=False, seed=seed,
                                            force_trimming=True, old=True)

            count = inod*nexp + iexp
            nid = base_nid + 100*idither + count +1
            obs_id = '{:s}-{:02d}'.format(base_obs_id, count+1)
            folder_name = 'NRS{:s}_1_{:d}_{:s}_{:s}_{:s}_{:s}'.format(obs_id, nid, jlab_id, pipeline_id,
                                                                      datetime_start, datetime_end)
            exp_folder = os.path.join(output_path, daily_folder, folder_name)

            try:
                os.makedirs(exp_folder)
            except Exception as error:
                if (error.args[0] != 17):
                    print('ERROR - Encountered an error when trying to create the exposure folder:')
                    print('ERROR - {:s}'.format(exp_folder))
                    print('ERROR - {}'.format(error.args))
                    raise ValueError

            print('# Created output folder ' + exp_folder)

            flag_arrays = []

            # Reading in quality flags from one of our darks
            hdu = fits.open(dark491fname)
            flag_arrays.append(hdu[3].data)
            print('# Reading in FLAG array 1')

            hdu = fits.open(dark492fname)
            flag_arrays.append(hdu[3].data)
            print('# Reading in FLAG array 2')
            hdu.close()

            gb_ctms[0].m_set_quality_arrays(np.swapaxes(flag_arrays[0][4:2044, 4:2044], 0, 1),
                                            np.swapaxes(flag_arrays[1][4:2044, 4:2044], 0, 1), overwrite=True)

            filename = 'NRS_{:s}.oldcts.fits'.format(obs_id)

            gb_ctms[0].m_write_to_fits(os.path.join(output_path, daily_folder, folder_name, filename))
            gb_ctms[0].m_display_full(title=obs_id, filename=os.path.join(output_path, daily_folder, folder_name,
                                                                          '{:s}.pdf'.format(filename)),
                                      display=False, log_scale=True, category='data', comment=folder_name,
                                      figsize=(10.8, 8.5),
                                      dpi=100, axisbg='k', close=True, contour=False, cbar='horizontal',
                                      transparent=False,
                                      axis_font_size=None, time_stamp=True)

            filename = 'NRS_{:s}.oldcts.fits'.format(obs_id)
            gb_ctms[0].m_write_to_fits(os.path.join(output_path, daily_folder, folder_name, filename))
            gb_ctms[0].m_display_full(title=obs_id, filename=os.path.join(output_path, daily_folder, folder_name, '{:s}.pdf'.format(filename)),
                                      display=False, log_scale=True, category='data', comment=folder_name, figsize=(10.8, 8.5),
                                      dpi=100, axisbg='k', close=True, contour=False, cbar='horizontal', transparent=False,
                                      axis_font_size=None, time_stamp=True)
            mode = 'unknown'
            slit = 'unknown'
            aperture = 'unknown'
            sca_ids = ['NRS1', 'NRS2']
            read_out = 'NRSIRS2'
            for sca_index in range(2):
                filename = 'NRS{:s}_1_{:d}_SE_{:s}.cts.fits'.format(obs_id, 491 + sca_index, datetime_file)
                gb_ctms[sca_index + 1].m_set_keywords(nid, 'IPS', sca_ids[sca_index], mode, slit, aperture,
                                                      FWA, GWA, read_out, False)
                gb_ctms[sca_index + 1].quality = flag_arrays[sca_index]
                gb_ctms[sca_index+1].m_write_to_fits(os.path.join(output_path, daily_folder, folder_name, filename))


print('# Done')
