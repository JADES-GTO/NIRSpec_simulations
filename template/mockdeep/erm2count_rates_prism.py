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

# ===============================================================
# Module-wide variables
# ===============================================================
_name = "erm2count_rates.py"
_version = "1.0.3"

# =======================================================================
# Script specific variables
# =======================================================================
input_path = 'erm/'
print("# Input path: {:s}".format(input_path))
output_path = '/Users/ggiardin/JWST/IPSWork/MockDEEP/CountRates'
print("# Output path: {:s}".format(output_path))
target_erm =[['CLEAR-PRISM_MOS_dither_00_n0_001.erm', 'CLEAR-PRISM_MOS_dither_00_n1_001.erm', 'CLEAR-PRISM_MOS_dither_00_n2_000.erm'],
             ['CLEAR-PRISM_MOS_dither_1a_n0_000.erm', 'CLEAR-PRISM_MOS_dither_1a_n1_000.erm', 'CLEAR-PRISM_MOS_dither_1a_n2_000.erm'],
             ['CLEAR-PRISM_MOS_dither_2a_n0_000.erm', 'CLEAR-PRISM_MOS_dither_2a_n1_000.erm', 'CLEAR-PRISM_MOS_dither_2a_n2_000.erm']]

bkg_mos_erm = ['CLEAR-PRISM_MOS_background_000.erm',
               'CLEAR-PRISM_MOS_background_001.erm',
               'CLEAR-PRISM_MOS_background_002.erm']

bkg_slit_erm = ['CLEAR-PRISM_SLIT_background_000.erm',
                'CLEAR-PRISM_SLIT_background_001.erm',
                'CLEAR-PRISM_SLIT_background_002.erm']

for d in range(3):
        for n in range(3):
                print(target_erm[d][n])


#background_erm = []


# ===============================================================
# Spectrograph config
# ===============================================================
filetype = 'prism'
FWA = 'CLEAR'
GWA = 'PRSIM'

# ===============================================================
# Expousre parameters / straibgs
# ===============================================================
print("# Type of input file: {:s}".format(filetype))
print("# Spectral configuration: {:s}/{:s}".format(FWA, GWA))
base_nid = 1000
print("# Base NID number: {:d}".format(base_nid))
env = 'IPS'
pipeline_id = archive_toolbox._dic_env[env][0]
jlab_id = archive_toolbox._dic_env[env][1]
print("# Environment: {:s} ({:s} , {:s})".format(env, pipeline_id, jlab_id))
base_obs_id = 'DEEP-PRM'
print("# Base observation ID: {:s}".format(base_obs_id))
datetime_start = '20180518T120000.000'
datetime_end = '20180518T124800.000'
print("# Datetime start and end: {} : {}".format(datetime_start, datetime_end))
datetime_file = '2018-03-03T12h50m00'
print("# Datetime file: {}".format(datetime_file))
daily_folder = 'Day2018138'
print("# Daily folder: {:s}".format(daily_folder))
print("# =====================================")
print("# Detector configuration:")
nexp = 4
print("# Number of exposures per scene: {:d}".format(nexp))
ng = 19
nint = 2
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
dark = fpa106_toolbox.f_generate_dark(uniform=False, seed=seed)
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
    # Looping over nodding positions
    for inod in range(3):
        # -------------------------------------------------------------------
        # Name of the various files
        # -------------------------------------------------------------------
        g_mos = target_erm[idither][inod]
        b_mos = bkg_mos_erm[idither]
        b_slit = bkg_slit_erm[idither]

        # Reading in targets e-/rate maps
        g_erm = c_electronratemap.ElectronRateMap()
        g_erm.m_read_from_fits(os.path.join(input_path, g_mos))
        print('# Processing '+g_mos )

        # Reading in background maps
        # -------------------------------------------------------------------
        # Background ERM
        # -------------------------------------------------------------------
        b_erm = c_electronratemap.ElectronRateMap()
        b_erm.m_read_from_fits(os.path.join(input_path, b_mos))
        b_erm1 = c_electronratemap.ElectronRateMap()
        b_erm1.m_read_from_fits(os.path.join(input_path, b_slit))
        b_erm.m_add(b_erm1)

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
                                            force_trimming=False, old=True)

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

            filename = 'NRS_{:s}.oldcts.fits'.format(obs_id)
            gb_ctms[0].m_write_to_fits(os.path.join(output_path, daily_folder, folder_name, filename))
            gb_ctms[0].m_display_full(title=obs_id, filename=os.path.join(output_path, daily_folder, folder_name, '{:s}.pdf'.format(filename)),
                                      display=False, log_scale=True, category='data', comment=folder_name, figsize=(10.8, 8.5),
                                      dpi=100, axisbg='k', close=True, contour=False, cbar='horizontal', transparent=False,
                                      axis_font_size=None, time_stamp=True)
            for sca_index in range(2):
                filename = 'NRS{:s}_1_{:d}_SE_{:s}.cts.fits'.format(obs_id, 491 + sca_index, datetime_file)
                gb_ctms[sca_index+1].m_write_to_fits(os.path.join(output_path, daily_folder, folder_name, filename))


print('# Done')
