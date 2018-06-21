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
import datetime
from nrspydet.model import configurationfpa as c_configurationfpa
from nrspydet.model import multiaccum as c_multiaccum
from nrspydet.misc import fpa106_toolbox as fpa106_toolbox
from nrspysim.products import electronratemap as c_electronratemap
from nips.archives import toolbox as archive_toolbox
from nrspysim.misc import erm_to_sci as erm_to_sci

# ===============================================================
# Module-wide variables
# ===============================================================
_name = "p3_process_scene_20936pzs_01.py"
_version = "1.0.3"
_h = 6.62606957e-34
_c = 299792458.

# =======================================================================
# Parsing the input arguments
# =======================================================================
parser = argparse.ArgumentParser()
parser.add_argument('inpath', type=str, help='Name of the folder where the input files are present.')
listOfValidPositionsFWA = ['CLEAR', 'F070LP', 'F100LP', 'F170LP', 'F290LP', 'F140X', 'F110W']
parser.add_argument('FWA', type=str, help='Filter wheel position string.', choices=listOfValidPositionsFWA)
listOfValidPositionsGWA = ['PRISM', 'G140M', 'G235M', 'G395M', 'G140H', 'G235H', 'G395H', 'MIRROR']
parser.add_argument('GWA', type=str, help='Filter wheel position string.', choices=listOfValidPositionsGWA)
parser.add_argument('nexp', type=int, help='Number of exposures.')
parser.add_argument('ng', type=int, help='Number of groups (NRSIRS2 readout mode).')
parser.add_argument('outpath', type=str,
                    help='Name of the folder where the output files and folders will be generated.')
parser.add_argument('obsid', type=str, help='Base for the generation of the observation ID.')

parser.add_argument('-s', '--seed', type=int,
                    help='Optional seed (integer) for the random number generator.', default=None)
parser.add_argument('-b', '--background', type=str,
                    help='Suffix of the background electron-rate map.', default='background')
parser.add_argument('-u', '--unity', type=str, help='Suffix of the unity electron-rate map.', default='unity')
parser.add_argument('-o', '--object', type=str, help='Suffix of the object electron-rate map.', default='galaxy')
parser.add_argument('-f', '--factor', type=float,
                    help='Optional scaling factor to be applied to the object cube.', default=1.0)

parser.add_argument('-dts', '--datetime_start', type=str,
                    help='Start date and time.', default='20180303T120000.000')
parser.add_argument('-dte', '--datetime_end', type=str,
                    help='End date and time.', default='20180303T120000.000')
parser.add_argument('-dtf', '--datetime_file', type=str,
                    help='End date and time.', default='2018-03-03T12h15m00')
parser.add_argument('-day', '--daily_folder', type=str,
                    help='Daily folder.', default='Day2018062')

# =======================================================================
# Loading the input arguments
# =======================================================================
args = parser.parse_args()
argv = sys.argv
narg = len(argv) - 1
print("# =====================================")
print("# Running :{:s}".format(argv[0]))
print("# _name   : {:s}".format(_name))
print("# _version: {:s}".format(_version))
print("# Date    : {:s}".format((datetime.datetime.now()).isoformat()))
print("# =====================================")

# =======================================================================
# Parsing the arguments
# =======================================================================
input_path = args.inpath
print("# Input path: {:s}".format(input_path))
output_path = args.outpath
print("# Output path: {:s}".format(output_path))
FWA = args.FWA
GWA = args.GWA
if (GWA == 'PRISM'):
    filetype = 'prism'
elif (GWA == 'MIRROR'):
    filetype = 'image'
else:
    filetype = 'grating'
print("# Type of input file: {:s}".format(filetype))
print("# Spectral configuration: {:s}/{:s}".format(FWA, GWA))
base_nid = 1000
print("# Base NID number: {:d}".format(base_nid))
env = 'IPS'
pipeline_id = archive_toolbox._dic_env[env][0]
jlab_id = archive_toolbox._dic_env[env][1]
print("# Environment: {:s} ({:s} , {:s})".format(env, pipeline_id, jlab_id))
base_obs_id = args.obsid
print("# Base observation ID: {:s}".format(base_obs_id))
datetime_start = args.datetime_start
datetime_end = args.datetime_end
print("# Datetime start and end: {} : {}".format(datetime_start, datetime_end))
datetime_file = args.datetime_file
print("# Datetime file: {}".format(datetime_file))
daily_folder = args.daily_folder
print("# Daily folder: {:s}".format(daily_folder))
print("# =====================================")
print("# Detector configuration:")
nexp = args.nexp
print("# Number of exposures: {:d}".format(nexp))
print("")
ng = args.ng
configfpa = c_configurationfpa.ConfigurationFPA()
configfpa.m_set_mode('IRS2', 'full-frame')
configfpa.m_set_array_parameters()
configfpa.m_set_header('none', 'none', 'none')

multiaccum = c_multiaccum.Multiaccum()
multiaccum.m_set_configuration(configfpa)
multiaccum.m_set_header('none', 'none', 'none')
multiaccum.m_set(1, ng, 1)
multiaccum.m_info()
print("# =====================================")
print("# Optional arguments:")
background_suffix = args.background
print("# Suffix for the background electron-rate map: {:s}".format(background_suffix))
unity_suffix = args.unity
print("# Suffix for the unity electron-rate map: {:s}".format(unity_suffix))
object_suffix = args.object
print("# Suffix for the object electron-rate map: {:s}".format(unity_suffix))
seed = args.seed
if (seed is None):
    print("# No random generator seed provided.")
else:
    print("# Random generator seed: {:d}".format(seed))
    numpy.random.seed(seed=seed)
factor = args.factor
print("# Normalisation factor: {:8.4e}".format(factor))
print("# =====================================")

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
# Generating the names of the electron-rate maps
# =======================================================================
print("# Generating the names of the electron-rate maps.")
list_galaxy = []
for dither_index in range(4):
    name = 'ifs-{:s}-{:s}-{:s}-{:02d}.erm'.format(FWA.lower(), GWA.lower(),
                                            object_suffix, dither_index+1)
    list_galaxy.append(name)
print("# Galaxy ERMs: {}".format(list_galaxy))

list_background = []
for dither_index in range(4):
    list_names = []
    name = 'ifs-{:s}-{:s}-{:s}.erm'.format(FWA.lower(), GWA.lower(),
                                            background_suffix)
    list_names.append(name)
    name = 'mos-{:s}-{:s}-{:s}.erm'.format(FWA.lower(), GWA.lower(),
                                            background_suffix)
    list_names.append(name)
    name = 'slit-{:s}-{:s}-{:s}.erm'.format(FWA.lower(), GWA.lower(),
                                            background_suffix)
    list_names.append(name)
    list_background.append(list_names)
print("# Background ERMs: {}".format(list_background))

list_unity = []
for dither_index in range(4):
    list_names = []
    name = 'ifs-{:s}-{:s}-{:s}.erm'.format(FWA.lower(), GWA.lower(),
                                            unity_suffix)
    list_names.append(name)
    name = 'mos-{:s}-{:s}-{:s}.erm'.format(FWA.lower(), GWA.lower(),
                                            unity_suffix)
    list_names.append(name)
    name = 'slit-{:s}-{:s}-{:s}.erm'.format(FWA.lower(), GWA.lower(),
                                            unity_suffix)
    list_names.append(name)
    list_unity.append(list_names)
print("# Unity ERMs: {}".format(list_unity))

# =======================================================================
# Generating the count-rate maps
# =======================================================================
work = os.path.join(output_path, daily_folder)
try:
    os.makedirs(work)
except Exception as error:
    if (error.args[0] != 17):
        print('ERROR - Encountered an error when trying to create the daily folder:')
        print('ERROR - {:s}'.format(work))
        print('ERROR - {}'.format(error.args))
        raise ValueError

print("# Generating the count-rate maps.")
b_ctms = None
for dither_index in range(4):
    # -------------------------------------------------------------------
    # Name of the various files
    # -------------------------------------------------------------------
    g_ifs = list_galaxy[dither_index]
    b_ifs = list_background[dither_index][0]
    b_mos = list_background[dither_index][1]
    b_slit = list_background[dither_index][2]
    # -------------------------------------------------------------------
    # Background ERM
    # -------------------------------------------------------------------
    b_erm = c_electronratemap.ElectronRateMap()
    b_erm.m_read_from_fits(os.path.join(input_path, b_ifs))
    # b_erm.m_display(491, display=False, filename=os.path.join(output_path, 'bsifs.pdf'))
    work = c_electronratemap.ElectronRateMap()
    work.m_read_from_fits(os.path.join(input_path, b_mos))
    # work.m_display(491, display=False, filename=os.path.join(output_path, 'bsmos.pdf'))
    b_erm.m_add(work, overwrite=False)
    # b_erm.m_display(491, display=False, filename=os.path.join(output_path, 'bsifsmos.pdf'))
    work = c_electronratemap.ElectronRateMap()
    work.m_read_from_fits(os.path.join(input_path, b_slit))
    # work.m_display(491, display=False, filename=os.path.join(output_path, 'bslit.pdf'))
    b_erm.m_add(work)
    # b_erm.m_display(491, display=False, filename=os.path.join(output_path, 'bsifsmosslit.pdf'))
    # -------------------------------------------------------------------
    # Galaxy ERM
    # -------------------------------------------------------------------
    g_erm = c_electronratemap.ElectronRateMap()
    g_erm.m_read_from_fits(os.path.join(input_path, g_ifs))
    g_erm.m_add(b_erm)
    # -------------------------------------------------------------------
    # Generating the count-rate maps
    # -------------------------------------------------------------------
    if (dither_index == 0):
        b_ctms = erm_to_sci.f_standard(b_erm, multiaccum, noise_category='extended', noise_model=noise_model,
                                       scale=1.0, nexp=nexp, planes=None, bias=None, use_superbias=False,
                                       dark=dark, dark_sub=True, readout=readout, gain=gain, ctm=ctm,
                                       verbose=True, seed=-1, force_trimming=False, old=True)
    gb_ctms = erm_to_sci.f_standard(g_erm, multiaccum, noise_category='extended', noise_model=noise_model, scale=1.0,
                                   nexp=nexp, planes=None, bias=None, use_superbias=False, dark=dark, dark_sub=True,
                                   readout=readout, gain=gain, ctm=ctm, verbose=True, seed=-1,
                                   force_trimming=False, old=True)
    g_ctms = [gb_ctms[0] - b_ctms[0], gb_ctms[1] - b_ctms[1], gb_ctms[2] - b_ctms[2]]

    # -------------------------------------------------------------------
    # Generating the background exposure folder and the associated
    # exposures
    # -------------------------------------------------------------------
    if (dither_index == 0):
        nid = base_nid
        obs_id = '{:s}-B-{:s}-{:s}'.format(base_obs_id, FWA, GWA)
        folder_name = 'NRS{:s}_1_{:d}_{:s}_{:s}_{:s}_{:s}'.format(obs_id, nid, jlab_id, pipeline_id,
                                                                  datetime_start, datetime_end)
        work = os.path.join(output_path, daily_folder, folder_name)
        try:
            os.makedirs(work)
        except Exception as error:
            if (error.args[0] != 17):
                print('ERROR - Encountered an error when trying to create the exposure folder:')
                print('ERROR - {:s}'.format(work))
                print('ERROR - {}'.format(error.args))
                raise ValueError

        filename = 'NRS_{:s}.oldcts.fits'.format(obs_id)
        b_ctms[0].m_write_to_fits(os.path.join(output_path, daily_folder, folder_name, filename))
        b_ctms[0].m_display_full(title=obs_id, filename=os.path.join(output_path, daily_folder, folder_name, '{:s}.pdf'.format(filename)),
                       display=False, log_scale=True, category='data', comment=folder_name, figsize=(10.8, 8.5),
                       dpi=100, axisbg='k', close=True, contour=False, cbar='horizontal', transparent=False,
                       axis_font_size=None, time_stamp=True)

        for sca_index in range(2):
            filename = 'NRS{:s}_1_{:d}_SE_{:s}.cts.fits'.format(obs_id, 491+sca_index, datetime_file)
            b_ctms[sca_index+1].m_write_to_fits(os.path.join(output_path, daily_folder, folder_name, filename))
    # -------------------------------------------------------------------
    # Generating the galaxy + background exposure folders and the
    # associated exposures
    # -------------------------------------------------------------------
    nid = base_nid + dither_index + 10
    obs_id = '{:s}-GB-{:s}-{:s}-{:02d}'.format(base_obs_id, FWA, GWA, dither_index+1)
    folder_name = 'NRS{:s}_1_{:d}_{:s}_{:s}_{:s}_{:s}'.format(obs_id, nid, jlab_id, pipeline_id,
                                                              datetime_start, datetime_end)
    work = os.path.join(output_path, daily_folder, folder_name)
    try:
        os.makedirs(work)
    except Exception as error:
        if (error.args[0] != 17):
            print('ERROR - Encountered an error when trying to create the exposure folder:')
            print('ERROR - {:s}'.format(work))
            print('ERROR - {}'.format(error.args))
            raise ValueError
    filename = 'NRS_{:s}.oldcts.fits'.format(obs_id)
    gb_ctms[0].m_write_to_fits(os.path.join(output_path, daily_folder, folder_name, filename))
    gb_ctms[0].m_display_full(title=obs_id, filename=os.path.join(output_path, daily_folder, folder_name, '{:s}.pdf'.format(filename)),
                             display=False, log_scale=True, category='data', comment=folder_name, figsize=(10.8, 8.5),
                             dpi=100, axisbg='k', close=True, contour=False, cbar='horizontal', transparent=False,
                             axis_font_size=None, time_stamp=True)
    for sca_index in range(2):
        filename = 'NRS{:s}_1_{:d}_SE_{:s}.cts.fits'.format(obs_id, 491 + sca_index, datetime_file)
        gb_ctms[sca_index+1].m_write_to_fits(os.path.join(output_path, daily_folder, folder_name, filename))

    # -------------------------------------------------------------------
    # Generating the galaxy-only exposure folders and the
    # associated exposures
    # -------------------------------------------------------------------
    nid = base_nid + dither_index + 20
    obs_id = '{:s}-G-{:s}-{:s}-{:02d}'.format(base_obs_id, FWA, GWA, dither_index+1)
    folder_name = 'NRS{:s}_1_{:d}_{:s}_{:s}_{:s}_{:s}'.format(obs_id, nid, jlab_id, pipeline_id,
                                                              datetime_start, datetime_end)
    work = os.path.join(output_path, daily_folder, folder_name)
    try:
        os.makedirs(work)
    except Exception as error:
        if (error.args[0] != 17):
            print('ERROR - Encountered an error when trying to create the exposure folder:')
            print('ERROR - {:s}'.format(work))
            print('ERROR - {}'.format(error.args))
            raise ValueError
    filename = 'NRS_{:s}.oldcts.fits'.format(obs_id)
    g_ctms[0].m_write_to_fits(os.path.join(output_path, daily_folder, folder_name, filename))
    g_ctms[0].m_display_full(title=obs_id, filename=os.path.join(output_path, daily_folder, folder_name, '{:s}.pdf'.format(filename)),
                             display=False, log_scale=True, category='data', comment=folder_name, figsize=(10.8, 8.5),
                             dpi=100, axisbg='k', close=True, contour=False, cbar='horizontal', transparent=False,
                             axis_font_size=None, time_stamp=True)
    for sca_index in range(2):
        filename = 'NRS{:s}_1_{:d}_SE_{:s}.cts.fits'.format(obs_id, 491 + sca_index, datetime_file)
        g_ctms[sca_index+1].m_write_to_fits(os.path.join(output_path, daily_folder, folder_name, filename))

