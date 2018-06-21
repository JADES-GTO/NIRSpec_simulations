#########################################################################
# IDENT			p3_create_scene_background.py
# LANGUAGE		Python
# AUTHOR		P.FERRUIT          
# PURPOSE		Create the input scene for a spatially uniform background
#				based on a user-provided input spectrum.
#
# VERSION
# 1.0.0		13.02.2018 PF Creation
# 1.0.1     14.02.2018 PF Update
#   1) Added the creation of the MSA configuration files (perfect MSA
#   and OTIS failed-open map).
# 1.0.2     14.02.2018 PF Bug correction
#   1) Corrected the name of the object file in the source list and
#   deleted some residual code that was creating a second object file.
#	
#########################################################################

# =======================================================================
# Imports
# =======================================================================
import argparse
import datetime
import errno
import os.path
import sys
import numpy
from matplotlib import pyplot as plt
from nrspylib.io import spectrum as c_spectrum
from nrspylib.nirspec import nirspec as c_nirspec
from nrspylib.ote import ote as c_ote
from nrspysim.misc import projtools as projtools
from nrspysim.scenes import continuum as c_continuum
from nrspysim.scenes import object as c_object
from nrspysim.scenes import objectlist as c_objectlist
from nrspylib.msa import msaconfig as c_msaconfig
from nrspylib.msa import msaoperability as c_msaoperability
from nrspylib.msa import msaimage as c_msaimage

# =======================================================================
# Module-level variables
# =======================================================================
_name = 'p3_create_scene_background.py'
_version = '1.0.0'

_h = 6.62606957e-34
_c = 299792458.

if __name__ == '__main__':
    # =======================================================================
    # Declaring and loading the input arguments
    # =======================================================================
    parser = argparse.ArgumentParser()
    parser.add_argument('inspectrum', type=str, help='Input spectrum in units of W m-2 m-1 arcsec-2.')
    parser.add_argument('outpath', type=str, help='Output path where the scene folder will be created.')
    parser.add_argument('outfolder', type=str, help='Name of the output folder that will contain the scene.')

    parser.add_argument('-n', '--normalisation', type=float,
                        help='Normalisation factor to be applied to the input spectrum.', default=1.0)
    parser.add_argument('-bb', '--bounding_box', type=float, nargs=4,
                        help='Bounding-box for the extension of the background source in V2V3 coordinates (units = arcminutes, [xmin xmax ymin ymax]).',
                        default=[3.5, 9.5, -10.0, -4.2])
    parser.add_argument('-p', '--pathmodel', type=str, help='Path to the models.',
                        default='/Users/pferruit/PycharmProjects/nirspec27/JWST_Python/data/IQLAC')
    parser.add_argument('-m', '--model', type=str, help='Name of the IQLAC instrument model.',
                        default='NIRS_FM2_05_CV3_FIT1')

    # =======================================================================
    # Loading the input arguments
    # =======================================================================
    args = parser.parse_args()
    argv = sys.argv
    narg = len(argv) - 1
    print("=====================================")
    print(argv[0])
    print("Version: {:s}".format(_version))
    print((datetime.datetime.now()).isoformat())
    print("=====================================")

    # =======================================================================
    # Parsing the input arguments - Paths and names
    # =======================================================================
    input_filename = args.inspectrum
    print("# Input spectrum: {:s}".format(input_filename))
    factor = args.normalisation
    print("# Normalisation factor to be applied to the input spectrum: {:5.3f}".format(factor))
    output_path = args.outpath
    print("# Output path: {:s}".format(output_path))
    output_folder = args.outfolder
    print("# Output scene folder: {:s}".format(output_folder))
    background_bounding_box = 60. * numpy.array(args.bounding_box, dtype=numpy.float)
    print("# A background source will be generated.")
    print("#  * Bounding box: [{:5.3f} , {:5.3f} , {:5.3f} , {:5.3f}] ([xmin,xmax,ymin,ymax], in arcmin.".format(
        background_bounding_box[0] / 60,
        background_bounding_box[1] / 60,
        background_bounding_box[2] / 60,
        background_bounding_box[3] / 60))
    FWA = 'CLEAR'
    GWA = 'MIRROR'
    print("# NIRSpec configuration: {:s}/{:s}".format(FWA, GWA))
    model_path = args.pathmodel
    print("# Path to the model folders: {:s}".format(model_path))
    model_name = args.model
    print("# Name of the  geometrical model of NIRSpec and the OTE: {:s}.".format(model_name))

    # =======================================================================
    # Geometrical model initialisation
    # =======================================================================
    print("# Initialising the geometrical models (OTE + NIRSpec MOS + NIRSpec IFU).")
    NIRSpecMOS = c_nirspec.NIRSpec()
    NIRSpecMOS.m_initializeMOS(FWA, GWA, path=model_path, model=model_name)
    NIRSpecIFU = c_nirspec.NIRSpec()
    slice_id = 0
    NIRSpecIFU.m_initializeIFU(FWA, GWA, slice_id, path=model_path, model=model_name)
    OTE = c_ote.OTE()
    OTE.m_initializeFromModel(model_path, model_name)

    # =======================================================================
    # Loading the input spectrum
    # =======================================================================
    print("# Loading the input spectrum.")
    spectrumObject = c_spectrum.Spectrum()
    spectrumObject.m_readFromSimpleFITS(input_filename)
    nblbda = spectrumObject.nbPixels
    wavelength = spectrumObject.m_getLbda()
    spectrum = spectrumObject.m_getValues()
    # No absorption lines
    flag = numpy.zeros(nblbda, dtype=numpy.bool)
    ew = numpy.zeros(nblbda, dtype=float)

    # =======================================================================
    # Preparing the output folder
    # =======================================================================
    print("# Preparing the output folder.")
    full_output_path = os.path.join(output_path, output_folder)
    try:
        os.makedirs(full_output_path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
        else:
            print("# CAUTION: the output folder already exists. Contents will be overwritten.")

    # =======================================================================
    # V2V3 coordinate grid - in arcsec, V2V3 plane
    # =======================================================================
    print("# Generating the grid of position used to define the background object.")
    n_x = 5
    n_y = 5
    min_x = background_bounding_box[0]
    max_x = background_bounding_box[1]
    min_y = background_bounding_box[2]
    max_y = background_bounding_box[3]
    grid_x = numpy.zeros(n_x, dtype=numpy.float)
    grid_y = numpy.zeros(n_y, dtype=numpy.float)
    for indexx in range(n_x):
        grid_x[indexx] = min_x + (max_x - min_x) * indexx / (n_x - 1)
    for indexy in range(n_y):
        grid_y[indexy] = min_y + (max_y - min_y) * indexy / (n_y - 1)

    # =======================================================================
    # Associated figure
    # =======================================================================
    print("# Generating the associated figure.")
    mesh_y, mesh_x = numpy.meshgrid(grid_y, grid_x)
    figure_name = os.path.join(output_path, output_folder, 'background_source.fits.pdf')
    figure = plt.figure()
    axes = figure.gca()
    axes.set_aspect('equal')
    ima = axes.scatter(mesh_x, mesh_y, s=30, alpha=0.5)
    projtools.f_display_quadrants(OTE, NIRSpecMOS, axes, quadrants='ALL', wavelength=2.5e-6,
                                          plane='V2V3', patches=True, arcsec=True)

    projtools.f_display_slits(OTE, NIRSpecMOS, axes, slits='ALL', wavelength=2.5e-6, plane='V2V3',
                              patches=False, arcsec=True)
    # plt.ylim([1e-21,1e-15])
    axes.set_xlabel('x-axis position (in arcsec, V2V3 plane)')
    axes.set_ylabel('x-axis position (in arcsec, V2V3 plane)')
    plt.title("Background source", fontsize=11)
    plt.grid(True)
    plt.figtext(0.99, 0.01, (datetime.datetime.now()).isoformat(), ha='right', style='italic', size=6)
    plt.savefig(figure_name)
    plt.close(figure)

    # =======================================================================
    # Creating the object file (BS+CONT)
    # =======================================================================
    print("# Creating the object.")
    input_plane = 'SKY'
    object_type = 'BS'
    continuum = c_continuum.Continuum()
    continuum.m_set_input_plane(input_plane)
    continuum.m_set_wave(wavelength)
    for index_x in range(n_x):
        for index_y in range(n_y):
            continuum.m_add_point(grid_x[index_x], grid_y[index_y], spectrum, flag, ew)
    current_object = c_object.Object()
    reference = input_filename
    author = '{:s} version {:s}'.format(argv[0], _version)
    description = 'user-provided data cube'
    current_object.m_set_header(reference, author, description)
    current_object.m_set_parameters(input_plane, object_type)
    current_object.m_add_continuum(continuum, overwrite=False)
    filename = os.path.join(output_path, output_folder, 'background.fits')
    current_object.m_write_to_fits(filename)

    # =======================================================================
    # Preparing the list of objects
    # =======================================================================
    print("# Creating the list of objects.")
    reference = input_filename
    author = '{:s} version {:s}'.format(argv[0], _version)
    description = 'BG+CONT uniform background scene'
    object_list = c_objectlist.ObjectList()
    object_list.m_set_header(reference, author, description)
    object_list.m_set_input_plane(input_plane)
    object_list.m_add_object(0.0, 0.0, 0.0, factor, 'BS', 'uniform background object', 'background.fits')
    filename = os.path.join(output_path, output_folder, '{:s}.list'.format(output_folder))
    object_list.m_write_to_fits(filename)

    # =======================================================================
    # Generating the MSA configuration file
    # =======================================================================
    print("# Generating an MSA configuration file #1 (all closed, perfect MSA).")
    description = "MSA configuration file for scene {:s} - perfect MSA".format(output_folder)
    author = '{:s} version {:s}'.format(argv[0], _version)
    reference = 'all closed - perfect micro-shutter array.'
    nq = 4
    nx = 365
    ny = 171
    configuration = c_msaconfig.MSAConfig()
    configuration.m_setParameters(description, author, reference, nq, nx, ny)
    configuration.m_setAll(0)
    configuration.m_writeToFITS(os.path.join(output_path, output_folder, 'all_closed_perfect_msa.mos'))

    # =======================================================================
    # Generating the MSA configuration report
    # =======================================================================
    print("# Generating the associated report.")
    # Dummy failure list wihtout any failure.
    description = "perfect MSA operability map"
    author = '{:s} version {:s}'.format(argv[0], _version)
    reference = "None"
    nq = 4
    nx = 365
    ny = 171
    list = c_msaoperability.MSAOperability()
    list.m_setParameters(description, author, reference, nq, nx, ny)
    list.m_setAllToFunctional(2)
    # Image including the dummy failure list
    image = c_msaimage.MSAImage()
    image.m_construct(configuration, list)
    image.m_generateReportPDF(output_folder, os.path.join(output_path, output_folder,
                                                          'all_closed_perfect_msa.mos.pdf'))

    # =======================================================================
    # Loading the OTIS MSA operability map
    # =======================================================================
    print("# Loading the OTIS MSA operability map.")
    filename = './data/nrs_msop_CHK_001_20170917.msl'
    description = "dummy description"
    author = '{:s} version {:s}'.format(argv[0], _version)
    reference = "dummy reference"
    nq = 4
    nx = 365
    ny = 171
    list = c_msaoperability.MSAOperability()
    list.m_setParameters(description, author, reference, nq, nx, ny)
    list.m_readFromFITS(filename)

    # =======================================================================
    # Generating the MSA configuration file
    # =======================================================================
    print("# Generating an MSA configuration file #2 (all closed, OTIS MSA).")
    description = "MSA configuration file for scene {:s} - OTIS operability map".format(output_folder)
    author = '{:s} version {:s}'.format(argv[0], _version)
    reference = filename
    nq = 4
    nx = 365
    ny = 171
    configuration = c_msaconfig.MSAConfig()
    configuration.m_setParameters(description, author, reference, nq, nx, ny)
    # All closed
    configuration.m_setAll(0)
    for current_q in range(4):
        list_i, list_j = list.m_getListOfFailedOpenShutters(current_q + 1)
        list_n = configuration.m_getNumber(list_i, list_j)
        for current_n in list_n:
            configuration.status[current_q, current_n] = 1
    configuration.m_writeToFITS(os.path.join(output_path, output_folder, 'all_closed_OTIS.mos'))

    # =======================================================================
    # Generating the MSA configuration report
    # =======================================================================
    print("# Generating the associated report.")
    # Dummy failure list wihtout any failure. Doing that to make sure
    # the failed-open shutters are resent in the configuration
    list = c_msaoperability.MSAOperability()
    list.m_setParameters(description, author, reference, nq, nx, ny)
    list.m_setAllToFunctional(2)
    # Image including the dummy failure list
    image = c_msaimage.MSAImage()
    image.m_construct(configuration, list)
    image.m_generateReportPDF(output_folder, os.path.join(output_path, output_folder,
                                                          'all_closed_OTIS.mos.pdf'))
