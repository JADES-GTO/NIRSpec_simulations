#########################################################################
# IDENT			p3_create_scene_20936pzs.py
# LANGUAGE		Python
# AUTHOR		P.FERRUIT          
# PURPOSE		Create the family of input scenes G20936pzs_01x
#               simulating a 4-point-nod dither sequence with
#               the galaxy 20936pzs cube prepared by Camilla.
#
# VERSION
# 1.0.0 13.02.2018 PF Creation
# 1.0.1 14.02.2018 PF Bug correction
#   1) The objects were declared as point-sources in the list instead
#   of extended sources.
#########################################################################
# ===============================================================
# Imports
# ===============================================================
import argparse
import datetime
import errno
import os.path
import sys
import numpy
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
from nrspylib.io import spectrum as c_spectrum
from nrspylib.nirspec import nirspec as c_nirspec
from nrspylib.ote import ote as c_ote
from nrspysim.misc import projtools as projtools
from nrspysim.scenes import continuum as c_continuum
from nrspysim.scenes import object as c_object
from nrspysim.scenes import objectlist as c_objectlist
from nips.references import siaf_reference as ref_siaf

# ===============================================================
# Module-wide variables
# ===============================================================
_name = "p3_create_scene_20936pzs.py"
_version = "1.0.1"
_h = 6.62606957e-34
_c = 299792458.

# ===============================================================
# Function to read the input data cube
# ===============================================================
def f_read(filename):
    # .......................................................
    # Accessing the data file containing the wavelength
    # scale
    # .......................................................
    inputfile = open('{:s}.dat'.format(filename), 'r')
    inputlines = inputfile.readlines()
    wavelength= 1e-10 * numpy.array(inputlines[1:], dtype=numpy.float)
    inputfile.close()
    nw = numpy.shape(wavelength)[0]
    inputfile = 0
    inputlines = 0
    # .......................................................
    # Opening the input FITS file
    # .......................................................
    fitsobj = pyfits.open(filename)
    # .......................................................
    # Accessing the data extension
    # .......................................................
    cube = numpy.swapaxes((fitsobj[0]).data, 0, 2)
    print(numpy.shape(cube))
    # .......................................................
    # Cleaning up
    # .......................................................
    fitsobj.close()
    fitsobj = 0
    # .......................................................
    # Returning the wavelength and cube arrays
    # .......................................................
    return wavelength,cube



if __name__ == '__main__':
    # =======================================================================
    # Preapring the list of input arguments
    # =======================================================================
    parser = argparse.ArgumentParser()
    parser.add_argument('inpath', type=str, help='Path to the folder containing the input file.')
    parser.add_argument('infile', type=str, help='Name of the input file.')
    parser.add_argument('intable', type=str, help='Name of the input SIAF table.')
    parser.add_argument('outpath', type=str, help='Path where the output files will be generated.')
    parser.add_argument('outprefix', type=str, help='Prefix used to generate the name of the output folder.')
    list_positions_fwa = ['CLEAR', 'F070LP', 'F100LP', 'F170LP', 'F290LP', 'F140X', 'F110W']
    parser.add_argument('fwa', type=str, help='Filter wheel position name.', choices=list_positions_fwa)
    parser.add_argument('-xs', '--xscale', type=float,
                         help='image scale conversion along the x-axis (spectral direction; arcsec/slice)',
                         default=0.1033)
    parser.add_argument('-ys', '--yscale', type=float,
                         help='image scale conversion along the y-axis (spatial direction; arcsec/pixel)',
                         default=0.1073)
    parser.add_argument('-f', '--factor', type=float,
                         help='multiplicative normalisation factor for the cube (> 0.0)',
                         default=0.1073)
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
    # Parsing the input arguments
    # =======================================================================
    input_path = args.inpath
    print("# Input path : {:s}".format(input_path))
    input_file = args.infile
    print("# Input file : {:s}".format(input_file))
    factor = args.factor
    if (factor <= 0.0):
        print("ERROR - The normalisation factor must be strictly positive.")
        print("ERROR - Got: {:12.6e}".format(factor))
    print("# Normalisation factor fro the cube: {:12.6e} (unitless, multiplicative)")
    step_x = 0.1
    step_y = 0.1
    spaxel = step_x * step_y
    print("# Spaxel size: {:3.2f} x {:3.2f} (arcsec x arcsec)".format(step_x, step_y))
    input_table = args.intable
    print("# Input SIAF table: {:s}".format(input_table))
    output_path = args.outpath
    print("# Output path: {:s}".format(output_path))
    output_prefix = args.outprefix
    print("# Prefix for the output scene folder name: {:s}".format(output_prefix))
    x_scale = args.xscale
    y_scale = args.yscale
    print("# Image scale conversion along the spectral direction (x-axis) : {:6.2f} mas/slicer".format(1e3*x_scale))
    print("# Image scale conversion along the spatial direction (x-axis)  : {:6.2f} mas/pixel".format(1e3*y_scale))
    FWA = args.fwa
    GWA = 'MIRROR'
    print("# NIRSpec configuration: {:s}/{:s}".format(FWA, GWA))
    model_path = args.pathmodel
    print("# Path to the model folders: {:s}".format(model_path))
    model_name = args.model
    print("# Name of the  geometrical model of NIRSpec and the OTE: {:s}.".format(model_name))

    # =======================================================================
    # Loading the input data
    # =======================================================================
    print("# Loading the input file: {:s}".format(input_file))
    filename = './data/20936/2D_gal20936_prior_z.fits'
    # f_readCubeFromFORTRAN(filename, 10201)
    wavelength,cube = f_read(os.path.join(input_path, input_file))
    n_x, n_y, n_w = numpy.shape(cube)
    print("#    * Cube spatial size: ( {:d} , {:d} )".format(n_x, n_y))
    print("#    * Cube spectral size: {:d}".format(n_w))
    start_x = - 0.5 * (n_x - 1) * step_x
    start_y = - 0.5 * (n_y - 1) * step_y
    print("#    * Cube step size: ( {:5.3f} , {:5.3f} )".format(step_x, step_y))
    print("#    * Cube start value: ( {:5.3f} , {:5.3f} )".format(start_x, start_y))
    # Converting from erg s-1 cm-2 A-1 spaxel-1 to W m-2 m-1 arcsec-2
    # and appying a user-provided factor if relevant.
    cube = 1e7 * factor * cube / spaxel**2
    # Creating additional associated arrays
    flag = numpy.zeros(n_w, dtype=numpy.bool)
    ew = numpy.zeros(n_w, dtype=float)

    # =======================================================================
    # Loading the input SIAF table
    # =======================================================================
    print("# Loading the SIAF table.")
    siaf = ref_siaf.SiafRef()
    siaf.m_read_from_fits(os.path.join(input_path, input_table))

    # =======================================================================
    # Preparing the dither sequence
    # =======================================================================
    print("# Preparing the dither sequence.")
    nb_points = 4
    # Dither sequence in units of slices and pixels - 4-point-dither
    # https://jwst-docs.stsci.edu/display/JTI/NIRSpec+IFU+Dither+and+Nod+Patterns
    gridx_relative = numpy.array([-2.25, 1.25, -0.25, -0.75], dtype=numpy.float)
    gridy_relative = numpy.array([-0.75, 1.75, -1.25, 2.25], dtype=numpy.float)
    # Dither sequence in arcsec and in the aperture ideal coordinate system
    gridx_ideal = gridx_relative * x_scale
    gridy_ideal = gridy_relative * y_scale
    # Dither sequence in the v2V3 coordinate system
    gridx_v2v3,gridy_v2v3 = siaf.m_ideal_to_v2v3('NRS_FULL_IFU', [gridx_ideal, gridy_ideal], absolute=True)
    for index_point in range(nb_points):
        print("# [{:02d}] ({:5.3f},{:5.3f}) ({:5.3f},{:5.3f}) ({:8.3f},{:8.3f})".format(index_point+1,
                                                                                      gridx_relative[index_point],
                                                                                      gridy_relative[index_point],
                                                                                      gridx_ideal[index_point],
                                                                                      gridy_ideal[index_point],
                                                                                      gridx_v2v3[index_point],
                                                                                      gridy_v2v3[index_point]))

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

    xcenter_v2v3,ycenter_v2v3 = projtools.f_project_ifu(OTE, NIRSpecIFU)
    # =======================================================================
    # Main loop on the scenes
    # =======================================================================
    w_min = 0.5e-6
    w_max = 6.0e-6
    for index_scene in range(nb_points):
        # =======================================================================
        # Preparing the output folder
        # =======================================================================
        output_name = '{:s}-{:02d}'.format(output_prefix, index_scene+1)
        print("# [{:02d}] Preparing the output folder: {:s}".format(index_scene+1, output_name))
        full_output_path = os.path.join(output_path, output_name)
        try:
            os.makedirs(full_output_path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
            else:
                print("# CAUTION: the output folder already exists. Contents will be overwritten.")

        # =======================================================================
        # Creating the object source file
        # =======================================================================
        print("# [{:02d}] Creating the object source file.".format(index_scene+1))
        input_plane = 'SKY'
        object_type = 'ES'
        current_continuum = c_continuum.Continuum()
        current_continuum.m_set_input_plane(input_plane)
        valid_indices = numpy.where((wavelength >= w_min) & (wavelength <= w_max))
        current_continuum.m_set_wave(wavelength[valid_indices])
        for index_x in range(n_x):
            pos_x = start_x + index_x * step_x
            for index_y in range(n_y):
                pos_y = start_y + index_y * step_y
                current_continuum.m_add_point(pos_x, pos_y, cube[index_x, index_y,:][valid_indices],
                                              flag[valid_indices], ew[valid_indices])
        current_object = c_object.Object()
        reference = input_file
        author = '{:s} version {:s}'.format(argv[0], _version)
        description = 'user-provided data cube'
        current_object.m_set_header(reference, author, description)
        current_object.m_set_parameters(input_plane, object_type)
        current_object.m_add_continuum(current_continuum, overwrite=False)
        filename = os.path.join(output_path, output_name, 'object.fits')
        current_object.m_write_to_fits(filename)

        # =======================================================================
        # Preparing the list of objects
        # =======================================================================
        print("# [{:02d}] Creating the list of objects.".format(index_scene+1))
        reference = input_file
        author = '{:s} version {:s}'.format(argv[0], _version)
        description = 'Extended source in the IFU field of view - 4-point-dither - {:d} position.'.format(index_scene+1)
        object_list = c_objectlist.ObjectList()
        object_list.m_set_header(reference, author, description)
        object_list.m_set_input_plane(input_plane)
        object_list.m_add_object(gridx_v2v3[index_scene] * 3600., gridy_v2v3[index_scene] * 3600., 0.0, 1.0, 'ES',
                                 'IFU 4-point-nod - #{:d}'.format(index_scene+1), 'object.fits')

        filename = os.path.join(output_path, output_name, '{:s}.list'.format(output_name))
        object_list.m_write_to_fits(filename)

        # =======================================================================
        # Generating the figure in the V2V3 plane
        # =======================================================================
        print("# Generating the figures in the V2V3 plane.")
        work_x = start_x + numpy.arange(n_x) * step_x
        work_y = start_y + numpy.arange(n_y) * step_y
        intensity = numpy.ravel(numpy.sum(cube, axis=2))
        intensity /= numpy.mean(intensity)
        grid_y, grid_x = numpy.meshgrid(work_y, work_x)
        vector_x = numpy.ravel(grid_x) / 3600 + gridx_v2v3[index_scene]
        vector_y = numpy.ravel(grid_y) / 3600 + gridy_v2v3[index_scene]

        figure_name = 'figure_V2V3.pdf'
        figure = plt.figure()
        axes = figure.gca()
        axes.set_aspect('equal')
        # siaf.m_display_aperture('NRS_FULL_IFU', axes, plane='V2V3', corners=True, reference=False, angle=False,
        #                         color='g', linewidth=4)
        projtools.f_display_slices(OTE, NIRSpecIFU, axes, wavelength=2.5e-6, plane='V2V3', arcsec=False, color='k')
        image = axes.scatter(vector_x, vector_y, c=intensity, s=30, alpha=0.5)
        axes.plot(gridx_v2v3, gridy_v2v3, linestyle='None', marker='o', color='r',
                  label='dither points', markersize=1.0)
        axes.plot(gridx_v2v3[index_scene], gridy_v2v3[index_scene], linestyle='None', marker='o', color='y',
                  label='source position', markersize=0.5)

        # axes.plot(3600 * xcenter_v2v3, 3600 * ycenter_v2v3, linestyle='None', marker='x', color='b', label='center_aperture')
        cbar = plt.colorbar(image, extend='neither', spacing='proportional', orientation='vertical', shrink=0.9,
                            format="%3.2f")
        cbar.set_label("mask value", size=8)
        cbar.ax.tick_params(labelsize=8)
        axes.set_xlabel('x-axis position (in degrees, V2 plane)')
        axes.set_ylabel('y-axis position (in degrees, V3 plane)')
        for label in (axes.get_xticklabels() + axes.get_yticklabels()):
            label.set_fontname('Arial')
            label.set_fontsize(8)

        axes.set_title("4-point-nod IFU dither position - position #{:d} - {:s}".format(index_scene+1, output_name))
        axes.legend(prop={'size':6})
        plt.grid(True)
        work_text = 'V2V3 coordinates (arcsec): ( {:7.3f} , {:7.3f} )'.format(3600. * gridx_v2v3[index_scene],
                                                                              3600. * gridy_v2v3[index_scene])
        plt.figtext(0.01,0.01, work_text, ha='left', style='italic', size=6)
        plt.figtext(0.99,0.01,(datetime.datetime.now()).isoformat(), ha='right', style='italic', size=6)
        plt.savefig(os.path.join(output_path, output_name, figure_name))
        plt.close(figure)

