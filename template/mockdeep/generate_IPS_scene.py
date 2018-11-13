#!/usr/bin/env python
import sys
import os
import errno
import datetime
import io
import shutil
import numpy as np
import json
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
from nrspysim.scenes import continuum
from nrspysim.scenes import specunres
from nrspysim.scenes import object
from nrspysim.scenes import objectlist
from nrspylib.nirspec import nirspec
from nrspylib.ote import ote
from nrspysim.misc import projtools
from nrspylib.msa import msaconfig
from nrspylib.msa import msaoperability
from nrspylib.msa import msaimage

name = "generate_IPS_scene"
_version = "1.0.0"

# Factor to convert Jacopo's SED's unig [ergs s-1 cm-2 A-1] to MKS i.e [J s-1 m-2 m-1] or [W m-2 m-1]
# 1.0e-7
_f_cgs2mks = 1.0e7

# Factor to go from cgs to mks (for emission lines)
_cgs2mks =  1.0e-3


class SceneGenerator:

    def __init__(self):
        # Configuration mock data - INPUT-OUTPUT folders
        basedir = '/Users/ggiardin/JWST/IPSWork/00MockDEEP/'
        self.input_path = os.path.join(basedir, 'data')
        self.output_path = os.path.join(basedir, 'test_scene')

        # Source spectrum from Beagle output file
        self.fname_spectrum = os.path.join(self.input_path, 'SF_and_quiescent_catalogue_mags.fits')
        jsonfile = os.path.join(self.input_path, 'cb2016_n2_mup300_N015_O01_deplO70_C100_Jan16_line_wavelengths_may2017.json')

        # Reference file of failed-shutters for ISIM-CV3 - currently used by Peter(eMPT) & MPT:
        self.inputFailureMap = 'nrs_msal_CHK_20151211.msl'

        # Configuration Instrument
        self.model_path = '/Users/ggiardin/JWST/Software/JWST_Python/data/IQLAC/'
        self.model_name = 'NIRS_FM2_05_CV3_FIT1'
        self.FWA = 'CLEAR'
        self.GWA = 'MIRROR'  # does not matter for scene creation

        # Opening json file that sepcify which (BEAGLE) lines to use
        with open(jsonfile) as data_file:
            data = json.loads(data_file.read())
            print('# JSON entries', len(data))

        # Building dictionary of lines to use
        self.line_dict = {}
        for i in range(len(data)):
            if(data[i]['use']):
                label = data[i]['label']
                wl = data[i]['wl']
                self.line_dict[label] = wl

    def read_eMPToutput(self, filename, nod_pos = 0):
        f = io.open(filename, 'r')
        print('Reading in '+filename)
        n = np.empty((0))
        src_id = np.empty((0))
        q = np.empty((0), dtype='int')
        i = np.empty((0), dtype='int')
        j = np.empty((0), dtype='int')
        offset_x = np.empty((0), dtype='float')
        offset_y = np.empty((0), dtype='float')

        # Finding the point in the  file where the relevant list start...
        for iline, line in enumerate(f):
            # Finding where to start reading in...
            if ('Accepted targets, assigned slitlets' in line):
                #print('# *** Found begin of list***')
                break

        # skipping header
        f.readline()
        f.readline()

        ind_q = 4+5*(nod_pos)
        ind_i = ind_q+1
        ind_j = ind_q+2
        ind_osx = ind_q+3
        ind_osy = ind_q+4

        # Reading placing of the sources in selected shutters:
        for iline, line in enumerate(f):
            # reading relevant output
            tks = line.split()
            if (len(tks) > 5):
                # print(tks[0], tks[1])
                n = np.append(n, tks[0])
                src_id = np.append(src_id, tks[2]) #We need ID in the orginal mock catalogue!
                q = np.append(q, int(tks[ind_q]))
                i = np.append(i, int(tks[ind_i]))
                j = np.append(j, int(tks[ind_j]))
                offset_x = np.append(offset_x, float(tks[ind_osx]))
                offset_y = np.append(offset_y, float(tks[ind_osy]))

            else:
                break
            # For debug: print(tks)

        return [n, src_id, q, i, j, offset_x, offset_y]

    def get_continuum_sed(self, filename, src_id):

        hdulist = pyfits.open(filename)
        src_ind = self.get_source_tabindex(hdulist, src_id)
        spec = self.read_Continuum(hdulist, src_ind)

        hdulist.close()
        return spec

    def get_source_tabindex(self, hdulist, src_id, verbose = False):
        #Get ID table
        ids = hdulist['IDs'].data['ID_cat']
        src_tabindex = np.where(ids == int(src_id))
        if(verbose):
            print('Source '+src_id+' found at index '+str(src_tabindex[0]))
        return src_tabindex[0]

    def read_Continuum(self, hdulist, src_tabindex, verbose=False):
        if (verbose):
            print('Getting spectra for source at index '+str(src_tabindex))

        # Get the wavelength array (units of Ang)
        wl = hdulist['continuum sed wl'].data['wl'][0, :]
        # converting wave Amstrong -> m
        wave = wl * 1e-10
        # Getting redshift information
        all_z = hdulist['galaxy properties'].data['redshift']
        z = all_z[src_tabindex]
        if (verbose):
            print('Source redshift is '+str(z))
        # Loading SED information from the image extension (units of erg s^-1 cm^-2 A^-1)
        fluxin = hdulist['continuum sed'].data[src_tabindex][0]

        # red-shfiting spectrum
        wave = (1 + z) * wave
        # Conserving the energy:
        flux_z = fluxin / (1 + z)
        # Converting from [ergs-1 s-1 cm-2 Amstrong-1] to [J s-1 m-2 m-1]
        flux = _f_cgs2mks * flux_z

        # Checking that there are no wavelength duplicate in input spectrum
        # Jacopo's spectra seems to contain some and IPS is fussy about that..
        previousWavelength = 0.0
        clean_wave = []
        clean_flux = []
        for i in range(wave.shape[0]):
            if (wave[i] == previousWavelength):
                if (verbose):
                        print("Found a duplication and skipping it (kept only the first occurence).")
                # do nothing
            else:
                clean_wave.append(wave[i])
                clean_flux.append(flux[i])

            previousWavelength = wave[i]

        c_wave = np.array(clean_wave)
        c_flux = np.array(clean_flux)

        if verbose:
            print('Read spectrum for MPT Target ' + str(src_tabindex))

        input_plane = 'SKY'
        spec = continuum.Continuum()
        spec.m_simple_point_source(input_plane, c_wave, c_flux, x_pos=0.0, y_pos=0.0)
        return spec

    def get_emi_lines(self, filename, src_id, verbose = False):

        hdulist = pyfits.open(filename)
        src_ind = self.get_source_tabindex(hdulist, src_id)
        #need to read emission lines
        if (verbose):
            print('Getting emission lines for source at index '+str(src_ind))


        # Getting redshift information
        all_z = hdulist['galaxy properties'].data['redshift']
        z = all_z[src_ind]
        if (verbose):
            print('Source redshift is ' + str(z))

        tab = hdulist['hii emission'].data
        wavelength = np.empty((0))
        flux = np.empty((0))

        for i, line in enumerate(self.line_dict):
            lname = line+'_flux'
            addline = True
            #getting line flux and converting from [ergs-1 s-1 cm-2] to [J s-1 m-2]
            try:
                l_flux = _cgs2mks * tab[lname][src_ind]
            except:
                print('# WARNING: line '+lname+' is not in fits file')
                addline = False
            if(addline):
                # here one could put a cut on l_flux for lines that are too faint...
                # It would make the memory inprint of the simulations smaller...
                flux = np.append(flux, l_flux)
                # getting line wavelength and converting Amstrong -> m
                wave = 1e-10 *self.line_dict[line]
                 # red-shfiting line
                wave = (1 + z) * wave
                wavelength = np.append(wavelength, wave)
        if(verbose):
            print('# Found '+str(wavelength.size)+' lines in FITS table')

        input_plane = 'SKY'
        emi_lines = specunres.SpecUnres()
        emi_lines.m_simple_point_source(input_plane, wavelength, flux, x_pos=0.0, y_pos=0.0)
        hdulist.close()
        return emi_lines

    def generate(self, eMPT_out, nod_pos=0, folder_root = 'PROVA_'):

        output_folder = folder_root+'_n'+str(nod_pos)

        print("# Preparing the output folder: "+output_folder)
        full_output_path = os.path.join(self.output_path, output_folder)
        try:
            os.makedirs(full_output_path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
            else:
                print("# CAUTION: the output folder already exists. Contents will be overwritten.")

        # eMPT out file for this scene
        filename_eMPTout = os.path.join(self.input_path, eMPT_out)

        n, src_id, q, i, j, offset_x, offset_y = self.read_eMPToutput(filename_eMPTout, nod_pos=nod_pos)
        print('# Read ' + str(n.shape[0]) + ' entries from ' + eMPT_out)
        print('# First entry ', src_id[0], q[0], i[0], j[0])  # First entry
        print('# Last entry', src_id[-1], q[-1], i[-1], j[-1])  # Last entry

        print("# Initialising the geometrical models (OTE + NIRSpec MOS + NIRSpec IFU).")
        NIRSpec = nirspec.NIRSpec()
        NIRSpec.m_initializeMOS(self.FWA, self.GWA, path=self.model_path, model=self.model_name)
        OTE = ote.OTE()
        OTE.m_initializeFromModel(self.model_path, self.model_name)
        print("# NIRSpec configuration: {:s}/{:s}".format(self.FWA, self.GWA))
        print("# Name of the  geometrical model of NIRSpec and the OTE: {:s}.".format(self.model_name))

        # Creating scene object
        reference = 'Scene for mock DEEP'
        description = 'Source placing as per eMPT mask ' + eMPT_out
        author = '{:s} version {:s}'.format(sys.argv[0], _version)
        object_list = objectlist.ObjectList()
        object_list.m_set_header(reference, author, description)

        input_plane = 'SKY'
        object_list.m_set_input_plane(input_plane)

        for k in range(50):
            contsp = self.get_continuum_sed(self.fname_spectrum, src_id[k])
            emi_lines = self.get_emi_lines(self.fname_spectrum, src_id[k])

            # Creating IPS spectrum for this source - Point Source
            object_type = 'PS'
            current_object = object.Object()
            current_object.m_set_parameters(input_plane, object_type)
            src_name = 'src_'+src_id[k]
            src_fitsfile = src_name+'.fits'
            current_object.m_set_header('Source spectrum', author, src_name)

            # Continuum part
            current_object.m_add_continuum(contsp, overwrite=False)

            # Emission lines
            current_object.m_add_specunres(emi_lines, overwrite=False)


            # Write file
            filename = os.path.join(self.output_path, output_folder , src_fitsfile)
            current_object.m_write_to_fits(filename)

            print('# Processing source '+src_id[k]+' - MOS: '+str(q[k])+' '+str(i[k])+' '+ str(j[k]))
            current_x, current_y = projtools.f_project_micro_shutter(OTE, NIRSpec,
                                                                     q[k], i[k], j[k], offset_x=offset_x[k], offset_y=offset_y[k],
                                                                     wavelength=2.5e-6)
            # Adding source to the IPS scene object
            object_list.m_add_object(3600 * current_x, 3600 * current_y, 0.0, 1.0, 'PS',
                                     src_id[k], src_fitsfile)

            # Generating the associated figure
            figure_name = src_id[k] + '_MOS_{:01d}_{:03d}_{:03d}.pdf'.format(q[k], i[k], j[k])
            figure = plt.figure()
            axes = figure.gca()
            axes.set_aspect('equal')
            # For nod_pos = 1, moving the slit let down
            if(nod_pos == 1):
                slit_offset = -1
            elif(nod_pos == 2):
                slit_offset = +1
            else:
                slit_offset=0
            projtools.f_display_mini_slit(OTE, NIRSpec, axes, q[k], i[k], j[k]+slit_offset, 1, wavelength=2.5e-6, plane='V2V3',
                                          patches=True, arcsec=True)
            ima = axes.plot(3600. * current_x, 3600. * current_y, linestyle='none', marker='o', color='r')
            axes.set_xlabel('x-axis position (in arcsec, V2V3 plane)')
            axes.set_ylabel('x-axis position (in arcsec, V2V3 plane)')
            plt.title(src_id[k]+"- Q{:01d} [{:03d}, {:03d}]".format(q[k], i[k], j[k]), fontsize=11)
            plt.grid(True)
            plt.figtext(0.99, 0.01, (datetime.datetime.now()).isoformat(), ha='right', style='italic', size=6)
            plt.savefig(os.path.join(self.output_path, output_folder, figure_name))
            plt.close(figure)


        filename_list = os.path.join(self.output_path, output_folder, '{:s}.list'.format(output_folder))
        print("# Creating the list of objects " + filename_list)
        object_list.m_write_to_fits(filename_list)

        # Generating the MSA configuration file
        ## NB THIS HAS TO BE DONE ONLY FOR ONE POINTING IN THE 3-NODDING!!!!
        if(nod_pos== 0):
            print("# Generating an MSA configuration file.")
            description = "MSA configuration file for scene {:s}".format(output_folder)
            author = '{:s} version {:s}'.format(sys.argv[0], _version)
            description = 'Q4[186,81] open - perfect micro-shutter array.'
            nq = 4
            nx = 365
            ny = 171
            configuration = msaconfig.MSAConfig()
            configuration.m_setParameters(description, author, reference, nq, nx, ny)
            configuration.m_setAll(0)
            for k in range(len(n)):
                configuration.m_setShutterOpen(q[k], i[k], j[k])
                if (j[k] < 171):
                    configuration.m_setShutterOpen(q[k], i[k], j[k] + 1)
                if (j[k] > 1):
                    configuration.m_setShutterOpen(q[k], i[k], j[k] - 1)

            # Here need to add failed-open!
            "# Opening input failed closed/opened map."
            flist = msaoperability.MSAOperability()
            print('# Opening msl file')
            flist.m_readFromFITS(os.path.join(self.input_path, self.inputFailureMap))
            print('# Opened file: ', self.inputFailureMap)

            for q in range(1, 5):
                fos = flist.m_getListOfFailedOpenShutters(q)
                configuration.m_setShuttersOpen(q, fos[0], fos[1])

            # Need to add field stop (one failed open is behind field-stop (do by hand?)
            mos_filename = os.path.join(self.output_path, output_folder, '{:s}.mos'.format(output_folder))
            configuration.m_writeToFITS(mos_filename)
            print('# Generated configuration file ' + mos_filename)

            # Generating the MSA configuration report
            image = msaimage.MSAImage()
            image.m_construct(configuration)
            pdf_filename = os.path.join(self.output_path, output_folder, '{:s}.mos.pdf'.format(output_folder))
            image.m_generateReportPDF(mos_filename, pdf_filename)
        else:
            # Need to copy the file MOS file generated for nod_pos=0 to the directory of
            # the scene of the other nod pointings
            print("# Copying MSA configuration file from nod position 0.")
            output_folder_n0 = folder_root+'_n0'
            mos_filename_n0 = os.path.join(self.output_path, output_folder_n0, '{:s}.mos'.format(output_folder_n0))
            mos_filename_ni = os.path.join(self.output_path, output_folder, '{:s}.mos'.format(output_folder_n0))
            shutil.copy(mos_filename_n0, mos_filename_ni)

    print('# Done')



if __name__ == '__main__':

    dithers_fname = ['output_dither_00.txt', 'output_dither_1a.txt', 'output_dither_2a.txt']
    #dithers_fname = ['output_dither_00.txt'] for testing

    sg = SceneGenerator()

    # Iterate over dither
    for dfile in dithers_fname[0:1]:
        print('# Processing '+dfile)
        # Iterate over noddings
        for i in range(0, 1):
            print('Generate IPS scene for mask '+dfile)
            sg.generate(dfile, nod_pos=i, folder_root=dfile[7:len(dfile)-4])
