#!/usr/bin/env python
import sys
import os
from nrspylib.nirspec import nirspec
from nrspylib.ote import ote
from generate_IPS_scene import SceneGenerator

name = "generate_DS9regions"
_version = "1.0.0"


class RegionGenerator:
    def __init__(self, fwa, gwa):
            self.sc_gen = SceneGenerator()
            self.fwa = fwa
            self.gwa = gwa

    def generate(self, eMPT_out, nod_pos=0, folder_root = 'PROVA_'):

        sg = self.sc_gen
        output_folder = folder_root+'_n'+str(nod_pos)

        full_output_path = os.path.join(sg.output_path, output_folder)

        # eMPT out file for this scene
        filename_eMPTout = os.path.join(sg.input_path, eMPT_out)

        n, src_id, q, i, j, offset_x, offset_y = self.sc_gen.read_eMPToutput(filename_eMPTout, nod_pos=nod_pos)
        print('# Read ' + str(n.shape[0]) + ' entries from ' + eMPT_out)
        print('# First entry ', src_id[0], q[0], i[0], j[0])  # First entry
        print('# Last entry', src_id[-1], q[-1], i[-1], j[-1])  # Last entry


        print("# Initialising the geometrical models (OTE + NIRSpec MOS + NIRSpec IFU).")
        NIRSpec = nirspec.NIRSpec()
        NIRSpec.m_initializeMOS(self.fwa, self.gwa, path=sg.model_path, model=sg.model_name)
        print("# NIRSpec configuration: {:s}/{:s}".format(sg.FWA, sg.GWA))
        print("# Name of the  geometrical model of NIRSpec and the OTE: {:s}.".format(sg.model_name))

        msa = NIRSpec.MSA
        spectrograph = NIRSpec.spectrograph
        fpa = NIRSpec.FPA
        wavelength = 2.5e-6
        order = 0

        list_posSCA491 = []
        list_posSCA492 = []
        for k in range(len(n)):
            # print(q[k], i[k], j[k], offset_x[k], offset_y[k])
            posMSA = msa.m_positionMS(q[k], i[k], j[k], offsetx=offset_x[k], offsety=offset_y[k])
            posFPA = spectrograph.m_applyForwardTransformToGrid(posMSA, wavelength, order)
            posSCA_491 = fpa.m_findFractionalPixelsNoTest(posFPA[0], posFPA[1], 491)
            posSCA_492 = fpa.m_findFractionalPixelsNoTest(posFPA[0], posFPA[1], 492)
            list_posSCA491.append([posSCA_491, q[k], i[k], j[k], src_id[k]])
            list_posSCA492.append([posSCA_492, q[k], i[k], j[k], src_id[k]])



        filename1 = os.path.join(sg.output_path, output_folder, '{:s}_SCA491.reg'.format(output_folder))
        filename2 = os.path.join(sg.output_path, output_folder, '{:s}_SCA492.reg'.format(output_folder))

        color = 'green'
        radius = 2.0

        # Writing output ds9 regions
        f = open(filename1, 'w')

        f.write('# Region file format: DS9 version 4.1\n')
        f.write('# Filename: ' + filename1 + '\n')
        f.write(
            'global color=' + color + ' dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
        f.write('image\n')

        for src in list_posSCA491:
            cx = src[0][0][0] - 4  # because Preocessed maps are 2040*2040 in size (i.e. there are no ref.pixels)
            cy = src[0][1][0] - 4  # because Preocessed maps are 2040*2040 in size (i.e. there are no ref.pixels)
            # If you want to write regions for .erm file uncomment these:
            #cx= src[0][0][0]
            #cy= src[0][1][0]
            q = src[1]
            i = src[2]
            j = src[3]
            name = src[4]
            info_text = 'MSA=[' + str(q) + ' ' + str(i) + ' ' + str(j) + ']' + ' ' + name
            #print('# '+info_text)
            f.write('circle(' + str(cx) + ',' + str(cy) + ',' + str(radius) + ') # text={(' + info_text + ')}\n')

        f.close()
        print('# '+filename1+' saved')

        f = open(filename2, 'w')

        f.write('# Region file format: DS9 version 4.1\n')
        f.write('# Filename: ' + filename2 + '\n')
        f.write(
            'global color=' + color + ' dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
        f.write('image\n')

        for src in list_posSCA492:
            cx = (2048 - src[0][0][0]) - 4  # Because in processed maps SCA 492 is swapped (in x & y) and there are no ref.pixels
            cy = (2048 - src[0][1][0]) - 4  # Because in processed maps SCA 492 is swapped (in x & y) and there are no ref.pixels
            # If you want to write regions for .erm file uncomment these:
            #cx= src[0][0][0]
            #cy= src[0][1][0]
            q = src[1]
            i = src[2]
            j = src[3]
            name = src[4]
            info_text = 'MSA=[' + str(q) + ' ' + str(i) + ' ' + str(j) + ']' + ' ' + name
            #print('# ' + info_text)
            f.write('circle(' + str(cx) + ',' + str(cy) + ',' + str(radius) + ') # text={(' + info_text + ')}\n')

        f.close()

        print('# '+filename2+' saved')
        print('# Done')



if __name__ == '__main__':

    dithers_fname = ['output_dither_00.txt', 'output_dither_1a.txt', 'output_dither_2a.txt']
    # dithers_fname = ['output_dither_00.txt'] # for testing

    rg = RegionGenerator('CLEAR', 'PRISM')

    # Iterate over dither
    for dfile in dithers_fname:
        print('# Processing '+dfile)
        # Iterate over noddings - only region for nod=0 generted for the moment
        for i in range(0, 1):
            print('# Generate DS9 regions for mask '+dfile)
            rg.generate(dfile, nod_pos=i, folder_root=dfile[7:len(dfile)-4])
