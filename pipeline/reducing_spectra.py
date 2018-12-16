#!/usr/bin/env python

#################### READ FITS HEADERS & MAKE A LIST OF SPECTRA
from pyraf import iraf
import os
import glob

iraf.noao()
iraf.oned()
iraf.rv()

import sys
from astropy.io import fits
import numpy as np
from astropy.io import ascii
from tqdm import tqdm

sigmas=1.

star_folders=glob.glob('../*')
star_folders.remove('../pipeline')

print 'Number of stars', len(star_folders)

for i in np.arange(0,len(star_folders)):
    name_star=star_folders[i][3:len(star_folders[i])]
    check_file=glob.glob(star_folders[i]+"/reduction/blue.dat")
    print 'Star '+name_star, i
    if len(check_file) == 1:                #### REDUCTION ALREADY DONE
        print 'Star already reduced'
    else:                                   #### START SPECTRAL REDUCTION
        data_file=glob.glob(star_folders[i]+"/*_prered.dat")
        if len(data_file) > 0:
            analyse = ascii.read(data_file[0])
            if len(analyse) > 1:
                print star_folders[i]+"/reduction/"
                norm_blue=[]
                norm_red=[]
                file_name=[]
                for k in tqdm(np.arange(len(analyse['NAME_S1D'])), desc="Reducing "+name_star):
                    hdu=fits.open(star_folders[i]+"/untar/"+analyse['NAME_S1D'][k])
                    file_name.append(analyse['NAME_S1D'][k][0:len(analyse['NAME_S1D'][k])-5])
                    print k, star_folders[i]+"/reduction/"+file_name[k]
                    scidata = hdu[0].data
                    length=len(scidata)-1
                    w_zero=np.where(scidata ==0)
                    low_gap=min(w_zero[0]-1)
                    up_gap=max(w_zero[0]+1)
                    iraf.scopy(star_folders[i]+"/untar/"+analyse['NAME_S1D'][k]+"[0:"+str(low_gap)+"]",star_folders[i]+"/reduction/"+file_name[k]+"_blue.fits")
                    iraf.scopy(star_folders[i]+"/untar/"+analyse['NAME_S1D'][k]+"["+str(up_gap)+":"+str(length)+"]",star_folders[i]+"/reduction/"+file_name[k]+"_red.fits")
                    iraf.dopcor(star_folders[i]+"/reduction/"+file_name[k]+"_blue.fits",star_folders[i]+"/reduction/"+file_name[k]+"_blue_d.fits", analyse['CCF_RVC'][k],isveloc="yes")
                    iraf.dopcor(star_folders[i]+"/reduction/"+file_name[k]+"_red.fits",star_folders[i]+"/reduction/"+file_name[k]+"_red_d.fits", analyse['CCF_RVC'][k],isveloc="yes")
                    iraf.continuum(star_folders[i]+"/reduction/"+file_name[k]+"_blue_d.fits",star_folders[i]+"/reduction/"+file_name[k]+"_blue_dn.fits",order=80,nit=7,low_rej=1.0,high_rej=3,interactive ="no", wavescale ="no")
                    iraf.continuum(star_folders[i]+"/reduction/"+file_name[k]+"_red_d.fits",star_folders[i]+"/reduction/"+file_name[k]+"_red_dn.fits",order=80,nit=7,low_rej=1.0,high_rej=3,interactive ="no", wavescale ="no")
                    norm_blue.append(file_name[k]+"_blue_dn.fits")
                    norm_red.append(file_name[k]+"_red_dn.fits")
                ascii.write([norm_blue], star_folders[i]+"/reduction/"+"blue.dat",format='no_header',overwrite=True)
                ascii.write([norm_red], star_folders[i]+"/reduction/"+"red.dat",format='no_header',overwrite=True)
            if len(analyse) == 1:
                print star_folders[i]+"/reduction/"
                os.chdir(star_folders[i]+"/reduction/") #### CHANGE DIRECTORY (DONT FORGET!)
                norm_blue=[]
                norm_red=[]
                file_name=[]
                hdu=fits.open(star_folders[i]+"/untar/"+analyse['NAME_S1D'][0])
                file_name.append(analyse['NAME_S1D'][0][0:len(analyse['NAME_S1D'][0])-5])
                scidata = hdu[0].data
                length=len(scidata)-1
                w_zero=np.where(scidata ==0)
                low_gap=min(w_zero[0]-1)
                up_gap=max(w_zero[0]+1)
                print file_name[0], name_star
                iraf.scopy(star_folders[i]+"/untar/"+analyse['NAME_S1D'][0]+"[0:"+str(low_gap)+"]",star_folders[i]+"/reduction/"+file_name[0]+"_blue.fits")
                iraf.scopy(star_folders[i]+"/untar/"+analyse['NAME_S1D'][0]+"["+str(up_gap)+":"+str(length)+"]",star_folders[i]+"/reduction/"+file_name[0]+"_red.fits")
                iraf.dopcor(star_folders[i]+"/reduction/"+file_name[0]+"_blue.fits",star_folders[i]+"/reduction/"+file_name[0]+"_blue_d.fits", analyse['CCF_RVC'][0],isveloc="yes")
                iraf.dopcor(star_folders[i]+"/reduction/"+file_name[0]+"_red.fits",star_folders[i]+"/reduction/"+file_name[0]+"_red_d.fits", analyse['CCF_RVC'][0],isveloc="yes")
                iraf.continuum(star_folders[i]+"/reduction/"+file_name[0]+"_blue_d.fits",star_folders[i]+"/reduction/"+name_star+"_blue.fits",order=80,nit=7,low_rej=1.0,high_rej=3,interactive ="no", wavescale ="no")
                iraf.continuum(star_folders[i]+"/reduction/"+file_name[0]+"_red_d.fits",star_folders[i]+"/reduction/"+name_star+"_red.fits",order=80,nit=7,low_rej=1.0,high_rej=3,interactive ="no", wavescale ="no")
                ascii.write([norm_blue], star_folders[i]+"/reduction/"+"blue.dat",format='no_header',overwrite=True)
                ascii.write([norm_red], star_folders[i]+"/reduction/"+"red.dat",format='no_header',overwrite=True)

