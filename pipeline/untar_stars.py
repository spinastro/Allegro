#!/usr/bin/env python

#################### READ FITS HEADERS & MAKE A LIST OF SPECTRA
import glob
import tarfile
import os
import sys
import shutil
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii

star_folders=glob.glob("../*")
star_folders.remove('../pipeline')

print 'Number of stars', len(star_folders)

for i in np.arange(len(star_folders)):
    check=glob.glob(star_folders[i]+'/untar')
    if len(check) == 0:
        print star_folders[i][3:len(star_folders[i])]
        tar_files=glob.glob(star_folders[i]+'/*tar')
        fits_files=glob.glob(star_folders[i]+'/*fits')
        print 'Number of tar', len(tar_files)
        print 'Number of fits', len(fits_files)
        os.mkdir(star_folders[i]+'/archive')
        os.mkdir(star_folders[i]+'/untar')
        os.mkdir(star_folders[i]+'/reduction')
        for k in np.arange(len(tar_files)):
            untar = tarfile.open(name=tar_files[k])
            untar.extractall(path=star_folders[i]+'/untar/')
            untar.close()
            shutil.move(tar_files[k],star_folders[i]+'/archive/')
        for k in np.arange(len(fits_files)):
            shutil.move(fits_files[k],star_folders[i]+'/archive/')
        nights=glob.glob(star_folders[i]+"/untar/data/reduced/*")
        if len(nights) > 0:
            for f in np.arange(len(nights)):
                files_to_move=glob.glob(nights[f]+'/*')
                for ff in np.arange(len(files_to_move)):
                    shutil.move(files_to_move[ff],star_folders[i]+'/untar/')
                os.rmdir(nights[f])
            os.rmdir(star_folders[i]+"/untar/data/reduced")
            os.rmdir(star_folders[i]+"/untar/data")
        files_to_move2=glob.glob(star_folders[i]+"/archive/HARPS*")
        s1d=glob.glob(star_folders[i]+"/untar/*s1d_A.fits")
        ccf=glob.glob(star_folders[i]+"/untar/*ccf_*_A.fits")
        bis=glob.glob(star_folders[i]+"/untar/*bis_*_A.fits")
        e2ds=glob.glob(star_folders[i]+"/untar/*e2ds_*_A.fits")
    else:
        print star_folders[i]+' already untarred'

