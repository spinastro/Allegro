#!/usr/bin/env python

#################### READ FITS HEADERS & MAKE A LIST OF SPECTRA
import glob
import tarfile
import os
import sys
import shutil
from astropy.io import fits
import numpy as np
from astropy.io import ascii
from astroplan import Observer
import astropy.coordinates
from pytz import timezone
from astropy.time import Time
from astropy import units as u
from pyraf import iraf
from astropy.table import Table
import time
import pandas as pd


iraf.noao()
iraf.rv()
iraf.oned()

LaSilla=Observer.at_site("lasilla", timezone="America/Santiago")

def info_s1d(s1d):
    with fits.open(s1d) as image:
        t_MJD=Time(image[0].header['MJD-OBS'], format='mjd')
        info = pd.DataFrame({"NAME_S1D":[s1d[s1d.find('untar/')+6:len(s1d)]],
                        "MJD":[image[0].header['MJD-OBS']],
                        "EXPTIME":[image[0].header['EXPTIME']],
                        "RA":[image[0].header['RA']],
                        "DEC":[image[0].header['DEC']],
                        "RADVEL":[image[0].header['HIERARCH ESO TEL TARG RADVEL']],
                        "NAME":[image[0].header['HIERARCH ESO OBS TARG NAME']],
                        "PROG_ID":[image[0].header['HIERARCH ESO OBS PROG ID']],
                        "OBS_START":[image[0].header['HIERARCH ESO OBS START']],
                        "OBS_TYPE":[image[0].header['HIERARCH ESO DPR TYPE']],
                        "SEEING":[image[0].header['HIERARCH ESO TEL AMBI FWHM START']],
                        "AIRMASS":[image[0].header['HIERARCH ESO TEL AIRM START']],
                        "SN65":[image[0].header['HIERARCH ESO DRS SPE EXT SN65']],
                        "TWILIGHT":[LaSilla.is_night(t_MJD, -12.*u.deg)]})
    return info

def info_we2ds(we2ds):
    if len(we2ds) > 0:
        image=fits.open(we2ds)
        flux=image[0].data[65]
        wsaturated=np.where(flux > flux_limit)
        info = pd.DataFrame({"sat_pixls":[len(wsaturated[0])]})
    else:
        info = pd.DataFrame({"sat_pixls":[np.nan]})
    return info

def info_bis(wbis):
    if len(wbis) > 0:
        with fits.open(wbis) as image:
            info = pd.DataFrame({"BJD":[image[0].header['HIERARCH ESO DRS BJD']],
                                "CCF_RVC":[image[0].header['HIERARCH ESO DRS CCF RVC']],
                                "FWHM":[image[0].header['HIERARCH ESO DRS CCF RVC']],
                                "MASK":[image[0].header['HIERARCH ESO DRS CCF MASK']]})
            key_check=image[0].header.get('HIERARCH ESO DRS DVRMS')
            if key_check == None:
                info["DVRMS"]=np.nan
            else:
                info["DVRMS"]=image[0].header['HIERARCH ESO DRS DVRMS'] #Estimated RV uncertainty [m/s]
    else:
        a=wbis.find('untar/')
        b=wbis.find('bis')
        iraf.fxcor(wbis[0:b]+"s1d_A.fits", "Sun.fits",output=wbis[0:a+6]+"RVs",interactive="no", observatory="esovlt")
        RVel=Table.read(wbis[0:a+6]+"RVs.txt", format='ascii')
        info = pd.DataFrame({"BJD":[np.nan],
                                "CCF_RVC":[RVel['col12'][0]],
                                "FWHM":[np.nan],
                                "MASK":[np.nan],
                                "DVRMS":[RVel['col14'][0]]})
        os.remove(wbis[0:a+6]+"RVs.txt")
        os.remove(wbis[0:a+6]+"RVs.gki")
        os.remove(wbis[0:a+6]+"RVs.log")
    return info




if __name__ == "__main__":

    flux_limit=150000

    star_folders=glob.glob('../*')
    star_folders.remove('../pipeline')

    print 'Number of stars', len(star_folders)

    for i in np.arange(len(star_folders)):
        star=star_folders[i][3:len(star_folders[i])]
        print star
        check_file=glob.glob(star_folders[i]+"/*info.csv")
        if len(check_file) == 0:
        
            s1d=glob.glob(star_folders[i]+"/untar/*s1d_A.fits")
            nfiles=len(s1d)
            arr_s1d = pd.DataFrame({})
            arr_s1d=arr_s1d.append([info_s1d(s1d[k]) for k in np.arange(0,len(s1d))], ignore_index=True)
        
            arr_we2ds = pd.DataFrame({})
            arr_we2ds=arr_we2ds.append([info_we2ds(glob.glob(s1d[k][0:len(star_folders[i]+'/untar/')+30]+"e2ds_A.fits")[0]) for k in np.arange(0,len(s1d))], ignore_index=True)

            arr_bis = pd.DataFrame({})
            arr_bis=arr_bis.append([info_bis(glob.glob(s1d[k][0:len(star_folders[i]+'/untar/')+30]+"bis_*_A.fits")[0]) for k in np.arange(0,len(s1d))], ignore_index=True)

            frames = [arr_s1d, arr_we2ds, arr_bis]
            result = pd.concat(frames, axis=1)
            result.to_csv(star_folders[i]+'/'+star+'_info.csv',index=False)
        else:
            print "Info already extracted"



