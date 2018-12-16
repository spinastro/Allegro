#!/usr/bin/env python

#################### READ FITS HEADERS & MAKE A LIST OF SPECTRA
import os
import glob
import sys
from astropy.io import fits
import numpy as np
from astropy.io import ascii
from pyraf import iraf
from astropy.table import Table

iraf.noao()
iraf.rv()
iraf.oned()

def selection(data):
    filt = (data['TWILIGHT'] == 'True') & (data['CCF_RVC'] == data['CCF_RVC']) & (data['DVRMS'] == data['DVRMS']) & (data['CCF_RVC'] != 'INDEF') & (data['DVRMS'] != 'INDEF') & (data['DVRMS'] < 10.) & (data['sat_pixls'] < lim_sat_pixels)
    num_filt=len(data['SN65'][filt])
    if num_filt > 2:
        mean_SN65=np.mean(data['SN65'][filt])
        std_SN65=np.std(data['SN65'][filt])
        low_lim_SN65=mean_SN65-sigmas*std_SN65
        w_good = filt & (data['SN65'] > low_lim_SN65)
        red_num = len(data['SN65'][w_good])
    if num_filt == 2:
        w_good = filt & (data['SN65'] > 50.)
        if len(data['SN65'][w_good]) == 2:
            red_num=len(data['SN65'][w_good])
        if len(data['SN65'][w_good]) < 2:
            SN_max=np.max(data['SN65'][filt])
            if SN_max > 30:
                w_good = filt & (data['SN65'] > SN_max-SN_max/3.)
                red_num = len(data['SN65'][w_good])
            else:
                print 'S/N < 30 - all spectra rejected'
                red_num = 0
    if num_filt == 1:
        if data['SN65'][filt] > 30:
            w_good= filt & (data['SN65'] > 30)
            red_num = len(data['SN65'][w_good])
        else:
            print 'S/N < 30 - all spectra rejected'
            red_num=0
    if num_filt == 0:
        print 'Can not read the data file'
        red_num=0
        
    if red_num > 0:
        return data[w_good]
    if red_num == 0:
        return False

if __name__ == "__main__":

    sigmas=1.
    lim_sat_pixels=1000

    star_folders=glob.glob('../*')
    star_folders.remove('../pipeline')

    print 'Number of stars', len(star_folders)

    for i in np.arange(0,len(star_folders)):
        name_star=star_folders[i][3:len(star_folders[i])]
        print name_star
        check_file=glob.glob(star_folders[i]+"/"+name_star+"_prered.dat")
        if len(check_file) == 0:
            data_file=glob.glob(star_folders[i]+"/*info.csv")
            data = ascii.read(data_file[0])
            print i,'star. Selection in progress. '+str(len(data))+' spectra.'
            analyse = selection(data)
        
            if analyse:
                a=[]
                for h in np.arange(len(analyse["CCF_RVC"])):
                    a.append(np.float(analyse["CCF_RVC"][h]))
                analyse["CCF_RVC"]=a

                w = np.where(np.abs(analyse["CCF_RVC"]) > 1000.)
                num=len(analyse["CCF_RVC"][w[0]])
                if num > 0:
                    for j in w[0]:
                        print analyse['NAME_S1D'][j], 'I am calculating the RV shift'
                        iraf.fxcor(star_folders[i]+"/untar/"+analyse['NAME_S1D'][j], "Sun.fits",output=star_folders[i]+"/untar/RVs",interactive="no", observatory="esovlt")
                        RVel=Table.read(star_folders[i]+"/untar/RVs.txt", format='ascii')
                        if RVel['col12'][0] == 'INDEF':
                            RVel['col12'][0]=-9999
                            RVel['col14'][0]=-9999
                        analyse["CCF_RVC"][j]=RVel['col12'][0]
                        analyse["FWHM"][j]=np.nan
                        analyse["DVRMS"][j]=RVel['col14'][0]
                        os.remove(star_folders[i]+"/untar/""RVs.txt")
                        os.remove(star_folders[i]+"/untar/""RVs.gki")
                        os.remove(star_folders[i]+"/untar/""RVs.log")

                f = analyse["CCF_RVC"] != -9999
                analyse=analyse[f]
                m=np.median(analyse["CCF_RVC"])
                s=np.std(analyse["CCF_RVC"])
                if s > 0.5:
                    w = np.where(np.abs(analyse["CCF_RVC"] - m) > 2*s)
                    num=len(analyse["CCF_RVC"][w[0]])
                    if num > 0:
                        for j in w[0]:
                            print analyse['NAME_S1D'][j], 'I am calculating the RV shift'
                            iraf.fxcor(star_folders[i]+"/untar/"+analyse['NAME_S1D'][j], "Sun.fits",output=star_folders[i]+"/untar/RVs",interactive="no", observatory="esovlt")
                            RVel=Table.read(star_folders[i]+"/untar/RVs.txt", format='ascii')
                            if RVel['col12'][0] == 'INDEF':
                                RVel['col12'][0]=-9999
                                RVel['col14'][0]=-9999
                            analyse["CCF_RVC"][j]=RVel['col12'][0]
                            analyse["FWHM"][j]=np.nan
                            analyse["DVRMS"][j]=RVel['col14'][0]
                            os.remove(star_folders[i]+"/untar/""RVs.txt")
                            os.remove(star_folders[i]+"/untar/""RVs.gki")
                            os.remove(star_folders[i]+"/untar/""RVs.log")

                        f = analyse["CCF_RVC"] != -9999
                        analyse=analyse[f]

                print 'Selected '+str(len(analyse["CCF_RVC"]))+' spectra.'
                ascii.write(analyse, star_folders[i]+'/'+name_star+"_prered.dat", format='csv')
