#!/usr/bin/env python
import sys
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import os
import glob
from astropy.io import fits
import numpy as np
from astropy.io import ascii

#ipython specvis.py star_name num_specs_per_plot

spectra=glob.glob("../"+str(sys.argv[1])+"/reduction/*_red_dn.fits")
num_spec=len(spectra)
num_spec_per_plot=int(sys.argv[2])

print 'Plot star:', str(sys.argv[1]), 'Num spec tot:', str(num_spec) , '; Num spec per plot:', str(sys.argv[2])

names=[]

interval=[65000,67000]

if num_spec < num_spec_per_plot:
    matrix_spectra=[]
    for i in np.arange(num_spec):
        sp=fits.open(spectra[i])
        names.append(spectra[i])
        flux=sp[0].data
        CRVAL1 = sp[0].header['CRVAL1']
        CDELT1 = sp[0].header['CDELT1']
        CRPIX1 = sp[0].header['CRPIX1']
        small_flux=flux[interval[0]:interval[1]]
        diff_wave=np.zeros(len(small_flux))
        for z in np.arange(len(small_flux)):
            diff_wave[z]=(CRVAL1+((-1)*CRPIX1+z+interval[0])*CDELT1) - 6000.
        w_start=np.where(np.abs(diff_wave) == min(np.abs(diff_wave)))[0]
        w_start=w_start[0]
        w_end=w_start+3000
        matrix_spectra.append(small_flux[w_start:w_end])
        wave=diff_wave[w_start:w_end] + 6000.
        plt.plot(wave,matrix_spectra[i], alpha=0.5, linewidth=0.5)
    plt.savefig(str(sys.argv[1])+'_plot_1.pdf')
    plt.clf()

else:
    n_plots=num_spec/num_spec_per_plot
    for k in np.arange(n_plots):
        print 'Plot', str(int(k+1))
        names.append('Plot '+str(int(k+1)))
        for i in np.arange(num_spec_per_plot):
            sp=fits.open(spectra[k*num_spec_per_plot+i])
            names.append(spectra[k*num_spec_per_plot+i])
            print spectra[k*num_spec_per_plot+i]
            flux=sp[0].data
            CRVAL1 = sp[0].header['CRVAL1']
            CDELT1 = sp[0].header['CDELT1']
            CRPIX1 = sp[0].header['CRPIX1']
            small_flux=flux[interval[0]:interval[1]]
            diff_wave=np.zeros(len(small_flux))
            for z in np.arange(len(small_flux)):
                diff_wave[z]=(CRVAL1+((-1)*CRPIX1+z+interval[0])*CDELT1) - 6000.
            w_start=np.where(np.abs(diff_wave) == min(np.abs(diff_wave)))[0]
            w_start=w_start[0]
            w_end=w_start+3000
            matrix_spectra=small_flux[w_start:w_end]
            wave=diff_wave[w_start:w_end] + 6000.
            plt.plot(wave,matrix_spectra, alpha=0.5, linewidth=0.5)
        plt.savefig(str(sys.argv[1])+'_plot_'+str(k)+'.pdf')
        plt.clf()

    spec_resid=num_spec-(n_plots)*num_spec_per_plot
    matrix_spectra=[]
    print 'Plot', str(int(n_plots+1))
    names.append('Plot '+str(int(n_plots+1)))
    for s in np.arange(spec_resid):
        sp=fits.open(spectra[spec_resid+1+s])
        print spectra[spec_resid+1+s]
        names.append(spectra[spec_resid+1+s])
        flux=sp[0].data
        CRVAL1 = sp[0].header['CRVAL1']
        CDELT1 = sp[0].header['CDELT1']
        CRPIX1 = sp[0].header['CRPIX1']
        small_flux=flux[interval[0]:interval[1]]
        diff_wave=np.zeros(len(small_flux))
        for z in np.arange(len(small_flux)):
            diff_wave[z]=(CRVAL1+((-1)*CRPIX1+z+interval[0])*CDELT1) - 6000.
        w_start=np.where(np.abs(diff_wave) == min(np.abs(diff_wave)))[0]
        w_start=w_start[0]
        w_end=w_start+3000
        matrix_spectra.append(small_flux[w_start:w_end])
        wave=diff_wave[w_start:w_end] + 6000.
        plt.plot(matrix_spectra[s], alpha=0.5, linewidth=0.5)
    plt.savefig(str(sys.argv[1])+'_plot_'+str(k+1)+'.pdf')
    plt.clf()




