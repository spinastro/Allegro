

import numpy as np
from glob import glob
from astropy.io import fits
#from astropy.stats import sigma_clipped_stats
#from astropy.stats import sigma_clipping
from time import time
from tempfile import NamedTemporaryFile
import os.path

from spectres import spectres
from tqdm import tqdm

DEFAULT_MEMMAP_DIR = "../"


def linear_dispersion_limits(metadata):
    """
    Return the minimum dispersion, the maximum dispersion, the dispersion
    spacing, and the number of pixels.
    """

    cdelt = metadata["CDELT1"]
    naxis = metadata["NAXIS1"]
    crval = metadata["CRVAL1"]
    crpix = metadata.get("CRPIX1", 0)

    min_dispersion = crval - (crpix - 1) * cdelt 
    max_dispersion = crval + (naxis - 1 - (crpix - 1)) * cdelt

    return (min_dispersion, max_dispersion, cdelt, naxis)


def linear_dispersion_map(metadata):
    min_dispersion, max_dispersion, delta, P = linear_dispersion_limits(metadata)
    return np.arange(min_dispersion, max_dispersion + delta, delta)[:P]


def sigma_clip(data, sigma, method="mean", sigma_lower=None, sigma_upper=None, 
               iters=5, axis=0):

    if method not in ("mean", "median"):
        raise ValueError("method unknown")

    sigma_lower = sigma_lower if sigma_lower is not None else sigma
    sigma_upper = sigma_upper if sigma_upper is not None else sigma
    
    #assert np.all(np.isfinite(data))

    for i in tqdm(range(iters), desc="Doing sigma clipping"):

        max_value = np.nanmedian(data, axis=axis)
        std = np.nanstd(data, axis=axis)
        min_value = max_value - std * sigma_lower
        max_value += std * sigma_upper

        # Set masked values to nan.
        data[(data > max_value) + (data < min_value)] = np.nan

    f = np.nanmedian if method == "median" else np.nanmean

    return f(data, axis=axis)


def scombine(paths, method="median", sigma=3, iters=5, resampling="linear",
             memmap=False, memmap_dir=None):

    tick = time()

    available_methods = ("mean", "median")
    if method not in available_methods:
        raise ValueError(f"method must be one of {available_methods}")

    available_resampling = ("linear", "spectres")
    if resampling not in available_resampling:
        raise ValueError(f"resampling must be one of {available_resampling}")

    N = len(paths)

    # Re-bin to common wavelengths.
    # (min_wavelength, max_wavelength, delta_wavelength, N_pixels)
    dispersion_values = np.zeros((N, 4))

    for i, path in enumerate(tqdm(paths, desc="Loading dispersion maps")):

        with fits.open(path) as image:
            dispersion_values[i] = linear_dispersion_limits(image[0].header)

    # Define common wavelength scale.
    common_min = np.max(dispersion_values[:, 0])
    common_max = np.min(dispersion_values[:, 1])
    common_delta = np.max(dispersion_values[:, 2])
    resampled_dispersion = np.arange(common_min, common_max, common_delta)[1:-1]

    P = resampled_dispersion.size

    resampling_kwds = dict()
    resample = np.interp if resampling == "linear" else spectres

    if resampling == "linear":
        resampling_kwds.update(left=np.nan, right=np.nan)    

    if memmap:
        if memmap_dir is None:
            memmap_dir = DEFAULT_MEMMAP_DIR 

        fn = NamedTemporaryFile(dir=memmap_dir)
        memmap_kwds = dict(filename=fn, dtype=float, shape=(N, P))
        print(f"Creating memory mapped array with kwds in {memmap_dir} ({fn.name}): {memmap_kwds}")

        resampled_fluxes = np.memmap(**memmap_kwds)
    else:
        resampled_fluxes = np.nan * np.ones((N, P))

    for i, path in enumerate(tqdm(paths, desc="Rebinning spectra")):
        with fits.open(path) as image:
            dispersion = linear_dispersion_map(image[0].header)
            flux = image[0].data

            resampled_fluxes[i] = resample(resampled_dispersion,
                                           dispersion, flux,
                                           **resampling_kwds)

    # Combine.
    resampled_flux = sigma_clip(resampled_fluxes, sigma=sigma, method=method, 
                                iters=iters, axis=0)

    image = fits.open(paths[0])
    image[0].data = resampled_flux
    image[0].header.update(CRVAL1=resampled_dispersion[0],
                           CDELT1=common_delta,
                           CRPIX1=1,
                           NAXIS1=resampled_dispersion.size)

    tock = time()
    print(f"Finito after {tock - tick:.0f} seconds")

    return image







if __name__ == "__main__":

    star_folders=glob('../*')
    star_folders.remove('../pipeline')

    print('Number of stars', len(star_folders))

    for i in np.arange(len(star_folders)):
        name_star=star_folders[i][3:len(star_folders[i])]
        print(i,"Star "+name_star)
        if os.path.isfile(star_folders[i]+"/reduction/"+name_star+"_red.fits"):                #### REDUCTION ALREADY DONE
            print("Star already reduced")
        else:                                   #### START SPECTRAL REDUCTION
            check_file=glob("../"+name_star+"/reduction/*_red_dn.fits")
            if len(check_file) > 0:
                red = scombine(glob("../"+name_star+"/reduction/*_red_dn.fits"),method="mean", sigma=10, iters=10, memmap=True)
                red[0].header=red[0].header[0:22]
                red[0].header['OBJECT']=name_star
                red.writeto("../"+name_star+"/reduction/"+name_star+"_red.fits")
                blue = scombine(glob("../"+name_star+"/reduction/*_blue_dn.fits"),method="mean", sigma=10, iters=10, memmap=True)
                blue[0].header=blue[0].header[0:22]
                blue[0].header['OBJECT']=name_star
                blue.writeto("../"+name_star+"/reduction/"+name_star+"_blue.fits")
            else:
                print("There are no files to combine")


