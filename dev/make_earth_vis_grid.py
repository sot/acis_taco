# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Functions to make the earth_vis_grid files used in the ACIS FP model.

Derived from notebook:
https://nbviewer.jupyter.org/url/cxc.cfa.harvard.edu/mta/ASPECT/ipynb/xija/make-vis-grid.ipynb
"""

import sys
import os
import time
from pathlib import Path

import numpy as np
import astropy_healpix
from astropy.coordinates import SphericalRepresentation
import astropy.units as u
from astropy.io import fits

sys.path.insert(0, str(Path(os.environ['HOME'], 'git', 'acis_taco')))

# Standard gridding
alt_min = 200 * u.km
alt_max = 300000 * u.km
n_alt = 100
alts = np.logspace(np.log10(alt_min.to_value(u.m)),
                   np.log10(alt_max.to_value(u.m)),
                   n_alt)


def make_earth_vis_grids(nside=32, n_reps=1):
    from acis_taco import calc_earth_vis, RAD_EARTH, acis_taco
    hp = astropy_healpix.HEALPix(nside=nside, order='nested')
    npix = astropy_healpix.nside_to_npix(nside)
    print(f'npix={npix}')
    lons, lats = hp.healpix_to_lonlat(np.arange(npix))
    time0 = time.time()

    # Allow randomization between altitudes
    for i_rep in range(n_reps):
        vis_arrays = []
        # Randomize the ray-trace points for each rep
        acis_taco._RANDOM_SALT = None
        acis_taco.SPHERE_XYZ = acis_taco.random_hemisphere(acis_taco.N_SPHERE)
        acis_taco.SPHERE_X = acis_taco.SPHERE_XYZ[:, 0].copy()

        for i_alt, alt in enumerate(alts):

            srs = SphericalRepresentation(lon=lons, lat=lats, distance=alt + RAD_EARTH)
            xyzs = srs.to_cartesian()
            peb_xs = xyzs.x.to_value()
            peb_ys = xyzs.y.to_value()
            peb_zs = xyzs.z.to_value()
            vis = []
            for peb_x, peb_y, peb_z in zip(peb_xs, peb_ys, peb_zs):
                _, illums, _ = calc_earth_vis(p_earth_body=[peb_x, peb_y, peb_z])
                vis.append(np.sum(illums))
            vis = np.array(vis)
            vis_arrays.append(vis)
            vis_grid = np.vstack(vis_arrays)
            if i_alt % 10 == 0:
                print(f'alt={alt / 1000:.2f} km at dt={time.time() - time0:.1f}')

        ii = 1
        while True:
            filename = Path(f'earth_vis_grid_nside{nside}_rep{ii}.npy')
            if filename.exists():
                ii += 1
                continue
            else:
                print(f'Saving {filename}')
                np.save(filename, vis_grid)
                break

    return vis_grid


def make_earth_vis_grid_fits(nside=32):
    """Collect the reps and average and write to FITS"""
    from acis_taco import RAD_EARTH
    ii = 1
    vis_list = []
    while True:
        filename = Path(f'earth_vis_grid_nside{nside}_rep{ii}.npy')
        if filename.exists():
            # print(f'Reading {filename}')
            vis = np.load(filename)
            vis_list.append(vis)
            ii += 1
        else:
            break

    npix = astropy_healpix.nside_to_npix(nside)

    visl = np.array(vis_list)
    print(visl.shape)

    vis = visl.mean(axis=0)

    # Some sanity checking
    assert vis.shape == (100, npix)
    assert np.min(vis) == 0.0
    assert np.max(vis) < 3.0

    scale = 2**16 / 3.0
    visi = np.round(vis * scale)
    assert np.max(visi) < 2**16

    # Convert to 16 bits.  Also transpose so that interpolation along
    # distance is contiguous in memory.
    # BTW see: https://github.com/astropy/astropy/issues/8726 for the
    # origin of the .copy()
    visi2 = visi.astype(np.uint16).transpose().copy()

    hdu = fits.PrimaryHDU(visi2)
    hdu.header['nside'] = nside
    hdu.header['alt_min'] = alt_min.to_value(u.m)
    hdu.header['alt_max'] = alt_max.to_value(u.m)
    hdu.header['n_alt'] = n_alt
    hdu.header['scale'] = scale
    hdu.header['earthrad'] = RAD_EARTH

    filename = Path(f'earth_vis_grid_nside{nside}.fits.gz')
    print(f'Writing {filename}')
    hdul = fits.HDUList([hdu])
    hdul.writeto(filename, overwrite=True)
