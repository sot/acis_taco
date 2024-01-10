#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Calculate Earth illumination on the ACIS radiator over a specified interval
of time.
"""
from __future__ import print_function

import os
from datetime import date

import cheta.fetch_sci as fetch
import numpy as np
import six.moves.cPickle as pickle
import ska_quatutil
import ska_sun
from chandra_time import DateTime

import acis_taco


def get_options():
    from optparse import OptionParser
    parser = OptionParser()
    parser.set_defaults()
    parser.add_option("--start",
                      help="Start time (default=NOW")
    parser.add_option("--stop",
                      help="Stop time (default=<NWEEK> weekly outputs")
    parser.add_option("--ngrid",
                      type='int',
                      default=50,
                      help="Number of grid points (one-axis) on map")
    parser.add_option("--nweeks",
                      type='int',
                      default=4,
                      help="Number of future weeks (default=3)")
    parser.add_option("--out",
                      default="illum.pkl",
                      help="Output file name (default=illum.pkl)")
    parser.add_option("--nproc",
                      type='int',
                      default=1,
                      help="Number of processors to use")
    parser.add_option("--verbose",
                      action='store_true',
                      default=False,
                      help="Print verbose output")
    parser.add_option("--data-root",
                      default='.',
                      help="Data root (default='.')")

    (opt, args) = parser.parse_args()
    return (opt, args)

def calc_illums(queue, chandra_ecis, att):
    """Calculate the illumination (Earth solid angle) at specified ephemeris
    points ``chandra_ecis`` and attitude ``att``.  Put the result into the
    multiprocess ``queue``.
    """
    illums = []
    for chandra_eci in chandra_ecis:
        vis, illum, rays = acis_taco.calc_earth_vis(chandra_eci, att, max_reflect=5)
        illums.append(sum(illum))
    queue.put(illums)


def get_antisun_grid(ngrid=10):
    xy0 = ngrid / 2.0
    max_pitch = 135.0
    antisun = acis_taco.AntiSun(xy0, xy0, max_pitch * 2 / ngrid)
    x = np.arange(ngrid)[np.newaxis, :] + np.zeros((ngrid, ngrid))
    x = x.flatten()
    y = np.arange(ngrid)[:, np.newaxis] + np.zeros((ngrid, ngrid))
    y = y.flatten()
    return antisun, x, y

def calc_perigee_map(start='2010:114:20:00:00', stop='2010:117:16:00:00', ngrid=50):
    # Get orbital ephemeris in requested time range
    print('Fetching ephemeris')
    objs = ('orbit', 'lunar', 'solar')
    axes = ('x', 'y', 'z')
    msids = ['{0}ephem0_{1}'.format(obj, axis)
             for obj in objs
             for axis in axes]
    ephems = fetch.MSIDset(msids, start, stop)
    times = ephems[msids[0]].times.copy()
    n_ephem = len(ephems[msids[0]].vals)
    if any(len(ephems[msid].vals) != n_ephem for msid in msids):
        raise ValueError('Unexpected mismatch in ephemeris telemetry sampling')

    ephem_xyzs = {}
    for obj in ('orbit', 'lunar', 'solar'):
        msids = ('{0}ephem0_{1}'.format(obj, axis) for axis in axes)
        ephem_xyzs[obj] = np.array([ephems[msid].vals for msid in msids])

    illum_maps = []
    illum_idxs = []
    antisun, xs, ys = get_antisun_grid(ngrid)
    antisun_pitches, phis = antisun.img2polar(xs, ys)

    for i_ephem, ephem_xyz in enumerate(ephem_xyzs['orbit'].transpose()):
        # Calculate a grid of attitude vectors centered on the anti-sun line and extending
        # from sun-pitch=180 (center) to sun-pitch=45 (edge).  
        orbit_r = np.sqrt(np.sum(ephem_xyzs['orbit'][:, i_ephem]**2))
        if orbit_r > 50000e3:
            continue
        sun_eci = ephem_xyzs['solar'][:, i_ephem] - ephem_xyzs['orbit'][:, i_ephem]
        att_vecs = antisun.img2eci(xs, ys, sun_eci)
        ras, decs = ska_quatutil.eci2radec(att_vecs)
        illum_map = np.zeros((ngrid, ngrid), dtype=np.float32)
        print(i_ephem, n_ephem, ephem_xyz)
        for iy in range(ngrid):
            for ix in range(ngrid):
                i_vec = iy * ngrid + ix
                if antisun_pitches[i_vec] < 135:
                    ra, dec = ras[i_vec], decs[i_vec]
                    roll = ska_sun.nominal_roll(ra, dec, times[i_ephem])
                    _, att_illums, _ = acis_taco.calc_earth_vis(ephem_xyz, [ra, dec, roll], max_reflect=5)
                    illum = sum(att_illums)
                else:
                    illum = 0.0
                illum_map[iy, ix] = illum
        illum_idxs.append(i_ephem)
        illum_maps.append(illum_map)

    # debug_here()
    return dict(illums=np.array(illum_maps),
                illum_idxs=np.array(illum_idxs),
                times=times,
                antisun=antisun,
                ephem_xyzs=ephem_xyzs,
                )

def get_intervals(start, nweeks, data_root):
    """Get ``nweeks`` worth of 9 day-intervals which go from Sunday 12am to
    Tuesday 12am starting from the ``start`` date.  These nominally cover the
    time for weekly loads.  Only return intervals for weeks which do not already
    have an output file (i.e. don't re-compute existing weeks).
    """
    monday = start - date.fromtimestamp(start.unix).weekday()   # weekday() Monday=0 .. Sunday=6
    intervals = []
    for week in range(nweeks):
        start = monday.day_start() - 1
        stop = monday.day_start() + 8
        # CALDATE: YYYYMonDD at hh:mm:ss.ss..
        weekname = monday.caldate[4:9].upper() + monday.caldate[2:4]
        outfile = os.path.join(os.path.abspath(data_root),
                               weekname + '.pkl')
        if not os.path.exists(outfile):
            intervals.append((start, stop, outfile))
        else:
            print(weekname, 'already available - skipping')
        monday += 7

    return intervals

if __name__ == '__main__':
    opt, args = get_options()

    start = DateTime(opt.start)
    if opt.stop is None:
        intervals = get_intervals(start, opt.nweeks, opt.data_root)
    else:
        intervals = [(start, opt.stop, opt.out)]

    for start, stop, outfile in intervals:
        print('Calculating illums for', outfile)
        out = calc_perigee_map(start, stop, ngrid=opt.ngrid)
        pickle.dump(out, open(outfile, 'wb'), protocol=2)
