#!/usr/bin/env /proj/sot/ska/bin/python
"""
Calculate Earth illumination on the ACIS radiatiator over a specified interval
of time.
"""

import sys, os
import itertools
import re
import time

import numpy as np
from multiprocessing import Process, Queue

import Quaternion
import Ska.engarchive.fetch_sci as fetch
import Ska.Sun
import Ska.quatutil
from antisun import AntiSun

import taco2
import cPickle as pickle

from IPython.Debugger import Tracer; debug_here = Tracer()

def get_options():
    from optparse import OptionParser
    parser = OptionParser()
    parser.set_defaults()
    parser.add_option("--start",
                      default='2010:225:19:00:00',
                      help="Start time")
    parser.add_option("--stop",
                      default='2010:228:22:00:00',
                      help="Stop time")
    parser.add_option("--radzone-dur",
                      type='float',
                      default=12,
                      help="Rad zone duration (hours)")
    parser.add_option("--out",
                      default='out',
                      help="Output root name")
    parser.add_option("--nproc",
                      type='int',
                      default=1,
                      help="Number of processors to use")
    parser.add_option("--verbose",
                      action='store_true',
                      default=False,
                      help="Print verbose output")
    
    (opt, args) = parser.parse_args()
    return (opt, args)

def calc_illums(queue, chandra_ecis, att):
    """Calculate the illumination (Earth solid angle) at specified ephemeris
    points ``chandra_ecis`` and attitude ``att``.  Put the result into the
    multiprocess ``queue``.
    """
    illums = []
    for chandra_eci in chandra_ecis:
        vis, illum, rays = taco2.calc_earth_vis(chandra_eci, att, max_reflect=5)
        illums.append(sum(illum))
    queue.put(illums)


def get_antisun_grid(ngrid=10):
    xy0 = ngrid / 2.0
    max_pitch = 135.0
    antisun = AntiSun(xy0, xy0, max_pitch * 2 / ngrid)
    x = np.arange(ngrid)[np.newaxis, :] + np.zeros((ngrid, ngrid))
    x = x.flatten()
    y = np.arange(ngrid)[:, np.newaxis] + np.zeros((ngrid, ngrid))
    y = y.flatten()
    return antisun, x, y

def get_att_vecs(ngrid=100):
    # Pitch angle from anti-sun direction
    max_pitch = 135.0
    x = np.linspace(-max_pitch, max_pitch, ngrid)[np.newaxis, :] + np.zeros((ngrid, ngrid))
    x = x.flatten()
    y = np.linspace(-max_pitch, max_pitch, ngrid)[:, np.newaxis] + np.zeros((ngrid, ngrid))
    y = y.flatten()
    theta = np.radians(np.sqrt(x**2 + y**2))
    phi = np.arctan2(y, x)
    sin_theta = np.sin(theta)
    att_vecs = np.array([np.cos(theta),
                         np.sin(phi) * sin_theta,
                         -np.cos(phi) * sin_theta])
    return theta, phi, att_vecs, x, y

def calc_perigee_map(start='2010:116:00:00:00', stop='2010:118:00:00:00', radzone_dur=12, ngrid=100):
    # Get orbital ephemeris in requested time range
    print ('Fetching ephemeris')
    objs = ('orbit', 'lunar', 'solar')
    axes = ('x', 'y', 'z')
    msids = ['{0}ephem1_{1}'.format(obj, axis)
             for obj in objs
             for axis in axes]
    ephems = fetch.MSIDset(msids, start, stop)
    times = ephems[msids[0]].times.copy()
    if any(len(ephems[msid].vals) != len(ephems[msids[0]].vals) for msid in msids):
        raise ValueError('Unexpected mismatch in ephemeris telemetry sampling')

    ephem_xyzs = {}
    for obj in ('orbit', 'lunar', 'solar'):
        msids = ('{0}ephem1_{1}'.format(obj, axis) for axis in axes)
        ephem_xyzs[obj] = np.array([ephems[msid].vals for msid in msids])

    # Find index and time of perigee, then take only times within "radzone" interval
    orbit_rs = np.sqrt((ephem_xyzs['orbit']**2).sum(0))
    i_perigee = np.argmin(orbit_rs)
    time_perigee = times[i_perigee]
    ok = abs(times - time_perigee) < radzone_dur * 3600.0 / 2
    times = times[ok]
    for obj in objs:
        ephem_xyzs[obj] = ephem_xyzs[obj][:, ok]
    n_ephem = len(times)
    
    # Get the sun position at perigee
    sun_ra, sun_dec = Ska.Sun.position(time_perigee)
    sun_eci = Ska.quatutil.radec2eci(sun_ra, sun_dec)

    # Calculate a grid of attitude vectors centered on the anti-sun line and extending
    # from sun-pitch=180 (center) to sun-pitch=45 (edge).  
    antisun, xs, ys = get_antisun_grid(ngrid)
    antisun_pitches, phis = antisun.img2polar(xs, ys)
    att_vecs = antisun.img2eci(xs, ys, sun_eci)
    ras, decs = Ska.quatutil.eci2radec(att_vecs)
    
    illums = np.zeros((n_ephem, ngrid, ngrid), dtype=np.float32)
    for i_ephem, ephem_xyz in enumerate(ephem_xyzs['orbit'].transpose()):
        print i_ephem, n_ephem, ephem_xyz
        for iy in range(ngrid):
            for ix in range(ngrid):
                i_vec = iy * ngrid + ix
                if antisun_pitches[i_vec] < 135:
                    ra, dec = ras[i_vec], decs[i_vec]
                    roll = Ska.Sun.nominal_roll(ra, dec, time_perigee)
                    _, att_illums, _ = taco2.calc_earth_vis(ephem_xyz, [ra, dec, roll], max_reflect=5)
                    illum = sum(att_illums)
                else:
                    illum = 0.0
                illums[i_ephem, iy, ix] = illum
        
    # debug_here()
    return dict(illums=illums,
                times=times,
                antisun=antisun,
                xs=xs,
                ys=ys,
                antisun_pitches=antisun_pitches,
                ras=ras,
                decs=decs,
                ephem_xyzs=ephem_xyzs,
                orbit_rs=orbit_rs,
                sun_eci=sun_eci,
                )

if 0 and __name__ == '__main__':
    opt, args = get_options()
    xs, ys, times, illums = calc_perigee_map(opt.start, opt.stop, opt.radzone_dur)
    
out = calc_perigee_map(ngrid=50,radzone_dur=4)
pickle.dump(out, open('cube3.pkl', 'w'))
