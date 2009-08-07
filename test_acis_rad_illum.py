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

import Ska.Table
import Quaternion
from Chandra.Time import DateTime
from Ska.TelemArchive.fetch import fetch
import Ska.Numpy
import matplotlib.pyplot as plt
from Ska.Matplotlib import plot_cxctime

import test_taco as taco

def get_options():
    from optparse import OptionParser
    parser = OptionParser()
    parser.set_defaults()
    parser.add_option("--tstart",
                      default='2009:060:00:00:00',
                      help="Start time")
    parser.add_option("--tstop",
                      default='2009:060:08:00:00',
                      help="Stop time")
    parser.add_option("--ephemfile",
                      default='orbitf352123502N001_eph0.fits',
                      help="Orbit ephemeris file containing start and stop times")
    parser.add_option("--out",
                      default='out',
                      help="Output root name")
    parser.add_option("--sample",
                      type='int',
                      default=2,
                      help="Sample")
    parser.add_option("--ngrid",
                      type='int',
                      default=100,
                      help="Number of rays to illuminate Earth disk")
    parser.add_option("--nproc",
                      type='int',
                      default=1,
                      help="Number of processors to use")
    parser.add_option("--verbose",
                      action='store_true',
                      default=False,
                      help="Print verbose output")
    parser.add_option("--movie",
                      action='store_true',
                      default=False,
                      help="Create visibility images for making a movie")
    
    (opt, args) = parser.parse_args()
    return (opt, args)

def calc_vis_values(iproc, ephem_Times, chandra_ecis, q1s, q2s, q3s, q4s):
    outvals = []
    for t, chandra_eci, q1, q2, q3, q4 in zip(ephem_Times, chandra_ecis, q1s, q2s, q3s, q4s):
        alt = np.sqrt(np.sum(chandra_eci**2))/1e3
        date = re.sub(r'\.000$', '', DateTime(t).date)
        q_att = Quaternion.normalize([q1,q2,q3,q4])
        vis, illum, rays = taco.calc_earth_vis(chandra_eci, q_att, ngrid=opt.ngrid)
        title = '%s alt=%6.0f illum=%6.4f' % (date, alt, illum)
        outvals.append((t, illum, alt, q1, q2, q3, q4))
        if opt.verbose:
            print title, taco.norm(chandra_eci), q1, q2, q3, q4
        elif iproc == 0:
            print 'Progress: %d%%\r' % int((100. * len(outvals)) / len(ephem_Times) + 0.5),
            sys.stdout.flush()

def main():
    global opt
    opt, args = get_options()

    # Get orbital ephemeris in requested time range
    print 'Reading orbital ephemeris file', opt.ephemfile
    ephem = Ska.Table.read_table(opt.ephemfile)
    tstart = DateTime(opt.tstart).secs
    tstop = DateTime(opt.tstop).secs
    ephem = ephem[(ephem.Time >= tstart)
                              & (ephem.Time <= tstop)][::opt.sample]
    if len(ephem) == 0:
        print 'Error: ephemeris file has no overlap with specified time interval'
        sys.exit(0)
    chandra_ecis = np.array([ephem.X, ephem.Y, ephem.Z]).copy().transpose()
    ephem_times = ephem.Time.copy()


    # Get spacecraft attitude in requested time range at the same sampling as ephemeris
    dt = ephem.Time[1] - ephem.Time[0]
    print ('Fetching attitude telemetry between %s and %s at %d second intervals'
           % (DateTime(ephem.Time[0]).date, DateTime(ephem.Time[-1]).date, dt))
    cols, atts = fetch(start=ephem.Time[0], stop=ephem.Time[-1], dt=dt, time_format='secs',
                 colspecs=['aoattqt1', 'aoattqt2',  'aoattqt3',  'aoattqt4'])
    atts = np.rec.fromrecords(atts, names=cols)
    q1s = Ska.Numpy.interpolate(atts.aoattqt1, atts.date, ephem.Time)
    q2s = Ska.Numpy.interpolate(atts.aoattqt2, atts.date, ephem.Time)
    q3s = Ska.Numpy.interpolate(atts.aoattqt3, atts.date, ephem.Time)
    q4s = Ska.Numpy.interpolate(atts.aoattqt4, atts.date, ephem.Time)

    # Divy up calculations amongst the n-processors
    i0s = range(0, len(q1s), len(q1s) // opt.nproc + 1)
    i1s = i0s[1:] + [len(q1s)]

    # Calculate illumination in a separate process over each sub-interval
    for iproc, i0, i1 in zip(itertools.count(), i0s, i1s):
        calc_vis_values(iproc, ephem_times[i0:i1], chandra_ecis[i0:i1],
                        q1s[i0:i1], q2s[i0:i1], q3s[i0:i1], q4s[i0:i1])

if __name__ == '__main__':
    main()
    print
    
