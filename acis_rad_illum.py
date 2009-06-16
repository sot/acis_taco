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

import Ska.Table
import Quaternion
from Chandra.Time import DateTime
from Ska.TelemArchive.fetch import fetch
import Ska.Numpy
import matplotlib.pyplot as plt
from Ska.Matplotlib import plot_cxctime

import taco

def get_options():
    from optparse import OptionParser
    parser = OptionParser()
    parser.set_defaults()
    parser.add_option("--tstart",
                      default='2009:060:00:00:00',
                      help="Start time")
    parser.add_option("--tstop",
                      default='2009:067:00:00:00',
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
                      default=4,
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

def plot_vis_image(title, date, rays, vis, iproc):
    plt.figure(iproc, figsize=(5,5))
    plt.clf()
    blocked = vis < 1e-5
    if len(blocked.nonzero()) > 0:
        pos = rays[blocked]
        plt.plot(pos[:,1], pos[:,2], '.k')

    for alpha, ray in zip(vis[~blocked], rays[~blocked]):
        plt.plot([ray[1]], [ray[2]], '.r', alpha=alpha)

    plt.xlim(-0.5,0.5)
    plt.ylim(-0.5,0.5)
    plt.title(title)
    filename = os.path.join(opt.out, date + '.png')
    if opt.verbose:
        print 'Writing image', filename
    plt.savefig(filename)

def calc_vis_values(queue, iproc, ephem_Times, chandra_ecis, q1s, q2s, q3s, q4s):
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
        if opt.movie:
            plot_vis_image(title, date, rays, vis, iproc)
    queue.put(outvals)

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

    if opt.movie:
        if len(atts) > 250:
            print "Error: movie option will produce more than 250 images.  Change code if needed."
            sys.exit(0)
        if not os.path.exists(opt.out):
            os.makedirs(opt.out)

    # Divy up calculations amongst the n-processors
    i0s = range(0, len(q1s), len(q1s) // opt.nproc + 1)
    i1s = i0s[1:] + [len(q1s)]

    # Calculate illumination in a separate process over each sub-interval
    queues = []
    procs = []
    for iproc, i0, i1 in zip(itertools.count(), i0s, i1s):
        queue = Queue()
        proc = Process(target=calc_vis_values, args=(queue, iproc, ephem_times[i0:i1], chandra_ecis[i0:i1],
                                                  q1s[i0:i1], q2s[i0:i1], q3s[i0:i1], q4s[i0:i1]))
        proc.start()
        procs.append(proc)
        queues.append(queue)

    # Join the results from each processor at the end
    outvals = []
    for proc, queue in zip(procs, queues):
        outvals.extend(queue.get())
        proc.join()

    print

    # Plot illumination versus date
    fig = plt.figure(1, figsize=(6,4))
    illum = np.rec.fromrecords(outvals, names=['time', 'illum', 'alt', 'q1', 'q2', 'q3', 'q4'])
    ticklocs, fig, ax = plot_cxctime(illum.time, illum.illum, fmt='-b')
    ax.set_title('ACIS radiator illumination')
    ax.set_ylabel('Illumination (steradians)')
    filename = opt.out + '.png'
    fig.savefig(filename)
    print 'Create image file', filename

    # Write results to FITS table
    filename = opt.out + '.fits'
    Ska.Table.write_fits_table(opt.out + '.fits', illum)
    print 'Created FITS table', filename

    if opt.movie:
        print 'To make a movie run the following command:'
        print 'convert -delay 30 %s/*.png -loop 0 %s.gif' % (opt.out, opt.out)

if __name__ == '__main__':
    main()
    print
    
