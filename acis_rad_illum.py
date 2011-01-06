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
import Ska.engarchive.fetch_sci as fetch
import Ska.Numpy
import matplotlib
if __name__ == '__main__':
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Ska.Matplotlib import plot_cxctime

import taco

import Chandra.acis_esa

def get_options():
    from optparse import OptionParser
    parser = OptionParser()
    parser.set_defaults()
    parser.add_option("--tstart",
                      default='2010:232:19:00:00',
                      help="Start time")
    parser.add_option("--tstop",
                      default='2010:237:22:00:00',
                      help="Stop time")
    parser.add_option("--out",
                      default='out',
                      help="Output root name")
    parser.add_option("--sample",
                      type='int',
                      default=600,
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

def calc_vis_values(queue, iproc, times, chandra_ecis, q1s, q2s, q3s, q4s):
    outvals = []
    for t, chandra_eci, q1, q2, q3, q4 in zip(times, chandra_ecis, q1s, q2s, q3s, q4s):
        alt = np.sqrt(np.sum(chandra_eci**2))/1e3
        date = re.sub(r'\.000$', '', DateTime(t).date)
        q_att = Quaternion.normalize([q1,q2,q3,q4])
        vis, illum, rays = taco.calc_earth_vis(chandra_eci, q_att, ngrid=opt.ngrid)
        title = '%s alt=%6.0f illum=%6.4f' % (date, alt, illum)
        outvals.append((t, illum['direct'], illum['reflect1'], illum['reflect2'], alt, q1, q2, q3, q4))
        if opt.verbose:
            print title, taco.norm(chandra_eci), q1, q2, q3, q4
        elif iproc == 0:
            print 'Progress: %d%%\r' % int((100. * len(outvals)) / len(times) + 0.5),
            sys.stdout.flush()
        if opt.movie:
            plot_vis_image(title, date, rays, vis, iproc)
    queue.put(outvals)

def main():
    global opt
    opt, args = get_options()
    tstart = DateTime(opt.tstart).secs
    tstop = DateTime(opt.tstop).secs

    # Get orbital ephemeris in requested time range
    print ('Fetching ephemeris')
    ephem_x = fetch.MSID('orbitephem1_x', opt.tstart, opt.tstop)
    ephem_y = fetch.MSID('orbitephem1_y', opt.tstart, opt.tstop)
    ephem_z = fetch.MSID('orbitephem1_z', opt.tstart, opt.tstop)
    ephem_times = ephem_x.times.copy()

    # Get spacecraft attitude in requested time range at the same sampling as ephemeris
    print ('Fetching attitude telemetry between {0} and {1}'.format(
            opt.tstart, opt.tstop))
    qatts = fetch.MSIDset(['aoattqt1', 'aoattqt2',  'aoattqt3',  'aoattqt4'],
                          opt.tstart, opt.tstop)
    #cols, atts = fetch(start=ephem.Time[0], stop=ephem.Time[-1], dt=dt, time_format='secs',
    #             colspecs=)
    # atts = np.rec.fromrecords(atts, names=cols)
    q1s = qatts['aoattqt1'].vals[::opt.sample]
    q2s = qatts['aoattqt2'].vals[::opt.sample]
    q3s = qatts['aoattqt3'].vals[::opt.sample]
    q4s = qatts['aoattqt4'].vals[::opt.sample]
    q_times = qatts['aoattqt1'].times[::opt.sample]

    ephem_x_vals = Ska.Numpy.interpolate(ephem_x.vals, ephem_times, q_times)
    ephem_y_vals = Ska.Numpy.interpolate(ephem_y.vals, ephem_times, q_times)
    ephem_z_vals = Ska.Numpy.interpolate(ephem_z.vals, ephem_times, q_times)

    chandra_ecis = np.array([ephem_x_vals, ephem_y_vals, ephem_z_vals]).copy().transpose()

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
        proc = Process(target=calc_vis_values, args=(queue, iproc, q_times[i0:i1], chandra_ecis[i0:i1],
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

    esa_directs = []
    esa_refls = []
    for t, q1, q2, q3, q4, x, y, z in zip(
        q_times, q1s, q2s, q3s, q4s, ephem_x_vals, ephem_y_vals, ephem_z_vals):
        direct, refl, total = Chandra.acis_esa.earth_solid_angle(
            Quaternion.Quat([q1, q2, q3, q4]), np.array([x, y, z]))

        esa_directs.append(direct)
        esa_refls.append(refl)

    # Plot illumination versus date
    fig = plt.figure(1, figsize=(6,4))
    plt.clf()
    illum = np.rec.fromrecords(outvals, names=['time', 'direct', 'reflect1', 'reflect2', 'alt', 'q1', 'q2', 'q3', 'q4'])
    ticklocs, fig, ax = plot_cxctime(illum.time, illum.direct, '-b')
    plot_cxctime(illum.time, illum.reflect1, '-r')
    plot_cxctime(q_times, esa_directs, '-c')
    plot_cxctime(q_times, esa_refls, '-m')
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
    
