#!/usr/bin/env python
# Create animation with:
# convert -delay 30 out*.png -loop 0 rad_illum_2009124.gif
#
import sys
import taco
import Quaternion
import numpy as np
import Ska.Table
from Chandra.Time import DateTime
from Ska.TelemArchive.fetch import fetch
import Ska.Numpy
import matplotlib.pyplot as plt
from multiprocessing import Process, Queue
from Ska.Matplotlib import plot_cxctime

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
    parser.add_option("--sample",
                      type='int',
                      default=1,
                      help="Sample")
    parser.add_option("--ephemfile",
                      default='orbitf352123502N001_eph0.fits',
                      help="Orbit ephemeris file")
    parser.add_option("--out",
                      default='out',
                      help="Output root name")
    parser.add_option("--size",
                      type='float',
                      default=4.0,
                      help="Output image size (inches)")
    parser.add_option("--nproc",
                      type='int',
                      default=4,
                      help="Number of processors to use")
    (opt, args) = parser.parse_args()
    return (opt, args)

opt, args = get_options()

# Get orbital ephemeris in requested time range
if 'chandra_ecis' not in globals():
    ephem = Ska.Table.read_table(opt.ephemfile)
    tstart = DateTime(opt.tstart).secs
    tstop = DateTime(opt.tstop).secs
    ephem = ephem[(ephem.Time >= tstart)
                              & (ephem.Time <= tstop)][::opt.sample]
    chandra_ecis = np.array([ephem.X, ephem.Y, ephem.Z]).transpose()
    ephem_times = ephem.Time.copy()

# Get spacecraft attitude in requested time range at the same sampling as ephemeris
if 'q1s' not in globals():
    dt = ephem.Time[1] - ephem.Time[0]
    cols, atts = fetch(start=ephem.Time[0], stop=ephem.Time[-1], dt=dt, time_format='secs',
                 colspecs=['aoattqt1', 'aoattqt2',  'aoattqt3',  'aoattqt4'])
    atts = np.rec.fromrecords(atts, names=cols)
    q1s = Ska.Numpy.interpolate(atts.aoattqt1, atts.date, ephem.Time)
    q2s = Ska.Numpy.interpolate(atts.aoattqt2, atts.date, ephem.Time)
    q3s = Ska.Numpy.interpolate(atts.aoattqt3, atts.date, ephem.Time)
    q4s = Ska.Numpy.interpolate(atts.aoattqt4, atts.date, ephem.Time)

def plot_vis_image():
    if t > DateTime('2009:124:03:08:55.000').secs and t < DateTime('2009:124:05:08:57').secs :
        plt.figure(i, figsize=(5,5))
        i += 1
        plt.clf()
        blocked = vis < 0.0001
        if len(blocked.nonzero()) > 0:
            pos = rays[blocked]
            plt.plot(pos[:,1], pos[:,2], '.k')

        for alpha, ray in zip(vis[~blocked], rays[~blocked]):
            plt.plot([ray[1]], [ray[2]], '.r', alpha=alpha)

        plt.xlim(-0.5,0.5)
        plt.ylim(-0.5,0.5)
        plt.title(title)
        plt.savefig('%s%02d.png' % (opt.out,i))

def calc_vis_values(queue, ephem_Times, chandra_ecis, q1s, q2s, q3s, q4s):
    outvals = []
    for t, chandra_eci, q1, q2, q3, q4 in zip(ephem_Times, chandra_ecis, q1s, q2s, q3s, q4s):
        alt = np.sqrt(np.sum(chandra_eci**2))/1e3
        q_att = Quaternion.normalize([q1,q2,q3,q4])
        vis, illum, rays = taco.calc_earth_vis(chandra_eci, q_att, ngrid=100)
        title = '%s %7.0f %7.4f' % (DateTime(t).date, alt, illum)
        outvals.append((t, illum, alt))
        print title, taco.norm(chandra_eci), q1, q2, q3, q4
    queue.put(outvals)

i0s = range(0, len(q1s), len(q1s) // opt.nproc + 1)
i1s = i0s[1:] + [len(q1s)]

del ephem
queues = []
procs = []
for i0, i1 in zip(i0s, i1s):
    print i0, i1
    queue = Queue()
    proc = Process(target=calc_vis_values, args=(queue, ephem_times[i0:i1], chandra_ecis[i0:i1],
                                              q1s[i0:i1], q2s[i0:i1], q3s[i0:i1], q4s[i0:i1]))
    proc.start()
    procs.append(proc)
    queues.append(queue)

outvals = []
for proc, queue in zip(procs, queues):
    outvals.extend(queue.get())
    proc.join()
    
fig = plt.figure(1, figsize=(6,4))
illum = np.rec.fromrecords(outvals, names=['time', 'illum', 'alt'])
ticklocs, fig, ax = plot_cxctime(illum.time, illum.illum, fmt='-b')
fig.savefig(opt.out + '.png')
ax.set_title('ACIS radiator illumination')
ax.set_ylabel('Illumination (steradians)')

Ska.Table.write_fits_table(opt.out + '.fits', illum)
