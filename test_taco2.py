#!/usr/bin/env /proj/sot/ska/bin/python
import sys
sys.path.append('/proj/sot/ska/share/taco')

import taco2
import numpy
from Quaternion import Quat
import matplotlib.pyplot as plt

plt.figure(1, figsize=(6,8))
plt.clf()
pitchs = numpy.linspace(0, 90.0, 46)

# Loop over position of Chandra relative to Earth
iplot = 1
for orbit_xyz in numpy.array([[11913349, -1600513,  -6787847], # Test case 1
                              [0., 0., -1000000e3], # Earth along +Z and far away
                              [0., 0., -7000e3]]):   # Earth along +Z and near
    print "Chandra position =", orbit_xyz
    print "RA, Dec, Roll, total, direct, reflected"
    directs = []
    reflects = []
    totals = []
    for pitch in pitchs:
        att = [0, pitch, 0]   # RA, Dec, Roll.  "Dec" set to pitch
        vis, illum, rays = taco2.calc_earth_vis(orbit_xyz, att, max_reflect=10)
        direct = illum[0]
        reflect = sum(illum[1:])
        total = sum(illum)
        directs.append(direct)   
        reflects.append(reflect)
        totals.append(total)
        print "%5.1f %5.1f %5.1f %8.4e %8.4e %8.4e" % (
            att[0], att[1], att[2], total, direct, reflect)

    plt.subplot(3, 1, iplot)
    plt.plot(pitchs, totals , '-k', label='Total' )
    plt.plot(pitchs, directs, '-b', label='Direct' )
    plt.plot(pitchs, reflects, '-r', label='Reflected'  )
    plt.grid()
    if iplot == 2:
        plt.legend()
    if iplot == 3:
        plt.xlabel('Earth Pitch (deg)')
    if iplot == 1:
        plt.title('Direct and reflected illum vs. Earth pitch')
    plt.savefig('test_taco2_pitch.png')
    iplot += 1

print 'Wrote plot file test_taco2_pitch.png'

