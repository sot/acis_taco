# Licensed under a 3-clause BSD style license - see LICENSE.rst
import taco2_r14
import taco2
import taco
import numpy
import Chandra.acis_esa
from Quaternion import Quat

p_earth_body = numpy.array([-11913349.37481491,   1600513.79810546,   6787847.04879577])
orbit_xyz = -p_earth_body
orbit_xyz = numpy.array([0., 0., -1000000e3])
orbit_xyz = numpy.array([0., 0., -7000e3])

illums1 = []
illums2 = []
pitchs = numpy.linspace(0, 90.0, 30)
esa_directs = []
esa_refls = []

for pitch in pitchs:
    att = [0, pitch, 0]
    vis1, illum1, rays1 = taco2_r14.calc_earth_vis(orbit_xyz, att, max_reflect=5)
    #vis1, illum1, rays1 = taco.calc_earth_vis(orbit_xyz, att, max_reflect=5)
    vis2, illum2, rays2 = taco2.calc_earth_vis(orbit_xyz, att, max_reflect=5)
    print pitch, illum1, illum2
    illums1.append(illum1)
    illums2.append(illum2)
print

if 1:
    figure(1)
    clf()
    plot(pitchs, [sum(x) for x in illums2] , '-k', label='tla direct-2' )
    plot(pitchs, [x[0] for x in illums2] , '-b', label='tla direct-2' )
    plot(pitchs, [x[1] for x in illums2] , '-r', label='tla refl1-2'  )
    plot(pitchs, [x[2] for x in illums2] , '-g', label='tla refl2-2'  )
    plot(pitchs, [x[3] for x in illums2] , '-c', label='tla refl2-2'  )
    plot(pitchs, [x[4] for x in illums2] , '-m', label='tla refl2-2'  )

    plot(pitchs, [sum(x) for x in illums1] , '--k', label='tla direct-2', linewidth=3)
    plot(pitchs, [x[0] for x in illums1] , '--b', label='tla direct-2', linewidth=3)
    plot(pitchs, [x[1] for x in illums1] , '--r', label='tla refl1-2', linewidth=3)
    plot(pitchs, [x[2] for x in illums1] , '--g', label='tla refl2-2', linewidth=3)
    plot(pitchs, [x[3] for x in illums1] , '--c', label='tla refl2-2', linewidth=3)
    plot(pitchs, [x[4] for x in illums1] , '--m', label='tla refl2-2', linewidth=3)

    # legend()
    grid()
    xlabel('Earth Pitch (deg)')
    title('Direct and reflected illum vs. Earth pitch')
