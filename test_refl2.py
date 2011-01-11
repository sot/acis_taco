import taco
import taco2
import numpy
import Chandra.acis_esa
from Quaternion import Quat

p_earth_body = numpy.array([-11913349.37481491,   1600513.79810546,   6787847.04879577])
orbit_xyz = -p_earth_body
orbit_xyz = numpy.array([0., 0., -1000000e3])

illums2 = []
pitchs = numpy.linspace(0, 90.0, 200)
esa_directs = []
esa_refls = []

for pitch in pitchs:
    print pitch
    att = [0, pitch, 0]
    vis2, illum2, rays2 = taco2.calc_earth_vis(orbit_xyz, att, max_reflect=4)
    illums2.append(illum2)
print

clf()
plot(pitchs, [x[0] for x in illums2] , '-b', label='tla direct-2', linewidth=2)
plot(pitchs, [x[1] for x in illums2] , '-r', label='tla refl1-2', linewidth=2)
plot(pitchs, [x[2] for x in illums2] , '-g', label='tla refl2-2', linewidth=2)
plot(pitchs, [x[3] for x in illums2] , '-c', label='tla refl2-2', linewidth=2)
plot(pitchs, [x[4] for x in illums2] , '-m', label='tla refl2-2', linewidth=2)

legend()
grid()
xlabel('Earth Pitch (deg)')
title('Direct and reflected illum vs. Earth pitch')
