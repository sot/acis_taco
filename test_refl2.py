import taco
import taco2
import numpy
import Chandra.acis_esa
from Quaternion import Quat

p_earth_body = np.array([-11913349.37481491,   1600513.79810546,   6787847.04879577])
orbit_xyz = -p_earth_body
orbit_xyz = np.array([0., 0., -100000e3])

illums = []
illums2 = []
pitchs = numpy.linspace(0, 90.0, 30)
esa_directs = []
esa_refls = []

for pitch in pitchs:
    print pitch
    att = [0, pitch, 0]
    vis, illum, rays = taco.calc_earth_vis(orbit_xyz, att, n_radiator_x=3, n_radiator_y=4, ngrid=100)
    vis2, illum2, rays2 = taco2.calc_earth_vis(orbit_xyz, att)
    illums.append(illum)
    illums2.append(illum2)
    direct, refl, total = Chandra.acis_esa.earth_solid_angle(Quat(att), orbit_xyz)
    esa_directs.append(direct)
    esa_refls.append(refl)


clf()
plot(pitchs, [x[0] for x in illums] , '-b', label='tla direct')
plot(pitchs, [x[0] for x in illums2] , '--b', label='tla direct2', linewidth=4)

plot(pitchs, esa_directs, '--k', label='nraw direct')
plot(pitchs, esa_refls, '--g', label='nraw refl')
plot(pitchs, [x+y for (x,y) in zip(esa_directs, esa_refls)], '--g', label='nraw total', linewidth=3)

legend()
grid()
xlabel('Earth Pitch (deg)')
title('Direct and reflected illum vs. Earth pitch')
