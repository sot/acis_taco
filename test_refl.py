import taco
import numpy
import Chandra.acis_esa
from Quaternion import Quat

p_earth_body = np.array([-11913349.37481491,   1600513.79810546,   6787847.04879577])
orbit_xyz = -p_earth_body
orbit_xyz = np.array([0., 0., -100000e3])
orbit_xyz = np.array([0., 0., -10000e3])

illums = []
pitchs = numpy.linspace(0, 90.0, 30)
esa_directs = []
esa_refls = []

for pitch in pitchs:
    print pitch
    att = [0, 0, pitch]
    vis, illum, rays = taco.calc_earth_vis(orbit_xyz, att, ngrid=100,
                                           n_radiator_x=5, n_radiator_y=8,
                                           reflect_atten=0.9,
                                           max_reflect=3)
    illums.append(illum)
    direct, refl, total = Chandra.acis_esa.earth_solid_angle(Quat(att), orbit_xyz)
    esa_directs.append(direct)
    esa_refls.append(refl)


clf()
plot(pitchs, [x[0] for x in illums] , '-b', label='tla direct')
plot(pitchs, [x[1] for x in illums] , '-r', label='tla reflect 1')
plot(pitchs, [x[2] for x in illums] , '-m', label='tla reflect 2')
plot(pitchs, [x[3] for x in illums] , '-c', label='tla reflect 3')
plot(pitchs, [sum(x) for x in illums] , '-b', label='tla total', linewidth=3)

plot(pitchs, esa_directs, '--k', label='nraw direct')
plot(pitchs, esa_refls, '--g', label='nraw refl')
plot(pitchs, [x+y for (x,y) in zip(esa_directs, esa_refls)], '--g', label='nraw total', linewidth=3)

legend()
grid()
xlabel('Earth Pitch (deg)')
title('Direct and reflected illum vs. Earth pitch')
