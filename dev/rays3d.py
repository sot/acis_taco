# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy
import taco
from enthought.mayavi import mlab

# source /data/cygob2/imgview/lopy/bin/activate.csh
# setenv LD_LIBRARY_PATH /usr/local/lib/vtk-5.6:/usr/local/Trolltech/Qt-current/lib

rad_corners =  numpy.array([[-0.209,  0.555, 0.036],
                            [ 0.209,  0.555, 0.036],
                            [-0.209, -0.389, 0.036],
                            [ 0.209, -0.389, 0.036]])

planes = taco.make_taco()
xs = []
ys = []
zs = []
triangles = []
for i, plane in enumerate(planes):
    xs.extend([plane.p0[0], plane.p1[0], plane.p2[0]])
    ys.extend([plane.p0[1], plane.p1[1], plane.p2[1]])
    zs.extend([plane.p0[2], plane.p1[2], plane.p2[2]])
    triangles.append((3*i, 3*i+1, 3*i+2))


p_earth_body = numpy.array([-20000e3,   0,   30000e3])
p_earth_body = numpy.array([-11913349.37481491,   1600513.79810546,   6787847.04879577])
att = [0, 0, 0]

vis, illum, outs = taco.calc_earth_vis(-p_earth_body, att, ngrid=30, to_earth=True, max_reflect=0)

mask_points = (len(outs) // 300) or None
rx = [out[0][0] for out in outs]
ry = [out[0][1] for out in outs]
rz = [out[0][2] for out in outs]
ru = [out[1][0] * (2 if out[3] else 1) for out in outs]
rv = [out[1][1] * (2 if out[3] else 1) for out in outs]
rw = [out[1][2] * (2 if out[3] else 1) for out in outs]

mlab.clf()
mlab.triangular_mesh(xs, ys, zs, triangles, color=(1, 1, 0), opacity=0.9)
mlab.triangular_mesh(xs, ys, zs, triangles, color=(1, .8, 0), representation='wireframe', line_width=4)
mlab.quiver3d(rx, ry, rz, ru, rv, rw, scale_factor=1.0, mask_points=mask_points, vmin=1, vmax=2)
mlab.show()
