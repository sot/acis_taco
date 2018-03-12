# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy
import taco
from enthought.mayavi import mlab
from itertools import repeat

# source /data/cygob2/imgview/lopy/bin/activate.csh
# setenv LD_LIBRARY_PATH /usr/local/lib/vtk-5.6:/usr/local/Trolltech/Qt-current/lib

rad_corners =  numpy.array([[-0.209,  0.555, 0.036],
                            [ 0.209,  0.555, 0.036],
                            [-0.209, -0.389, 0.036],
                            [ 0.209, -0.389, 0.036]])

def make_taco(refl=0):
    """Define geometry for ACIS radiator sun-shade (aka the space taco)."""
    x = [-0.25 * (refl+1),  0.25 * (refl+1)]
    y = [-0.689,  -0.333,  -0.0635,  0.553, 0.6896]
    z = [0.0, 0.777, 0.777, 0.421, 0.0]
    p = numpy.array([zip(repeat(x[0]), y, z),
                  zip(repeat(x[1]), y, z)])
    planes = [taco.Plane(p[0,0], p[0,1], p[0,2]),
              taco.Plane(p[0,0], p[0,2], p[0,3]),
              taco.Plane(p[0,0], p[0,3], p[0,4]),
              taco.Plane(p[1,0], p[1,1], p[1,2]),
              taco.Plane(p[1,0], p[1,2], p[1,3]),
              taco.Plane(p[1,0], p[1,3], p[1,4]),
              taco.Plane(p[0,0], p[0,4], p[1,4]),
              taco.Plane(p[0,0], p[1,0], p[1,4]),
              ]
    return planes

xs = []
ys = []
zs = []
triangles = []

for i, plane in enumerate(make_taco(0)):
    xs.extend([plane.p0[0], plane.p1[0], plane.p2[0]])
    ys.extend([plane.p0[1], plane.p1[1], plane.p2[1]])
    zs.extend([plane.p0[2], plane.p1[2], plane.p2[2]])
    triangles.append((3*i, 3*i+1, 3*i+2))

pts = [[-0.209,  0.555, 0.036],
       [ 0.209,  0.555, 0.036],
       [-0.209, -0.389, 0.036],
       [ 0.209, -0.389, 0.036]]

rad_xs = [pts[0][0], pts[1][0], pts[3][0],
          pts[0][0], pts[2][0], pts[3][0]]
rad_ys = [pts[0][1], pts[1][1], pts[3][1],
          pts[0][1], pts[2][1], pts[3][1]]
rad_zs = [pts[0][2], pts[1][2], pts[3][2],
          pts[0][2], pts[2][2], pts[3][2]]
rad_triangles = [(0,1,2), (3,4,5)]

p_earth_body = numpy.array([-20000e3,   0,   30000e3])
p_earth_body = numpy.array([-11913349.37481491,   1600513.79810546,   6787847.04879577])
pos = p_earth_body / numpy.sqrt(numpy.sum(p_earth_body**2)) * 2
att = [0, 0, 0]

mlab.clf()
mlab.triangular_mesh(rad_xs, rad_ys, rad_zs, rad_triangles, color=(0, 0, 1), opacity=1)
mlab.triangular_mesh(xs, ys, zs, triangles, color=(1, 1, 0), opacity=0.9)
mlab.triangular_mesh(xs, ys, zs, triangles, color=(1, .8, 0), representation='wireframe', line_width=4)
mlab.points3d(pos[0], pos[1], pos[2], scale_factor=1.0)

if 'show_rays' in globals() and show_rays: 
    vis, illum, outs = taco.calc_earth_vis(-p_earth_body, att, ngrid=30, to_earth=True, max_reflect=0)

    mask_points = (len(outs) // 300) or None
    rx = [out[0][0] for out in outs]
    ry = [out[0][1] for out in outs]
    rz = [out[0][2] for out in outs]
    ru = [out[1][0] * (2 if out[3] else 1) for out in outs]
    rv = [out[1][1] * (2 if out[3] else 1) for out in outs]
    rw = [out[1][2] * (2 if out[3] else 1) for out in outs]

    mlab.quiver3d(rx, ry, rz, ru, rv, rw, scale_factor=1.0, mask_points=mask_points, vmin=1, vmax=2)


mlab.show()
