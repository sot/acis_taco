"""
Define geometry for ACIS radiator and sunshade and perform raytrace
calculation of Earth illumination on the radiator surface.
"""
import os

from itertools import repeat
from Quaternion import Quat
import numpy

def make_taco():
    """Define geometry for ACIS radiator sun-shade (aka the space taco) in mm."""
    y_pnts = numpy.array([-689.0,  -333,  -63,  553, 689]) + TACO_Y_OFF
    z_pnts = numpy.array([0.0, 777, 777, 421, 0]) - RAD_Z_OFF

    y_edge = numpy.arange(max(y_pnts))
    z_edge = interpolate(z_pnts, y_pnts, y_edge)

    return z_edge

def make_radiator(n_radiator_x=3, n_radiator_y=4):
    """Specify points on the ACIS radiator surface.
    Corners of radiator at: (lengths in meters)
    [[-0.209,  0.555, 0.036],
     [ 0.209,  0.555, 0.036],
     [-0.209, -0.389, 0.036],
     [ 0.209, -0.389, 0.036]]
     Center-y is 0.083, length is 0.944.
     Z standoff height is 0.36 
    """
    rad_x = numpy.linspace(-200, 200, n_radiator_x)
    rad_y = numpy.linspace(-389, 555, n_radiator_y) + TACO_Y_OFF

    #rad_x = numpy.array([0.0])
    #rad_y = numpy.array([0.0]) + TACO_Y_OFF

    return rad_x, rad_y

def calc_earth_vis(p_chandra_eci,
                   chandra_att,
                   max_reflect=0,
                   calc_vis=True):
    """Calculate the relative Earth visibility for the ACIS radiator given
    the Chandra orbit position ``p_chandra_eci`` and attitude ``chandra_att``.

    The relative visibility is normalized so that 1.0 represents the entire
    radiator having visibility toward the surface point and being exactly
    normal toward the surface.

    Total illumination gives the effective solid angle (steradians) of the
    visible Earth from the radiator perspective.  It is averaged over the
    different points on the radiator.

    :param p_chandra_eci: Chandra orbital position [x, y, z] (meters)
    :param chandra_att: Chandra attitude [ra, dec, roll] (deg)
    
    :returns: relative visibility, total illumination, projected rays
    """

    # Calculate position of earth in ECI and Chandra body coords.  
    q_att = Quat(chandra_att) 
    
    # For T = attitude transformation matrix then p_body = T^-1 p_eci
    p_earth_body = numpy.dot(q_att.transform.transpose(), -numpy.array(p_chandra_eci))

    # Quaternion to transform x-axis to the body-earth vector
    q_earth = quat_x_v2(p_earth_body)
    open_angle = numpy.arcsin(RAD_EARTH / numpy.sqrt(numpy.sum(p_earth_body**2)))
    rays, earth_solid_angle = sphere_grid(N_GRID, open_angle)
    n_rays = len(rays)
    # or          = numpy.dot(rays, q_earth.transform.transpose())
    rays_to_earth = numpy.dot(q_earth.transform, rays.transpose()).transpose()  # shape (n_rays, 3)

    # Accept only rays with a positive Z component and make sure no X component is < 1e-6
    rays_to_earth = rays_to_earth[rays_to_earth[:, 2] > 0.0, :]
    rays_x = numpy.abs(rays_to_earth[:, 0])
    rays_x[rays_x < 1e-6] = 1e-6
    rays_x_all = rays_x

    # Initialize outputs
    vis = numpy.zeros((max_reflect+1, n_rays))
    illum = numpy.zeros(max_reflect+1)

    # Main ray-trace loop.  Calculate ray visibility for an increasing number
    # of reflections.  Rays that get blocked after N reflections are candidates
    # for getting out after N+1 reflections.
    i_rays = numpy.arange(n_rays)
    out_rays = []

    for rad_x in RAD_XS:
        for rad_y in RAD_YS:
            rays_x = rays_x_all
            rays_y = rays_to_earth[:, 1]
            rays_z = rays_to_earth[:, 2]
            for refl in range(max_reflect + 1):
                taco_x = TACO_X_OFF * (refl + 1)
                dx = (taco_x - rad_x) / rays_x

                # Find rays that intersect shade within Y limits of the shade (0, N_TACO mm)
                ray_taco_iys = numpy.array(rad_y + rays_y * dx, dtype=numpy.int)
                y_ok = (ray_taco_iys >= 0) & (ray_taco_iys < N_TACO)
                y_not_ok = ~y_ok
                ray_taco_iys[y_not_ok] = 0

                # From those rays find ones below the TACO_Z_EDGES curve
                ray_taco_zs = rays_z * dx
                z_ok = ray_taco_zs > TACO_Z_EDGES[ray_taco_iys]
                y_z_ok = (y_ok & z_ok) | y_not_ok

                # Calculate the delta-visibility for good rays (cos(incidence_angle) = Z component)
                d_vis = rays_z[y_z_ok] * REFLECT_ATTEN**refl / N_RAD_POINTS
                illum[refl] += earth_solid_angle * numpy.sum(d_vis) / n_rays

                blocked = ~y_z_ok
                rays_x = rays_x[blocked]
                rays_y = rays_y[blocked]
                rays_z = rays_z[blocked]

    return vis, illum, out_rays

def interpolate(yin, xin, xout):
    """
    Interpolate the curve defined by (xin, yin) at points xout.  The array
    xin must be monotonically increasing.  The output has the same data type as
    the input yin.

    :param yin: y values of input curve
    :param xin: x values of input curve
    :param xout: x values of output interpolated curve

    @:rtype: numpy array with interpolated curve
    """
    lenxin = len(xin)

    i1 = numpy.searchsorted(xin, xout)
    i1[ i1==0 ] = 1
    i1[ i1==lenxin ] = lenxin-1

    x0 = xin[i1-1]
    x1 = xin[i1]
    y0 = yin[i1-1]
    y1 = yin[i1]

    return (xout - x0) / (x1 - x0) * (y1 - y0) + y0

def norm(vec):
    return vec / numpy.sqrt(numpy.sum(vec**2))

def quat_x_v2(v2):
    """Generate quaternion that rotates X-axis into v2 by the shortest path"""
    x = numpy.array([1.,0,0])
    v2 = norm(numpy.array(v2))
    dot = numpy.dot(x, v2)
    if abs(dot) > 1-1e-8:
        x = norm(numpy.array([1., 0., 1e-5]))
        dot = numpy.dot(v2, x)
    angle = numpy.arccos(dot)
    axis = norm(numpy.cross(x, v2))
    sin_a = numpy.sin(angle / 2)
    cos_a = numpy.cos(angle / 2)
    return Quat([axis[0] * sin_a,
                            axis[1] * sin_a,
                            axis[2] * sin_a,
                            cos_a])

def sphere_grid(ngrid, open_angle):
    """Calculate approximately uniform spherical grid of rays containing
    ``ngrid`` points and extending over the opening angle ``open_angle``
    (radians).

    :returns: numpy array of unit length rays, grid area (steradians)
    """
    from math import sin, cos, radians, pi, sqrt
    
    grid_area = 2*pi*(1-cos(open_angle))
    gridsize = sqrt(grid_area / ngrid)

    grid = []
    theta0 = pi/2-open_angle
    n_d = int(round(open_angle / gridsize))
    d_d = open_angle / n_d

    for i_d in range(0, n_d+1):
        dec = i_d * d_d + pi/2 - open_angle
        if abs(i_d) != n_d:
            n_r = int(round( 2*pi * cos(dec) / d_d))
            d_r = 2*pi / n_r
        else:
            n_r = 1
            d_r = 1
        for i_r in range(0, n_r):
            ra = i_r * d_r
            # This has x <=> z (switched) from normal formulation to make the
            # grid centered about the x-axis
            grid.append((sin(dec), sin(ra) * cos(dec), cos(ra) * cos(dec)))

    # (x, y, z) = zip(*grid)  (python magic)
    return numpy.array(grid), grid_area

class Line(object):
    """Line from p0 to p1"""
    def __init__(self, p0, p1):
        self.p0 = numpy.array(p0)
        self.p1 = numpy.array(p1)
        self.u = self.p1 - self.p0
        self.len = numpy.sqrt(numpy.sum(self.u**2))  # line length
        self.u /= self.len              # unit vector in line direction

    def __str__(self):
        return "%s %s" % (str(self.p0), str(self.p1))

    def __repr__(self):
        return str(self)

def py_plane_line_intersect(p, l):
    """Determine if the line ``l`` intersects the plane ``p``.  Pure python.

    :rtype: boolean
    """
    mat = numpy.array([l.p0 - l.p1,
                    p.p1 - p.p0,
                    p.p2 - p.p0]).transpose()
    try:
        t, u, v = numpy.dot(numpy.linalg.inv(mat), l.p0 - p.p0)
        intersect = 0 <= u <= 1 and 0 <= v <= 1 and u + v <= 1 and t > 0
    except numpy.linalg.LinAlgError, msg:
        # This presumably occurred because is parallel to plane
        intersect = False

    return intersect

RAD_EARTH = 6371e3
RAD_Z_OFF = 36
TACO_X_OFF = 250
TACO_Y_OFF = 689
N_GRID = 100
REFLECT_ATTEN=0.9
TACO_Z_EDGES = make_taco()
N_TACO = len(TACO_Z_EDGES)
RAD_XS, RAD_YS = make_radiator(3, 4)
N_RAD_POINTS = len(RAD_XS) * len(RAD_YS)

