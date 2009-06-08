#!/usr/bin/env python
import numpy as np
from Quaternion import Quat
from itertools import repeat
import matplotlib.pyplot as plt

Rad_Earth = 6371e3

def make_taco():
    x = [-0.2,  0.2]
    y = [-0.5,  -0.3,  0.0,  0.3, 0.5]
    z = [0.0, 0.6, 0.5, 0.3, 0.0]
    p = np.array([zip(repeat(x[0]), y, z),
                  zip(repeat(x[1]), y, z)])
    planes = [Plane(p[0,0], p[0,1], p[0,2]),
              Plane(p[0,0], p[0,2], p[0,3]),
              Plane(p[0,0], p[0,3], p[0,4]),
              Plane(p[1,0], p[1,1], p[1,2]),
              Plane(p[1,0], p[1,2], p[1,3]),
              Plane(p[1,0], p[1,3], p[1,4]),
              Plane(p[0,0], p[0,4], p[1,4]),
              Plane(p[0,0], p[1,0], p[1,4]),
              ]
    return planes

def make_radiator():
    return np.array([[0.0, 0.0, 0.1]])

class Line(object):
    """Line from p0 to p1"""
    def __init__(self, p0, p1):
        self.p0 = np.array(p0)
        self.p1 = np.array(p1)
        self.u = self.p1 - self.p0
        self.len = np.sqrt(np.sum(self.u**2))  # line length
        self.u /= self.len              # unit vector in line direction

    def __str__(self):
        return "%s %s" % (str(self.p0), str(self.p1))

    def __repr__(self):
        return str(self)

class Plane(object):
    """Plane containing p0, p1, p2"""
    def __init__(self, p0, p1, p2):
        self.p0 = np.array(p0)
        self.p1 = np.array(p1)
        self.p2 = np.array(p2)
    
    def __str__(self):
        return "%s %s %s" % (str(self.p0), str(self.p1), str(self.p2))

    def __repr__(self):
        return str(self)

def plane_line_intersect(p, l):
    mat = np.array([(l.p0 - l.p1),
                    p.p1 - p.p0,
                    p.p2 - p.p0]).transpose()
    try:
        t, u, v = np.dot(np.linalg.inv(mat), l.p0 - p.p0)
        intersect = 0 <= u <= 1 and 0 <= v <= 1 and u + v <= 1 and t > 0
    except np.linalg.LinAlgError, msg:
        # This presumably occurred because is parallel to plane
        intersect = False

    return intersect

def sphere_grid(ngrid, open_angle):
    """Calculate approximately uniform spherical grid with spacing ``gridsize``"""
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
    return np.array(grid), grid_area

def norm(vec):
    return vec / np.sqrt(np.sum(vec**2))

def quat_x_v2(v2):
    """Generate quaternion that rotates X-axis into v2 by the shortest path"""
    x = np.array([1.,0,0])
    v2 = norm(np.array(v2))
    dot = np.dot(x, v2)
    if abs(dot) > 1-1e-8:
        x = norm(np.array([1., 0., 1e-5]))
        dot = np.dot(v2, x)
    angle = np.arccos(dot)
    axis = norm(np.cross(x, v2))
    sin_a = np.sin(angle / 2)
    cos_a = np.cos(angle / 2)
    return Quat([axis[0] * sin_a,
                            axis[1] * sin_a,
                            axis[2] * sin_a,
                            cos_a])

def calc_earth_vis(p_chandra_eci,
                   chandra_att,
                   ngrid=1000,
                   planes=make_taco(),
                   p_radiators=make_radiator(),):
    """Calculate the relative Earth visibility for the ACIS radiator given
    the Chandra orbit position ``p_chandra_eci`` and attitude ``chandra_att``.

    The relative visibility is normalized so that 1.0 represents the entire
    radiator seeing the full Earth at 100000 km.

    :param p_chandra_eci: Chandra orbital position [x, y, z] (meters)
    :param chandra_att: Chandra attitude [ra, dec, roll] (deg)
    :param ngrid: number of points on visible Earth to raytrace
    :param planes: list of Plane objects defining structure
    :param p_radiators: points on radiator surface

    :returns: relative visibility
    """
    # Calculate position of earth in ECI and Chandra body coords.  
    # Quat([1,0,0,0]) is a 180 deg roll that puts -Z "up"
    q_att = Quat(chandra_att) * Quat([1,0,0,0.])
    # For T = attitude transformation matrix then p_body = T^-1 p_eci
    p_earth_body = np.dot(q_att.transform.transpose(), -p_chandra_eci)
    # Quaternion to transform x-axis to the body-earth vector
    q_earth = quat_x_v2(p_earth_body)
    open_angle = np.arcsin(Rad_Earth / np.sqrt(np.sum(p_earth_body**2)))
    rays, grid_area = sphere_grid(ngrid, open_angle)
    rays_to_earth = np.dot(q_earth.transform, rays.transpose()).transpose()

    # points that are visible to ACIS radiator
    illum = 0.
    visible = np.zeros((len(rays),), int)

    for p_radiator in p_radiators:
        for i_ray, ray in enumerate(rays_to_earth + p_radiator):
            line = Line(p_radiator, ray)
            for plane in planes:
                if plane_line_intersect(plane, line):
                    # Blocked by SIM structure (aka space taco)
                    break
            else:
                # calculate projection of radiator normal [0,0,1] onto the vector
                # from radiator to ray.  
                # visible[i_ray] += abs(line.u[2])
                visible[i_ray] += 1

    illum = grid_area * np.sum(visible) / (len(p_radiators) * len(rays))
    return visible, illum, rays

