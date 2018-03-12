#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np
import Quaternion
from itertools import repeat
import matplotlib.pyplot as plt

Rad_Earth = 6371. 

def make_taco():
    x = [0.0, -0.6, -0.5, -0.3, 0.0]
    y = [0.0,  0.2,  0.5,  0.8, 1.0]
    z = [0.1,  0.5]
    p = np.array([zip(x, y, repeat(z[0])),
                  zip(x, y, repeat(z[1]))])
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
    return np.array([[-.1, 0.5, 0.3]])

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

def sphere_grid(pix_size):
    """Calculate approximately uniform spherical grid with spacing ``pix_size``"""
    from math import sin, cos, radians, pi
    
    grid = []
    n_d = int(round(pi/2 / radians(pix_size)))
    d_d = pi/2 / n_d

    for i_d in range(-n_d, n_d+1):
        dec = i_d * d_d
        if abs(i_d) != n_d:
            n_r = int(round( 2*pi * cos(dec) / d_d))
            d_r = 2*pi / n_r
        else:
            n_r = 1
            d_r = 1
        for i_r in range(0, n_r):
            ra = i_r * d_r
            grid.append((cos(ra) * cos(dec), sin(ra) * cos(dec), sin(dec)))

    # (x, y, z) = zip(*grid)  (python magic)
    return np.array(grid)

def nearest_sphere_intersect(l, c, r, len):
    l_dot_c = np.dot(l, c)
    det = l_dot_c**2 - np.sum(c**2) + r**2
    if det < 0:
        return False
    sqrtdet = np.sqrt(det)
    d1 = l_dot_c - sqrtdet
    d2 = l_dot_c + sqrtdet
    return abs(d1 - len) <= abs(d2 - len)
    

def calc_earth_vis(p_chandra_eci, chandra_att,
                   planes=make_taco(),
                   p_radiators=make_radiator(),
                   earth_surface_grid=sphere_grid(4)):
    """Calculate the relative Earth visibility for the ACIS radiator given
    the Chandra orbit position ``p_chandra_eci`` and attitude ``chandra_att``.

    The relative visibility is normalized so that 1.0 represents the entire
    radiator seeing the full Earth at 100000 km.

    :param p_chandra_eci: Chandra orbital position [x, y, z] (meters)
    :param chandra_att: Chandra attitude [ra, dec, roll] (deg)

    :returns: relative visibility
    """
    # Calculate position of earth in ECI and Chandra body coords
    # For T = attitude transformation matrix then p_body = T^-1 p_eci
    q_att = Quaternion.Quat(chandra_att)
    p_earth_body = np.dot(q_att.transform.transpose(), -p_chandra_eci)
    p_earth_surfaces = (p_earth_body + earth_surface_grid * Rad_Earth)

    # points that are visible to ACIS radiator
    behind_earth = []
    visible = []
    blocked = []
    illum = 0.

    for p_radiator in p_radiators:
        for p_earth_surface in p_earth_surfaces:
            line = Line(p_radiator, p_earth_surface)
            if nearest_sphere_intersect(line.u, p_earth_body, Rad_Earth, line.len):
                for plane in planes:
                    if plane_line_intersect(plane, line):
                        # Blocked by SIM structure (aka space taco)
                        blocked.append(p_earth_surface)
                        break
                else:
                    visible.append(p_earth_surface)
                    # calculate projection of radiator normal [0,0,1] onto the vector
                    # from radiator to p_earth_surface.  
                    illum += abs(line.u[2]) / line.len**2
                    # Do I need another cos(theta) term for Lambert's law?? YES
            else:
                behind_earth.append(p_earth_surface)

    illum *= 10000.**2 / (len(p_radiators) * len(p_earth_surfaces)) / 0.038
    return visible, blocked, behind_earth, illum



                    
                    

