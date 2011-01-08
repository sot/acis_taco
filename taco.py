"""
Define geometry for ACIS radiator and sunshade and perform raytrace
calculation of Earth illumination on the radiator surface.
"""
import os

from itertools import repeat
import numpy as np
from Quaternion import Quat
import scipy.weave

Rad_Earth = 6371e3

#double sheild_y[5] = {-0.689600,   -0.333000,   -0.0635000,     0.553000,     0.689600};
#double sheild_z[5] = {0.00000,     0.777000,     0.777000,     0.421000,      0.00000};

def make_taco(refl=0):
    """Define geometry for ACIS radiator sun-shade (aka the space taco)."""
    x = [-0.25 * (refl+1),  0.25 * (refl+1)]
    y = [-0.689,  -0.333,  -0.0635,  0.553, 0.6896]
    z = [0.0, 0.777, 0.777, 0.421, 0.0]
    p = np.array([zip(repeat(x[0]), y, z),
                  zip(repeat(x[1]), y, z)])
    planes = [Plane(p[0,0], p[0,1], p[0,2]),
              Plane(p[0,0], p[0,2], p[0,3]),
              Plane(p[0,0], p[0,3], p[0,4]),
              Plane(p[1,0], p[1,1], p[1,2]),
              Plane(p[1,0], p[1,2], p[1,3]),
              Plane(p[1,0], p[1,3], p[1,4]),
              ]
    return planes

def make_radiator(n_radiator_x=2, n_radiator_y=3, nraw=False):
    """Specify points on the ACIS radiator surface.
    Corners of radiator at: (lengths in meters)
    [[-0.209,  0.555, 0.036],
     [ 0.209,  0.555, 0.036],
     [-0.209, -0.389, 0.036],
     [ 0.209, -0.389, 0.036]]
     Center-y is 0.083, length is 0.944.
     Z standoff height is 0.36 
    """
    rad_z = 0.036
    if nraw:
        rad_pts = [[0.0, 0.555, rad_z],
                   [0.0, 0.1, rad_z],
                   [0.0, -0.389, rad_z]]
        rad_pts = [[0.0, 0.0, rad_z]]
    else:
        rad_pts = []
        for x in np.linspace(-0.2, 0.2, n_radiator_x):
            for y in np.linspace(-0.389, 0.555, n_radiator_y):
                rad_pts.append([x, y, rad_z])
    return rad_pts

class Line(object):
    """Line from p0 to p1"""
    def __init__(self, p0, p1):
        self.p0 = np.array(p0)
        self.p1 = np.array(p1)
        self.u = self.p1 - self.p0
        self.len = np.sqrt(np.sum(self.u**2))  # line length
        self.u /= self.len              # unit vector in line direction
        self.p1 = self.p0 + self.u * 10.0

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

def py_plane_line_intersect(p, l):
    """Determine if the line ``l`` intersects the plane ``p``.  Pure python.

    :rtype: boolean
    """
    mat = np.array([l.p0 - l.p1,
                    p.p1 - p.p0,
                    p.p2 - p.p0]).transpose()
    try:
        t, u, v = np.dot(np.linalg.inv(mat), l.p0 - p.p0)
        intersect = 0 <= u <= 1 and 0 <= v <= 1 and u + v <= 1 and t > 0
    except np.linalg.LinAlgError, msg:
        # This presumably occurred because is parallel to plane
        intersect = False

    return intersect

def plane_line_intersect(p, l):
    """Determine if the line ``l`` intersects the plane ``p``.  C code via scipy.weave.inline.

    :rtype: int (0 or 1)
    """
    j = l.p0
    k = l.p1
    A = p.p0
    B = p.p1
    C = p.p2

    return scipy.weave.inline(intersect_code, ['A', 'B', 'C', 'j', 'k'])

def sphere_grid(ngrid, open_angle):
    """Calculate approximately uniform spherical grid of rays containing
    ``ngrid`` points and extending over the opening angle ``open_angle``
    (radians).

    :returns: numpy array of unit length rays, grid area (steradians)
    """
    from math import sin, cos, radians, pi, sqrt
    
    grid_area = 2*pi*(1-cos(open_angle))
    if ngrid <= 1:
        return np.array([[1., 0., 0.]]), grid_area

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
                   ngrid=100,
                   p_radiators=None,
                   to_earth=False,
                   reflect_atten=0.9,
                   n_radiator_x=3,
                   n_radiator_y=4,
                   max_reflect=10):
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
    :param ngrid: number of points on visible Earth to raytrace
    :param p_radiators: points on radiator surface
    
    :returns: relative visibility, total illumination, projected rays
    """
    planes = [make_taco(refl) for refl in range(max_reflect+1)]
    if p_radiators is None:
        p_radiators = make_radiator(n_radiator_x, n_radiator_y)

    # Calculate position of earth in ECI and Chandra body coords.  
    q_att = Quat(chandra_att) 
    
    # For T = attitude transformation matrix then p_body = T^-1 p_eci
    p_earth_body = np.dot(q_att.transform.transpose(), -np.array(p_chandra_eci))

    # Quaternion to transform x-axis to the body-earth vector
    q_earth = quat_x_v2(p_earth_body)
    open_angle = np.arcsin(Rad_Earth / np.sqrt(np.sum(p_earth_body**2)))
    rays, earth_solid_angle = sphere_grid(ngrid, open_angle)
    # or          = np.dot(rays, q_earth.transform.transpose())
    rays_to_earth = np.dot(q_earth.transform, rays.transpose()).transpose()

    vis = np.zeros((max_reflect+1, len(rays)))
    illum = np.zeros(max_reflect+1)
    out_rays = []
    candidates = []
    next_candidates = []

    # Make all the rays
    for p_radiator in p_radiators:
        for i_ray, ray in enumerate(rays_to_earth + p_radiator):
            if rays_to_earth[i_ray][2] > 0:
                candidates.append((i_ray, p_radiator, ray, Line(p_radiator, ray)))
            
    # Main ray-trace loop.  Calculate ray visibility for an increasing number
    # of reflections.  Rays that get blocked after N reflections are candidates
    # for getting out after N+1 reflections.
    for refl in range(max_reflect+1):
        for i_ray, p_radiator, ray, line in candidates:
            for plane in planes[refl]:
                if line.u[0] * plane.p0[0] > 0 and plane_line_intersect(plane, line):  # Blocked
                    out_rays.append((p_radiator, rays_to_earth[i_ray], refl, False))
                    next_candidates.append((i_ray, p_radiator, ray, line))
                    break
            else:
                # Projection of radiator normal [0,0,1] onto vector from radiator to ray.  
                vis[refl, i_ray] += abs(line.u[2])
                out_rays.append((p_radiator, rays_to_earth[i_ray], refl, True))

        vis[refl, :] *= reflect_atten**refl / len(p_radiators)
        illum[refl] = earth_solid_angle * np.sum(vis[refl, :]) / len(rays)

        candidates = next_candidates
        next_candidates = []

    return vis, illum, out_rays

def plot_vis_image(rays, vis, title=''):
    import matplotlib.pyplot as plt
    plt.figure(1, figsize=(5,5))
    plt.clf()
    blocked = vis < 1e-5
    if len(blocked.nonzero()) > 0:
        pos = rays[blocked]
        plt.plot(pos[:,1], pos[:,2], '.k')

    for alpha, ray in zip(vis[~blocked], rays[~blocked]):
        plt.plot([ray[1]], [ray[2]], '.r', alpha=alpha)

    plt.xlim(-1.1, 1.1)
    plt.ylim(-1.1, 1.1)
    plt.title(title)

intersect_code = """
double a[3], b[3], c[3], in_det;
static double R[4][4], Rp[4];

/* a = B - A */
a[0] = B[0] - A[0]; 
a[1] = B[1] - A[1]; 
a[2] = B[2] - A[2];
/* b = C - B */
b[0] = C[0] - A[0];
b[1] = C[1] - A[1];
b[2] = C[2] - A[2];
/* c = a &times; b */
c[0] = a[1] * b[2] - a[2] * b[1];
c[1] = a[2] * b[0] - a[0] * b[2];
c[2] = a[0] * b[1] - a[1] * b[0];
 
/* M^(-1) = (1/det(M)) * adj(M) */
in_det = 1 / (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
R[0][0] = (b[1] * c[2] - b[2] * c[1]) * in_det;
R[0][1] = (b[2] * c[0] - b[0] * c[2]) * in_det;
R[0][2] = (b[0] * c[1] - b[1] * c[0]) * in_det;
R[1][0] = (c[1] * a[2] - c[2] * a[1]) * in_det;
R[1][1] = (c[2] * a[0] - c[0] * a[2]) * in_det;
R[1][2] = (c[0] * a[1] - c[1] * a[0]) * in_det;
R[2][0] = (c[0]) * in_det;
R[2][1] = (c[1]) * in_det;
R[2][2] = (c[2]) * in_det;

/* O = M^(-1) * A */
R[0][3] = -(R[0][0] * A[0] + R[0][1] * A[1] + R[0][2] * A[2]);
R[1][3] = -(R[1][0] * A[0] + R[1][1] * A[1] + R[1][2] * A[2]);
R[2][3] = -(R[2][0] * A[0] + R[2][1] * A[1] + R[2][2] * A[2]);
 
/* fill in last row of 4x4 matrix */
R[3][0] = R[3][1] = R[3][2] = 0;
R[3][3] = 1;

double J[3], K[3];
static double i[2];
J[2] = R[2][0] * j[0] + R[2][1] * j[1] + R[2][2] * j[2] + R[2][3];
K[2] = R[2][0] * k[0] + R[2][1] * k[1] + R[2][2] * k[2] + R[2][3];

if (J[2] * K[2] >= 0)
  {
    return_val = 0;
  } else
  {
    J[0] = R[0][0] * j[0] + R[0][1] * j[1] + R[0][2] * j[2] + R[0][3];
    K[0] = R[0][0] * k[0] + R[0][1] * k[1] + R[0][2] * k[2] + R[0][3];
    i[0] = J[0] + J[2] * ((K[0] - J[0]) / (J[2] - K[2]));
    if (i[0] < 0 || i[0] > 1)
      {
        return_val = 0;
      } else
      {
        J[1] = R[1][0] * j[0] + R[1][1] * j[1] + R[1][2] * j[2] + R[1][3];
        K[1] = R[1][0] * k[0] + R[1][1] * k[1] + R[1][2] * k[2] + R[1][3];
        i[1] = J[1] + J[2] * ((K[1] - J[1]) / (J[2] - K[2]));
        if (i[1] < 0 || i[1] > 1 || i[0] + i[1] > 1)
          {
            return_val = 0;
          } else 
          {
            return_val = 1;
          }
      }
  }
"""
