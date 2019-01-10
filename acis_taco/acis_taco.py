# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Define geometry for ACIS radiator and sunshade and perform raytrace
calculation of Earth illumination on the radiator surface.
"""
import hashlib
import functools

from Quaternion import Quat
import numpy as np
from numpy.random import RandomState

prng = RandomState()

_RANDOM_SALT = 1

def set_random_salt(salt):
    """
    Set the internal ``_RANDOM_SALT`` module attribute to ``salt``.

    If set to ``None`` then the system random seed will be used everywhere and
    the results of ``calc_earth_vis`` will not be deterministic (i.e. they will
    vary slightly from one run to the next).

    If set otherwise then the ``calc_earth_vis`` function will return precisely
    the same results for the same input arguments and same value of ``salt``.
    Different ``salt`` values will generate different random variations of the
    results.  This is done by making an MD5 digest of random salt and the
    function arguments and using this to generate an integer which sets the
    random seed.

    Any object type is allowed for ``salt``, the only thing that matters is
    is that it has a unique string representation.

    :param salt: any object with a unique str representation, repr(salt)
    """
    global _RANDOM_SALT
    _RANDOM_SALT = salt

def _encode_data(val):
    return np.array_str(np.array(val, ndmin=1)).encode("utf8")

def _make_reproducible(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        if _RANDOM_SALT is not None:
            md5 = hashlib.md5()
            md5.update(repr(_RANDOM_SALT).encode("utf8"))
            for val in args:
                md5.update(_encode_data(val))
            for val in kwargs.values():
                md5.update(_encode_data(val))
            digest = md5.hexdigest()
            seed = int(digest[:7], 16)
            prng.seed(seed)

        out = func(*args, **kwargs)
        return out
    return wrapper


def make_taco():
    """Define geometry for ACIS radiator sun-shade (aka the space taco) in mm."""
    y_pnts = np.array([-689.0, -333, -63, 553, 689]) + TACO_Y_OFF
    z_pnts = np.array([0.0, 777, 777, 421, 0]) - RAD_Z_OFF

    y_edge = np.arange(max(y_pnts))
    z_edge = interpolate(z_pnts, y_pnts, y_edge)

    return z_edge


@_make_reproducible
def calc_earth_vis(p_chandra_eci,
                   chandra_att,
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
    
    :returns: relative visibility, total illumination, projected rays
    """


    # Calculate position of earth in ECI and Chandra body coords.  
    q_att = Quat(chandra_att) 
    
    # For T = attitude transformation matrix then p_body = T^-1 p_eci
    p_earth_body = np.dot(q_att.transform.transpose(), -np.array(p_chandra_eci))

    illum = np.zeros(max_reflect+1)
    out_rays = []

    # If Earth disk is entirely below the "horizon" of the ACIS radiator then illum=0
    if p_earth_body[2] < -RAD_EARTH:
        return [], illum, out_rays

    # Quaternion to transform x-axis to the body-earth vector
    q_earth = quat_x_v2(p_earth_body)
    open_angle = np.arcsin(RAD_EARTH / np.sqrt(np.sum(p_earth_body**2)))
    
    rays, earth_solid_angle = sphere_rand(open_angle)
    n_rays = len(rays)
    rays_to_earth = np.dot(q_earth.transform, rays.transpose()).transpose()  # shape (n_rays, 3)

    # Accept only rays with a positive Z component and make sure no X component is < 1e-6
    rays_to_earth = rays_to_earth[rays_to_earth[:, 2] > 0.0, :]
    n_rays_to_earth = len(rays_to_earth)
    if n_rays_to_earth % 2 == 1:  # Ugh, make sure this list has even length
        rays_to_earth = rays_to_earth[:-1, :]
        n_rays_to_earth -= 1

    vis = np.zeros((max_reflect+1, n_rays))
    if len(rays_to_earth) == 0:
        return vis, illum, out_rays

    rays_to_earth_x = np.abs(rays_to_earth[:, 0])
    rays_to_earth_x[rays_to_earth_x < 1e-6] = 1e-6

    # Main ray-trace loop.  Calculate ray visibility for an increasing number
    # of reflections.  Rays that get blocked after N reflections are candidates
    # for getting out after N+1 reflections.
    
    # Radiator size: 17 inches wide (X) by 19 inches long (Y), centered within
    # SIM shaded structure.

    rad_x = prng.uniform(low=0.0, high=215.9, size=n_rays_to_earth // 2)
    rad_y = prng.uniform(low=-241.3, high=241.3, size=n_rays_to_earth // 2) + TACO_Y_OFF
    rad_x = np.append(rad_x, -rad_x)
    rad_y = np.append(rad_y, rad_y)
    rad_z = np.zeros(n_rays_to_earth)
    rays_x = rays_to_earth_x
    rays_y = rays_to_earth[:, 1]
    rays_z = rays_to_earth[:, 2]

    for refl in range(max_reflect + 1):
        if len(rays_x) == 0:
            break
        taco_x = TACO_X_OFF * (refl + 1)
        dx = (taco_x - rad_x) / rays_x

        # Find rays that intersect shade within Y limits of the shade (0, N_TACO mm)
        ray_taco_iys = np.array(rad_y + rays_y * dx, dtype=np.int)
        y_ok = (ray_taco_iys >= 0) & (ray_taco_iys < N_TACO)
        y_not_ok = ~y_ok
        ray_taco_iys[y_not_ok] = 0

        # From those rays find ones below the TACO_Z_EDGES curve
        ray_taco_zs = rays_z * dx
        z_ok = ray_taco_zs > TACO_Z_EDGES[ray_taco_iys]
        y_z_ok = (y_ok & z_ok) | y_not_ok

        # Calculate the delta-visibility for good rays (cos(incidence_angle) = Z component)
        d_vis = rays_z[y_z_ok] * REFLECT_ATTEN**refl
        # vis[refl, rays_i[y_z_ok]] += d_vis
        illum[refl] += earth_solid_angle * np.sum(d_vis) / n_rays

        blocked = ~y_z_ok
        rays_x = rays_x[blocked]
        rays_y = rays_y[blocked]
        rays_z = rays_z[blocked]
        # rays_i = rays_i[blocked]
        rad_x = rad_x[blocked]
        rad_y = rad_y[blocked]

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

    i1 = np.searchsorted(xin, xout)
    i1[ i1==0 ] = 1
    i1[ i1==lenxin ] = lenxin-1

    x0 = xin[i1-1]
    x1 = xin[i1]
    y0 = yin[i1-1]
    y1 = yin[i1]

    return (xout - x0) / (x1 - x0) * (y1 - y0) + y0

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


@_make_reproducible
def sphere_rand(open_angle, min_ngrid=100, max_ngrid=10000):
    """Calculate approximately uniform spherical grid of rays containing
    ``ngrid`` points and extending over the opening angle ``open_angle``
    (radians).

    :returns: numpy array of unit length rays, grid area (steradians)
    """
    from math import sin, cos, radians, pi, sqrt
    
    xmin = cos(open_angle)
    grid_area = 2*pi*(1-xmin)

    idx_xmin = np.searchsorted(SPHERE_X, xmin)
    if N_SPHERE - idx_xmin < min_ngrid:
        idx_xmin = N_SPHERE - min_ngrid

    ngrid = int(min_ngrid + (max_ngrid-min_ngrid) * (1-xmin))
    idx_sphere = prng.randint(idx_xmin, N_SPHERE, ngrid)

    return SPHERE_XYZ[idx_sphere, :], grid_area
    

@_make_reproducible
def random_hemisphere(nsample):
    x = prng.uniform(low=0.3, high=1.0, size=nsample)
    x.sort()                    # x is not random
    t = 2*np.pi * prng.uniform(size=nsample)
    r = np.sqrt(1-x**2)
    z = r * np.cos(t)
    y = r * np.sin(t)
    return np.array([x, y, z]).transpose()


RAD_EARTH = 6371e3
RAD_Z_OFF = 36
TACO_X_OFF = 230
TACO_Y_OFF = 689

REFLECT_ATTEN = 0.9
TACO_Z_EDGES = make_taco()
N_TACO = len(TACO_Z_EDGES)

N_SPHERE = 1500000
SPHERE_XYZ = random_hemisphere(N_SPHERE)
SPHERE_X = SPHERE_XYZ[:, 0].copy()
