# Licensed under a 3-clause BSD style license - see LICENSE.rst
def make_slab():
    x = [0.1,    0.1,   -1.0, -1.0]
    y = [10.0,  -10.,    10,   -10 ]
    z = [-0.5,  0.5]
    p = np.array([zip(x, y, repeat(z[0])),
                  zip(x, y, repeat(z[1]))])
    planes = [Plane(p[0,2], p[0,0], p[0,1]),
              Plane(p[0,2], p[0,3], p[0,1]),
              Plane(p[1,2], p[1,0], p[1,1]),
              Plane(p[1,2], p[1,3], p[1,1]),
              Plane(p[0,0], p[1,0], p[1,1]),
              Plane(p[0,0], p[0,1], p[1,1]),
              ]
    return planes

def make_radiator2():
    return np.array([[-0.4, 0., 0.]])

