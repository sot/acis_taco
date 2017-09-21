# Licensed under a 3-clause BSD style license - see LICENSE.rst
if 0:
    v_ras = []
    v_decs = []
    nv_ras = []
    nv_decs = []

    for ra in range(0, 360, 5):
        for dec in range(-90, 90, 5):
            for roll in range(0, 360, 360):
                q_att = Quat([ra, dec, roll*1.0])
                p_earth_body = np.dot(q_att.transform.transpose(), p_earth_eci)
                taco.calc_earth_vis(p_earth_body, q_att,
                                    planes=slab_planes,
                                    p_radiators=slab_rads)
    plt.clf()
    plt.plot(v_ras, v_decs, '.b')
    plt.plot(nv_ras, nv_decs, '.r')
