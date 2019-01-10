# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np
import acis_taco
from numpy.testing import assert_allclose

p_chandra_eci = np.array([5000., 5000., 0.]) * 1000
att = [0, 90, 0]

illum_gold = np.array([7.91359573e-01, 4.86349158e-01, 2.39867755e-01,
                       1.00615970e-01, 4.47565800e-02, 1.20143256e-02, 
                       2.88178250e-03, 1.21989573e-03, 1.07290480e-04,   
                       0.00000000e+00, 0.00000000e+00])

def test_illum_regress():
    acis_taco.set_random_salt(1)
    illum = acis_taco.calc_earth_vis(p_chandra_eci, att)[1]
    assert_allclose(illum, illum_gold)


