# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np
import acis_taco
from numpy.testing import assert_allclose

p_chandra_eci = np.array([5000., 5000., 0.]) * 1000
att = [0, 90, 0]


def test_illum_regress(save=False):
    acis_taco.set_random_salt(1)
    illum = acis_taco.calc_earth_vis(p_chandra_eci, att)[1]
    if save:
        np.save("illum.npy", illum)
    else:
        illum_gold = np.load("illum.npy")
    assert_allclose(illum, illum_gold)


if __name__ == "__main__":
    test_illum_regress(save=True)
