from __future__ import print_function

# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np
import taco

p_chandra_eci = np.array([5000., 5000., 0.]) * 1000
att = [0, 90, 0]

taco.set_random_salt(None)
illum1 = taco.calc_earth_vis(p_chandra_eci, att)[1]
illum2 = taco.calc_earth_vis(p_chandra_eci, att)[1]
assert not np.all(illum1 == illum2)

taco.set_random_salt(1)
illum1 = taco.calc_earth_vis(p_chandra_eci, att)[1]
illum2 = taco.calc_earth_vis(p_chandra_eci, att)[1]
assert np.all(illum1 == illum2)

taco.set_random_salt(2)
illum2 = taco.calc_earth_vis(p_chandra_eci, att)[1]
assert not np.all(illum1 == illum2)
