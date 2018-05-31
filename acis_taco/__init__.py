# Licensed under a 3-clause BSD style license - see LICENSE.rst
from .acis_taco import *
from .antisun import AntiSun
from .version import version as __version__

def test(*args, **kwargs):
    """
    Run py.test unit tests.
    """
    import testr
    return testr.test(*args, **kwargs)