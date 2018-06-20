# Licensed under a 3-clause BSD style license - see LICENSE.rst
from setuptools import setup
from acis_taco import __version__
try:
    from testr.setup_helper import cmdclass
except ImportError:
    cmdclass = {}
import glob
import sys
import os

scripts = glob.glob("scripts/*")

if "--user" not in sys.argv:
    share_path = os.path.join(sys.prefix, "share", "acis_taco")
    data_path = os.path.join(sys.prefix, "data", "acis_taco")

    data_files = [(share_path, ['share/make_esaview_data.py']),
                  (data_path, ['task_schedule.cfg', 'task_schedule_occ.cfg'])]
else:
    data_files = None

setup(name='acis_taco',
      author='Tom Aldcroft',
      description='ACIS Earth Solid Angle package',
      author_email='aldcroft@head.cfa.harvard.edu',
      url='http://cxc.cfa.harvard.edu/mta/ASPECT/tool_doc/taco/',
      version=__version__,
      zip_safe=False,
      packages=['acis_taco', 'acis_taco.tests'],
      scripts=scripts,
      include_package_data=True,
      cmdclass=cmdclass,
      data_files=data_files,
      tests_require=['pytest'],
      )


