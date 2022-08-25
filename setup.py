# Licensed under a 3-clause BSD style license - see LICENSE.rst
from setuptools import setup
try:
    from testr.setup_helper import cmdclass
except ImportError:
    cmdclass = {}
import glob
import sys
import os

scripts = glob.glob("scripts/*")

if "--user" not in sys.argv:
    share_path = os.path.join("share", "acis_taco")
    data_files = [(share_path, ['share/make_esaview_data.py',
                                'task_schedule.cfg', 'task_schedule_occ.cfg'])]
else:
    data_files = None

setup(name='acis_taco',
      author='Tom Aldcroft',
      description='ACIS Earth Solid Angle package',
      author_email='taldcroft@cfa.harvard.edu',
      url='http://cxc.cfa.harvard.edu/mta/ASPECT/tool_doc/acis_taco/',
      use_scm_version=True,
      setup_requires=['setuptools_scm', 'setuptools_scm_git_archive'],
      zip_safe=False,
      packages=['acis_taco', 'acis_taco.tests'],
      scripts=scripts,
      include_package_data=True,
      cmdclass=cmdclass,
      data_files=data_files,
      tests_require=['pytest'],
      )


