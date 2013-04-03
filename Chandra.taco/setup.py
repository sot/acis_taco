from setuptools import setup

from Chandra.taco.version import version

setup(name='Chandra.taco',
      author = 'Tom Aldcroft',
      description='ACIS Earth Solid Angle package',
      author_email = 'aldcroft@head.cfa.harvard.edu',
      url = 'http://cxc.cfa.harvard.edu/mta/ASPECT/tool_doc/taco/',
      version=version,
      zip_safe=False,
      py_modules=['Chandra.taco'],
      namespace_packages=['Chandra'],
      packages=['Chandra', 'Chandra.taco'],
      package_dir={'Chandra':'Chandra'},
      )


