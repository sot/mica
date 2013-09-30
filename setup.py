from setuptools import setup, find_packages
from mica.version import version
setup(name='mica',
      version=str(version),
      author='Jean Connelly',
      author_email='jconnelly@cfa.harvard.edu',
      license="""New BSD/3-clause BSD License
Copyright (c) 2012 Smithsonian Astrophysical Observatory
All rights reserved.""",
      zip_safe=False,
      packages=['mica', 'mica.archive', 'mica.vv',
                'mica.starcheck', 'mica.catalog',
                'mica.report'],
      )
