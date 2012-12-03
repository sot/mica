from setuptools import setup, find_packages
setup(name='mica',
      version="0.1",
      author='Jean Connelly',
      author_email='jconnelly@cfa.harvard.edu',
      license="""New BSD/3-clause BSD License
Copyright (c) 2012 Smithsonian Astrophysical Observatory
All rights reserved.""",
      zip_safe=False,
      packages=['mica', 'mica.archive', 'mica.quaternion', 'mica.vv'],
      )
