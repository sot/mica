from distutils.core import setup

from mica.version import version

license = """\
New BSD/3-clause BSD License
Copyright (c) 2012 Smithsonian Astrophysical Observatory
All rights reserved."""

setup(name='mica',
      description='Mica aspects archive',
      version=str(version),
      author='Jean Connelly',
      author_email='jconnelly@cfa.harvard.edu',
      license=license,
      zip_safe=False,
      packages=['mica', 'mica.archive', 'mica.archive.aca_dark', 'mica.vv',
                'mica.starcheck', 'mica.catalog', 'mica.report'],
      )
