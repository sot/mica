# Licensed under a 3-clause BSD style license - see LICENSE.rst
from setuptools import setup

from mica.version import version

try:
    from testr.setup_helper import cmdclass
except ImportError:
    cmdclass = {}

license = """\
New BSD/3-clause BSD License
Copyright (c) 2012 Smithsonian Astrophysical Observatory
All rights reserved."""

try:
    from testr.setup_helper import cmdclass
except ImportError:
    cmdclass = {}

setup(name='mica',
      description='Mica aspects archive',
      version=version,
      author='Jean Connelly',
      author_email='jconnelly@cfa.harvard.edu',
      license=license,
      zip_safe=False,
      packages=['mica', 'mica.archive', 'mica.archive.tests', 'mica.archive.aca_dark',
                'mica.vv', 'mica.vv.tests',
                'mica.starcheck', 'mica.starcheck.tests', 'mica.catalog', 'mica.report', 'mica.report.tests', 'mica.web',
                'mica.stats', 'mica.stats.tests'],
      package_data={'mica.web': ['templates/*/*.html', 'templates/*.html'],
                    'mica.report': ['templates/*.html']},
      tests_require=['pytest'],
      cmdclass=cmdclass,
      )
