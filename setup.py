# Licensed under a 3-clause BSD style license - see LICENSE.rst
import sys
import os
from setuptools import setup
from pathlib import Path

license = """\
New BSD/3-clause BSD License
Copyright (c) 2012 Smithsonian Astrophysical Observatory
All rights reserved."""

try:
    from testr.setup_helper import cmdclass
except ImportError:
    cmdclass = {}

entry_points = {
    'console_scripts': [
        'mica_update_ocat_target_table=mica.archive.cda.services:main_update_ocat_local'
    ]
}

# Borrowed from acis_taco setup.py.  This will install all of the following
# into sys.prefix/share/mica/ via the data_files directive.
if "--user" not in sys.argv:
    share_path = os.path.join("share", "mica")
    scripts = [str(script) for script in Path('scripts').glob('update_*.py')]
    data_files = [(share_path, scripts + ['task_schedule.cfg'])]
else:
    data_files = None

setup(
    name='mica',
    description='Mica aspects archive',
    use_scm_version=True,
    setup_requires=['setuptools_scm', 'setuptools_scm_git_archive'],
    author='Jean Connelly',
    author_email='jconnelly@cfa.harvard.edu',
    license=license,
    zip_safe=False,
    include_package_data=True,
    data_files=data_files,
    packages=[
        'mica',
        'mica.archive',
        'mica.archive.tests',
        'mica.archive.aca_dark',
        'mica.vv',
        'mica.vv.tests',
        'mica.starcheck',
        'mica.starcheck.tests',
        'mica.catalog',
        'mica.report',
        'mica.report.tests',
        'mica.web',
        'mica.stats',
        'mica.stats.tests',
    ],
    package_data={
        'mica.web': ['templates/*/*.html', 'templates/*.html'],
        'mica.report': ['templates/*.html', '*.sql'],
        'mica.archive': ['*.sql'],
        'mica.starcheck': ['*.sql'],
    },
    entry_points=entry_points,
    tests_require=['pytest'],
    cmdclass=cmdclass,
)
