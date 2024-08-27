# Licensed under a 3-clause BSD style license - see LICENSE.rst
import re
import sys

sys.path.insert(0, '..')

from mica.archive import aca_hdr3
from astropy.table import Table
from astropy.table.pprint import _pformat_table

hdr3 = aca_hdr3.HDR3_DEF
m_ids = sorted(hdr3)
msids = [
    '<pre style="font-family:serif">{}</pre>'.format(hdr3[m_id]['msid'])
    for m_id in m_ids
]
descs = [
    '<pre style="font-family:serif">{}</pre>'.format(hdr3[m_id]['longdesc'].strip())
    for m_id in m_ids
]

t = Table([msids, descs], names=['MSID', 'Description'])
ftable, nlines = _pformat_table(t, max_lines=-1, max_width=-1, html=True)

ftable[0] = re.sub(r'<table>', '<table border="1">', ftable[0])

f = open('hdr3_only_msids.html', 'w')
f.writelines(ftable)
f.close()
