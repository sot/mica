from mica.archive import aca_hdr3
from astropy.table import Table

mem_loc = sorted(aca_hdr3.HDR3_DEF)
desc = ['{}'.format(aca_hdr3.HDR3_DEF[m_id]['desc'])
        for m_id in mem_loc]
longdesc = ['<pre>{}</pre>'.format(aca_hdr3.HDR3_DEF[m_id]['longdesc'].strip())
            for m_id in mem_loc]

t = Table([mem_loc, desc, longdesc], names=['location', 'desc', 'longdesc'])

from astropy.table.pprint import _pformat_table

ftable, nlines =  _pformat_table(t, max_lines=-1, max_width=-1, html=True)
f = open('hdr3_only_msids.html')
f.writelines(ftable)
f.close()
