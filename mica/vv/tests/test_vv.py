import os
import numpy as np
from .. import vv
from .. import process
from ... import common


def test_get_vv_dir():
    obsdir = vv.get_vv_dir(16504)
    assert obsdir == os.path.abspath(os.path.join(common.MICA_ARCHIVE, 'vv/16/16504_v01'))


def test_get_vv_files():
    obsfiles = vv.get_vv_files(16504)
    assert sorted(obsfiles)[-1] == os.path.abspath(os.path.join(common.MICA_ARCHIVE,
                                                                'vv/16/16504_v01/vv_report.pkl'))

def test_get_rms_data():
    data = vv.get_rms_data()
    dz_rms = data[(data['obsid'] == 16505) & (data['slot'] == 4) & (data['isdefault'] == 1)]['dz_rms'][0]
    np.testing.assert_allclose(dz_rms, 0.047886185719034906)


def test_get_vv():
    obs = vv.get_vv(16504)
    np.testing.assert_allclose(obs['slots']['7']['dz_rms'], 0.11610256063309182)


def test_run_vv():
    obi = process.get_arch_vv(2121)
    np.testing.assert_allclose(obi.info()['sim']['max_d_dy'], 0.002197265625)

