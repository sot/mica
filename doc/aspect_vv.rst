.. _vv-label:

Post-facto aspect V&V
=====================

The :mod:`mica.vv` module provides tools to create and inspect V&V-type data.

Access Mica V&V Data
--------------------

Locate mica v&v products:

   >>> from mica.vv import get_vv_dir
   >>> get_vv_dir(16504)
   '/data/aca/archive/vv/16/16504_v01'

List files:

   >>> from mica.vv import get_vv_files
   >>> get_vv_files(16504)
   ['/data/aca/archive/vv/16/16504_v01/vv_report.json',
    '/data/aca/archive/vv/16/16504_v01/vv_slot.pkl',
    '/data/aca/archive/vv/16/16504_v01/vv_report.pkl',
    '/data/aca/archive/vv/16/16504_v01/slot_0_yz.png',
    ...

Retrieve all residuals:

   >>> from mica.vv import get_rms_data
   >>> data = get_rms_data()
   >>> data[data['obsid'] == 16505]['dz_rms']
   array([  2.62251067e-05,   5.16913551e-02,   5.32668958e-02,
            4.78861857e-02])

Retrieve mica v&v values for an already-mica-processed obsid.  The dictionaries of
values still need more documentation at this time.

   >>> from mica.vv import get_vv
   >>> obs = get_vv(16504)
   >>> obs['slots']['7']['dz_rms']
   0.11610256063309182

Run mica V&V tools directly
---------------------------

Using the mica-archived aspect solution and obspar, run the mica obsid tools

   >>> from mica.vv import get_arch_vv
   >>> obi = get_arch_vv(2121)

Plot slot 4 residuals

   >>> obi.plot_slot(4)

Look at the SIM drift values

   >>> obi.info()['sim']
   {'max_d_dy': 0.002197265625,
    'max_d_dz': 0.0018472671508789062,
    'max_medf_dy': 3.3234403133392334,
    'max_medf_dz': 7.3717021942138672,
    'min_medf_dy': 3.0496618747711182,
    'min_medf_dz': 6.3389387130737305}

The Obi class can also be called directly on data that isn't in the mica archive. For
example, on c3po-v, something like this could be used to plot residuals on un-ingested
data (these directories will likely not exist to run this example in the future):

   >>> proc_dir = '/dsops/ap/sdp.10/opus/prs_run/done/ASP_L1____502323245n674/'
   >>> aspect_dir = proc_dir  + 'output'
   >>> obspar = proc_dir + 'input/axaff14565_000N001_obs0a.par'
   >>> import mica.vv
   >>> obi = mica.vv.Obi(obspar, aspect_dir)

