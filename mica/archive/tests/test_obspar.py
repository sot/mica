from .. import obspar

def test_get_obsids():
    """
    Test that archive.obspar.get_obsids gets a reasonable set of obsids.
    """
    start = "2016:001"
    stop = "2016:005"
    mica_obs = obspar.get_obsids(start, stop)
    # from kadi.events.obsids.filter(start, stop)
    kadi_obsid_list = [51365, 18736, 18036, 18203, 18119,
                       18737, 51364, 51363, 51362, 51361,
                       18216, 16821, 17023, 18293, 18623,
                       18248]
    assert len(kadi_obsid_list) == len(mica_obs)
    for o1, o2 in zip(kadi_obsid_list, mica_obs):
        assert o1 == o2
