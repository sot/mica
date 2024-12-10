import mica.report.report


def test_target_summary_or():
    """
    Test the target_summary method for an OR.

    This test is for obsid 2121 which is quite historical at this point and should
    not change."""
    obsid = 2121
    summary = mica.report.report.target_summary(obsid)
    assert summary is not None
    assert isinstance(summary, dict)
    assert summary["prop_num"] == 2700413
    assert summary["lts_lt_plan"] is None
    assert summary["soe_st_sched_date"] == "Nov 14 2000 12:49AM"

def test_target_summary_er():
    """
    Test that target_summary for an ER obsid returns None"""
    obsid = 54000
    summary = mica.report.report.target_summary(obsid)
    assert summary is None

