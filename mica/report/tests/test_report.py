import getpass

import pytest

import mica.report.report

user = getpass.getuser()

try:
    import ska_dbi.sqsh

    with ska_dbi.sqsh.Sqsh(
        server="sqlsao", dbi="sybase", user=user, database="axafvv"
    ) as db:
        HAS_SYBASE_ACCESS = True
except Exception:
    HAS_SYBASE_ACCESS = False


@pytest.mark.skipif(
    "not HAS_SYBASE_ACCESS", reason="Report test requires Sybase VV access"
)
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


@pytest.mark.skipif(
    "not HAS_SYBASE_ACCESS", reason="Report test requires Sybase VV access"
)
def test_target_summary_er():
    """
    Test that target_summary for an ER obsid returns None"""
    obsid = 54000
    summary = mica.report.report.target_summary(obsid)
    assert summary is None
