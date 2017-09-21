#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Update the Mica ACA dark calibrations database.  This is normally run as a cron job.
"""
if __name__ == '__main__':
    from mica.archive.aca_dark import update_aca_dark
    update_aca_dark.main()
