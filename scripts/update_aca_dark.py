#!/usr/bin/env python
"""
Update the Mica ACA dark calibrations database.  This is normally run as a cron job.
"""
if __name__ == '__main__':
    from mica.aca_dark import update_aca_dark
    update_aca_dark.main()
