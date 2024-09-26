# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Update the agasc_supplement.h5.

This file is a supplement to the stable AGASC to inform star selection
and star catalog checking.

Currently this script only has the capability to add a bad star to the
bad star table.  It might end up including functionality to automatically
update another table with effective mags based on acq / guide history.

For process instructions see: https://github.com/sot/mica/wiki/AGASC-supplement
"""

import argparse
import os
from pathlib import Path

import pyyaks.logger
from astropy.table import Table

SKA = Path(os.environ["SKA"])
logger = None  # Set via global in main()


def get_options(args=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-root",
        default=".",
        help=("Directory containing agasc_supplement.h5 (default='.')"),
    )
    parser.add_argument(
        "--bad-star-id", type=int, help="AGASC ID of star to add to bad-star list"
    )
    parser.add_argument(
        "--bad-star-source",
        type=int,
        help=(
            "Source identifier indicating provenance (default=max "
            "existing source + 1)"
        ),
    )
    parser.add_argument(
        "--log-level", default=20, help="Logging level (default=20 (info))"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Dry run (no actual file or database updates)",
    )

    opt = parser.parse_args(args)
    return opt


def main(args=None):
    global logger

    # Setup for updating the sync repository
    opt = get_options(args)

    # Set up logging
    loglevel = int(opt.log_level)
    logger = pyyaks.logger.get_logger(
        name="mica_update_agasc_supplement", level=loglevel, format="%(message)s"
    )

    data_root = Path(opt.data_root)
    suppl_file = data_root / "agasc_supplement.h5"
    if suppl_file.exists():
        logger.info(f"Updating agasc_supplement at {suppl_file}")
    else:
        raise IOError(f"file {suppl_file.absolute()} not found")

    if opt.bad_star_id:
        add_bad_star(opt.bad_star_id, opt.bad_star_source, suppl_file, opt.dry_run)


def add_bad_star(bad_star_id, bad_star_source, suppl_file, dry_run):
    bad_star_id = int(bad_star_id)
    dat = Table.read(str(suppl_file), format="hdf5", path="bad")

    if bad_star_source is None:
        bad_star_source = dat["source"].max() + 1
    else:
        bad_star_source = int(bad_star_source)

    dat.add_row((bad_star_id, bad_star_source))

    logger.info(
        f"Appending {bad_star_id} with source={bad_star_source} to {suppl_file}"
    )
    logger.info("")
    logger.info("IMPORTANT:")
    logger.info("Edit following if source ID is new:")
    logger.info("   https://github.com/sot/mica/wiki/AGASC-supplement")
    logger.info("")
    logger.info("The wiki page also includes instructions for test, review, approval")
    logger.info("and installation.")
    if not dry_run:
        dat.write(
            str(suppl_file), format="hdf5", path="bad", append=True, overwrite=True
        )


if __name__ == "__main__":
    main()
