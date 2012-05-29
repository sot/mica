#!/usr/bin/env python

import logging
import obsid_archive


# borrowed from eng_archive
archfiles_hdr_cols = ('tstart', 'tstop', 'caldbver', 'content',
                      'ascdsver', 'revision', 'date')

#config = ConfigObj('asp1.conf')
config = dict(data_root='/data/aca/archive/asp1',
              temp_root='/data/aca/archive/temp',
              sql_def='archfiles_asp_l1_def.sql',
              apstat_table='aspect_1',
              apstat_id='aspect_1_id',
              label='asp_l1',
              small='asp1{fidprops}',
              small_glob='*fidpr*',
              small_ver_regex='pcadf\d+N(\d{3})_',
              full='asp1')


def get_options():
    import argparse
    parser = argparse.ArgumentParser(
        description="Fetch aspect level 1 products and make a file archive")
    defaults = dict(config)
    parser.set_defaults(**defaults)
    parser.add_argument("--obsid",
                        type=int,
                        help="specific obsid to process")
    parser.add_argument("--version",
                        default='last',
                        help="specific processing version to retrieve")
    parser.add_argument("--firstrun",
                        action='store_true',
                        help="for archive init., ignore rev in aspect_1 table")
    parser.add_argument("--data-root",
                        help="parent directory for all data")
    parser.add_argument("--temp-root",
                        help="parent temp directory")
    parser.add_argument("--proctype",
                        default="asp1")
    opt = parser.parse_args()
    return opt


def get_dir(obsid):
    archive = obsid_archive.ObsArchive(dict(config))
    return archive.get_dir(obsid)


def get_obs_dirs(obsid):
    archive = obsid_archive.ObsArchive(dict(config))
    return archive.get_obs_dirs(obsid)


def main():
    opt = get_options()
    config = dict(opt.__dict__, cols=archfiles_hdr_cols)
    archive = obsid_archive.ObsArchive(config)
    archive.logger.setLevel(logging.INFO)
    archive.logger.addHandler(logging.StreamHandler())
    archive.update()

if __name__ == '__main__':
    main()
