# script to grab aspect_1_ids for asp1 products in the database and
# make a table to assist in cron job V&V processing

import os
import re
import logging
from astropy.io import fits
import Ska.DBI
from glob import glob
import gzip

mica_archive = os.environ.get('MICA_ARCHIVE') or '/data/aca/archive'

def update(obsids):
    from mica.archive import asp_l1
    """
    For a list of obsids, update the mica table of aspect 1 processing
    """
    logger = logging.getLogger('asp 1 proc table')
    logger.setLevel(logging.INFO)
    logger.addHandler(logging.StreamHandler())

    apstat_db = Ska.DBI.DBI(dbi='sybase',
                            server='sqlsao',
                            database='axafapstat')
    proc_db = Ska.DBI.DBI(dbi='sqlite',
                          server=os.path.join(mica_archive, 'processing_asp_l1.db3'))
    archdb = Ska.DBI.DBI(dbi='sqlite',
                         server=os.path.join(mica_archive, 'asp1', 'archfiles.db3'))
    for obs in obsids:
        logger.info("Adding asp1 processing for obs {}".format(obs))
        asols = asp_l1.get_files(obsid=obs, content='ASPSOL', revision='all')
        for sol in asols:
            logger.info("\tprocessing {}".format(sol))
            procdir = os.path.dirname(sol)
            logfile = glob(os.path.join(procdir, "*log*"))[0]
            aspect_log = gzip.open(logfile).read()
            # read the obspar version with a regex from the log
            obspar_version = int(
                re.search("axaff\d{5}_\d{3}N(\d{3})_obs0a\.par", aspect_log).group(1))
            hdus = fits.open(sol)
            obi = hdus[1].header['OBI_NUM']
            revision = hdus[1].header['REVISION']
            aspect_1 = apstat_db.fetchall("""SELECT * FROM aspect_1
                                             WHERE obsid = {obsid}
                                             AND obi = {obi}
                                             AND revision = {revision}
                                          """.format(obsid=obs,
                                                     obi=obi,
                                                     revision=revision))
            if len(aspect_1) > 1:
                raise ValueError
            indb = proc_db.fetchall("""SELECT * from aspect_1_proc
                                       WHERE aspect_1_id = {}
                                    """.format(aspect_1[0]['aspect_1_id']))
            if len(indb):
                logger.info("\tSkipping; already in table")
                continue
            archrec = archdb.fetchall("""select * from archfiles
                                     where obsid = {obsid}
                                     and revision = {revision}
                                     and content = 'ASPSOL'
                                  """.format(obsid=obs,
                                             revision=revision))
            proc_db.insert(dict(aspect_1_id=aspect_1[0]['aspect_1_id'],
                                obsid=obs,
                                obi=obi,
                                revision=revision,
                                obspar_version=obspar_version,
                                ap_date=str(aspect_1[0]['ap_date'])),
                           'aspect_1_proc')
            isdefault = archrec[0]['isdefault']
            if isdefault is not None:
                proc_db.execute("""UPDATE aspect_1_proc SET isdefault = 1
                                   WHERE obsid = {obsid}
                                   AND revision = {revision}
                                """.format(obsid=obs,
                                           revision=revision))
            logger.info("\tUpdated table for {}".format(obs))
