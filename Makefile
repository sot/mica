# Set the task name
TASK = mica

# Uncomment the correct choice indicating either SKA or TST flight environment
FLIGHT_ENV = SKA

include /proj/sot/ska/include/Makefile.FLIGHT

# Documentation build directory for installation to INSTALL_DOC
DOC = doc/_build/html

SHARE = scripts/update_aca_l0.py scripts/update_asp_l1.py \
	scripts/update_obspar.py scripts/update_vv.py scripts/update_starchecks.py \
	scripts/update_reports.py scripts/update_aca_dark.py \
	scripts/update_acq_stats.py scripts/update_guide_stats.py
DATA = task_schedule.cfg \
	mica/archive/obspar_def.sql mica/archive/archfiles_asp_l1_def.sql \
	mica/archive/archfiles_aca_l0_def.sql \
	mica/archive/processing_asp_l1_def.sql \
	mica/starcheck/starcheck.sql mica/report/report_processing.sql

install: install_doc
	mkdir -p $(INSTALL_SHARE)
	mkdir -p $(INSTALL_DATA)
	rsync --times --cvs-exclude $(DATA) $(INSTALL_DATA)/
	rsync --times --cvs-exclude $(SHARE) $(INSTALL_SHARE)/

install_doc:
	mkdir -p $(INSTALL_DOC)
	rsync --archive --times $(DOC)/ $(INSTALL_DOC)/
