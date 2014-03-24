# Set the task name
TASK = mica

# Uncomment the correct choice indicating either SKA or TST flight environment
FLIGHT_ENV = SKA

include /proj/sot/ska/include/Makefile.FLIGHT

# Documentation build directory for installation to INSTALL_DOC
DOC = doc/_build/html

SHARE = scripts/update_aca_l0.py scripts/update_asp_l1.py \
	scripts/update_obspar.py scripts/update_vv.py scripts/update_starchecks.py \
	scripts/update_reports.py scripts/update_aca_dark.py
DATA = task_schedule.cfg \
	mica/archive/obspar_def.sql mica/archive/archfiles_asp_l1_def.sql \
	mica/archive/archfiles_aca_l0_def.sql \
	mica/archive/processing_asp_l1_def.sql \
	mica/starcheck/starcheck.sql mica/report/report_processing.sql

TEMPLATES = templates/report.html templates/aiprops.html \
	templates/props.html templates/vv.html \
	templates/vv_slots_single.html \
	templates/star.html

install: install_doc
	mkdir -p $(INSTALL_SHARE)
	mkdir -p $(INSTALL_DATA)
	mkdir -p $(INSTALL_DATA)/templates/
	rsync --times --cvs-exclude $(DATA) $(INSTALL_DATA)/
	rsync --times --cvs-exclude $(SHARE) $(INSTALL_SHARE)/
	rsync --times --cvs-exclude $(TEMPLATES) $(INSTALL_DATA)/templates/

install_doc:
	mkdir -p $(INSTALL_DOC)
	rsync --archive --times $(DOC)/ $(INSTALL_DOC)/
