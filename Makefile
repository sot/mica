# Set the task name
TASK = mica

# Uncomment the correct choice indicating either SKA or TST flight environment
FLIGHT_ENV = SKA

include /proj/sot/ska/include/Makefile.FLIGHT

SHARE = scripts/update_aca_l0.py scripts/update_asp_l1.py scripts/update_obspar.py
DATA = task_schedule.cfg \
	mica/archive/obspar_def.sql mica/archive/archfiles_asp_l1_def.sql \
	mica/archive/archfiles_aca_l0_def.sql


install:
	mkdir -p $(INSTALL_SHARE)
	mkdir -p $(INSTALL_DATA)
	rsync --times --cvs-exclude $(DATA) $(INSTALL_DATA)/
	rsync --times --cvs-exclude $(SHARE) $(INSTALL_SHARE)/

