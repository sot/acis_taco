# Set the task name
TASK = taco

FLIGHT_ENV = SKA

DATA = task_schedule.cfg task_schedule_occ.cfg
SRC = esaview.py antisun.py make_esaview_data.py
include /proj/sot/ska/include/Makefile.FLIGHT

.PHONY: dist doc install

# Make a versioned distribution.  Could also use an EXCLUDE_MANIFEST
dist:
	python setup.py sdist

doc:
	cd doc; \
	make html

install:
#  Uncomment the lines which apply for this task
	mkdir -p $(INSTALL_DATA)
	mkdir -p $(INSTALL_DOC)
	mkdir -p $(INSTALL_SHARE)
	rsync -av doc/_build/html/ $(INSTALL_DOC)/
	cp -p dist/$(TASK)-*.tar.gz $(SKA)/export/
	cp -p $(DATA) $(INSTALL_DATA)/
	cp esaview $(INSTALL_BIN)/
	cp -p $(SRC) $(INSTALL_SHARE)/
