# Set the task name
TASK = taco

FLIGHT_ENV = SKA

SRC = esaview.py antisun.py
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
	mkdir -p $(INSTALL_DOC)
	mkdir -p $(INSTALL_SHARE)
	rsync -av doc/_build/html/ $(INSTALL_DOC)/
	cp -p dist/$(TASK)-*.tar.gz $(SKA)/export/
	cp esaview $(INSTALL_BIN)/
	cp $(SRC) $(INSTALL_SHARE)/
