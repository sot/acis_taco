# Set the task name
TASK = taco

FLIGHT_ENV = SKA

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
	rsync -av doc/_build/html/ $(INSTALL_DOC)/
	cp -p dist/$(TASK)-*.tar.gz $(SKA)/export/

