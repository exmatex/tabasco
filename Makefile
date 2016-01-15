.PHONY: all clean clean-all libcm tabasco test test_CoEVP

FLANN=yes
ifeq ($(FLANN),yes)
FLANN_LOC=../CoEVP/flann/flann/src/cpp
endif
REDIS=yes
ifeq ($(REDIS),yes)
REDIS_LOC=../CoEVP/redis/hiredis
endif
SILO=yes
ifeq ($(SILO),yes)
SILO_LOC=../CoEVP/silo/silo
SILODIFF=CoEVP/silo/silo/bin/silodiff
endif

all: tabasco

libcm:
	${MAKE} -C CoEVP FLANN=$(FLANN) REDIS=$(REDIS) SILO=$(SILO) CHARM=yes

tabasco: charm++/charm_bin/TabaSCo

charm++/charm_bin/TabaSCo: libcm
	${MAKE} -C charm++ FLANN_LOC=$(FLANN_LOC) SILO_LOC=$(SILO_LOC) REDIS_LOC=$(REDIS_LOC)

clean:
	${MAKE} -C CoEVP clean
	${MAKE} -C charm++ clean

clean-all: clean
	${MAKE} -C CoEVP clean-all

get_reference:
	${MAKE} -C CoEVP $@

reference:
	${MAKE} -C CoEVP $@

test_CoEVP:
	${MAKE} -C CoEVP FLANN=$(FLANN) REDIS=$(REDIS) SILO=$(SILO) CHARM=yes test

dummy: ;

TESTFILE=reftest.json

CHARMRUN=../charm++/charm_bin/charmrun ++local ++p 1
test/.charmflags: dummy
	mkdir -p test
	@[ -f $@ ] || touch $@
	@echo "CHARMRUN=$(CHARMRUN)" | cmp -s $@ - || echo "CHARMRUN=$(CHARMRUN)" > $@

TABASCO_OPTS=../charm++/input/$(TESTFILE)
test/.tabascoopts: dummy
	mkdir -p test
	@[ -f $@ ] || touch $@
	@echo "TABASCO_OPTS=$(TABASCO_OPTS)" | cmp -s $@ - || echo "TABASCO_OPTS=$(TABASCO_OPTS)" > $@

STEPS=0500
#bit hackish, but let's assume we have $(STEPS) steps
test/taylor_$(STEPS).silo: charm++/charm_bin/TabaSCo test/.charmflags test/.tabascoopts
	@[ "$(SILO)" = "yes" ] || { echo "make test needs SILO=yes" && exit 1; }
	mkdir -p test
	cd test && $(CHARMRUN) ../charm++/charm_bin/TabaSCo $(TABASCO_OPTS)

SILODIFF_OPTS=-A 1e-8 -E _hdf5libinfo
test: test/taylor_$(STEPS).silo 
	@[ -x "$(SILODIFF)" ] || { echo "SILODIFF=$(SILODIFF) seems to be wrong" && exit 1; }
	$(SILODIFF) ${SILODIFF_OPTS} CoEVP/test/reference test > test/diff
	@[ ! -s test/diff ] || { echo "Difference in files" && head -n 50 test/diff && exit 1; }
