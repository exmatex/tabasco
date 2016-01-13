.PHONY: all clean clean-all libcm tabasco

FLANN=no
REDIS=no
SILO=yes
ifeq ($(SILO),yes)
SILO_LOC=../CoEVP/silo/silo
endif

all: tabasco

libcm:
	${MAKE} -C CoEVP FLANN=$(FLANN) REDIS=$(REDIS) SILO=$(SILO) CHARM=yes

tabasco: libcm
	${MAKE} -C charm++ SILO_LOC=$(SILO_LOC)

clean:
	${MAKE} -C CoEVP clean
	${MAKE} -C charm++ clean

clean-all: clean
	${MAKE} -C CoEVP clean-all

get_reference:
	${MAKE} -C CoEVP $@

test:
	${MAKE} -C CoEVP FLANN=$(FLANN) REDIS=$(REDIS) SILO=$(SILO) CHARM=yes test
