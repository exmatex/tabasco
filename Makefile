.PHONY: all clean clean-all libcm tabasco

FLANN=no
REDIS=no
SILO=no
ifeq ($(SILO),yes)
SILO_LOC=../CoEVP/silo/silo
endif

all: tabasco

libcm:
	${MAKE} -C CoEVP FLANN=$(FLANN) REDIS=$(REDIS) SILO=$(SILO)

tabasco: libcm
	${MAKE} -C charm++ SILO_LOC=$(SILO_LOC)

clean:
	${MAKE} -C CoEVP clean
	${MAKE} -C charm++ clean

clean-all: clean
	${MAKE} -C CoEVP clean-all
