.PHONY: all clean clean-all libcm tabasco

all: tabasco

libcm:
	${MAKE} -C CoEVP $@

tabasco: libcm
	${MAKE} -C charm++

clean:
	${MAKE} -C CoEVP clean
	${MAKE} -C charm++ clean

clean-all: clean
	${MAKE} -C CoEVP clean-all
