SRCDIR=src

ifneq ($(shell which charmc), )
CHARM_ROOT=$(dir $(shell which charmc))..
endif

##############CHARM##########
ifeq (,$(CHARM_ROOT))
$(info Please establish Charmm++ environment variables before using this Makefile.)
$(info E.g. by running setting CHARM_ROOT or run 'module load charm++')
$(error CHARM_ROOT is not set)
else
$(info using CHARM_ROOT=$(CHARM_ROOT))
endif

ifneq ($(wildcard $(CHARM_ROOT)/include/*/charm.h), )
CHARMINC=$(dir $(wildcard $(CHARM_ROOT)/include/*/charm.h))
$(info using CHARMINC=$(CHARMINC))
else
CHARMINC=$(CHARM_ROOT)/include/
endif

CHARMBIN=$(CHARM_ROOT)/bin
CXX=$(CHARMBIN)/charmc
OBJDIR=charm_obj
BINDIR=charm_bin

OPTFLAGS= -g -O3 -I../CoEVP -I../CoEVP/CM/include -std=c++0x
#OPTFLAGS= -g -O3 -I../CoEVP -I../CoEVP/CM/include -std=c++11
#other features
#DEFS=-DTRACE

ifneq ($(wildcard /usr/include/boost/foreach.hpp), )
BOOST_INCLUDES=/usr/include
endif

ifeq ($(BOOST_INCLUDES), )
$(error Set BOOST_INCLUDES or run 'module load boost' first)
else
BOOST_CFLAG=-I$(BOOST_INCLUDES)
endif
#We use boost header only so far
BOOSTLIB=$(BOOST_INCLUDES)/../lib
BOOST_LDFLAG=-L$(BOOSTLIB)

LINALG_LDFLAG=$(shell pkg-config --silence-errors --libs lapack blas || true)
ifeq ($(LINALG_LDFLAG), )
LINALG_LDFLAG= -llapack -lblas
endif

OBJS:=$(addprefix $(OBJDIR)/, CoarseScaleModel.o FineScaleModel.o NearestNeighborSearch.o Interpolate.o DBInterface.o DBVecMessage.o DBMap.o ModelDBInterface.o input.o Main.o)
# ../CoEVP/LULESH/lulesh.o)
CXXFLAGS+=$(OPTFLAGS) $(BOOST_CFLAG) $(DEFS)
ifeq ($(NNS_AS_CHARE), yes)
CXXFLAGS += -DNNS_AS_CHARE
endif
ifeq ($(INTERP_AS_CHARE), yes)
CXXFLAGS += -DINTERP_AS_CHARE
endif
ifeq ($(DB_AS_CHARE), yes)
CXXFLAGS += -DDB_AS_CHARE
endif

ifneq ($(strip $(FLANN_LOC)),)
CXXFLAGS += -DFLANN -I$(realpath $(FLANN_LOC)) -I$(realpath $(FLANN_LOC))/flann
endif

ifneq ($(strip $(REDIS_LOC)),)
CXXFLAGS += -DREDIS -I$(realpath $(REDIS_LOC))
endif

#LDFLAGS= -tracemode summary -L../CoEVP/CM/lib -lcm $(LINALG_LDFLAG) -lm -lrt 
#LDFLAGS= -tracemode projections -L../CoEVP/CM/lib -lcm $(LINALG_LDFLAG) -lm -lrt 
LDFLAGS= -L../CoEVP/CM/lib -lcm $(LINALG_LDFLAG) -lm -lrt 
ifneq ($(strip $(SILO_LOC)),)
LDFLAGS += -Wl,-rpath,$(realpath $(SILO_LOC))/lib -L$(realpath $(SILO_LOC))/lib -lsiloh5
endif

LULESH_OBJ=../CoEVP/LULESH/lulesh.o

#target
NAME=$(BINDIR)/TabaSCo
default: all
all: $(SUBDIRS) $(OBJDIR) $(NAME)

dummy: ;

$(OBJDIR)/.cxxflags: dummy
	@[ -f $@ ] || touch $@
	@echo "CXXFLAGS=$(CXXFLAGS)" | cmp -s $@ - || echo "CXXFLAGS=$(CXXFLAGS)" > $@

$(OBJDIR)/%.d: $(SRCDIR)/%.cc
	@#1. sed:  put one file per line * * -> *\\\n*
	@#2. sed gcc -MG does not know that missing files will be in $(SRCDIR)
	@# no path -> SRCDIR
	@echo g++ -MM -MG -MT $(OBJDIR)/$*.o $(CXXFLAGS) -I$(SRCDIR) -I$(CHARMINC) $< \> $@
	@g++ -MM -MG -MT $(OBJDIR)/$*.o $(CXXFLAGS) -I$(SRCDIR) -I$(CHARMINC) $< | \
	sed 's/\([^[:space:]]\) \([^\\[:space:]]\)/\1 \\\n \2/g' | \
	sed '/^[^/]*.\(def\|decl\)\./s@[^[:space:]]@$(SRCDIR)/&@' > $@

#Charmm++ ci files
$(SRCDIR)/%.decl.h $(SRCDIR)/%.def.h: $(SRCDIR)/%.ci
	@#charmc writes to pwd only
	cd $(<D) && $(CXX) $(<F)

DEPS=$(OBJS:.o=.d)
ifneq "$(MAKECMDGOALS)" "clean"
-include $(DEPS)
endif

#make subdirs for objects and executable
$(NAME): | $(BINDIR)

$(BINDIR):
	@mkdir -p $(BINDIR)

$(OBJS): | $(OBJDIR)
$(DEPS): | $(OBJDIR)

$(OBJDIR):
	@mkdir -p $(OBJDIR)

##--- Executable ---##

$(NAME): $(OBJS) $(LULESH_OBJ)
	$(CXX) -o $@ $^ $(LDFLAGS)
	@mv charmrun $(BINDIR)/

#GNU make implicit rule
$(OBJDIR)/%.o: $(SRCDIR)/%.cc $(OBJDIR)/.cxxflags
	$(CXX) -c $(CXXFLAGS) $< -o $@

.PHONY: $(SUBDIRS)
$(SUBDIRS):
	$(MAKE) $(MFLAGS) -C $@

subdirclean:
	@for i in $(SUBDIRS); do \
	  echo $(MAKE) $(MFLAGS) -C $$i clean; \
	  $(MAKE) $(MFLAGS) -C $$i clean || exit 1; \
	done

clean: subdirclean
	rm -f $(SRCDIR)/*.decl.h $(SRCDIR)/*.def.h $(OBJDIR)/charmrun
	rm -f *.vtk *.dat core.* 
	rm -f $(OBJS) $(DEPS) $(NAME) $(OBJDIR)/main_*.[od]

