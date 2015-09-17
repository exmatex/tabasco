SRCDIR=src

##############CHARM##########
ifeq (,$(CHARM_ROOT))
$(info Please establish Charmm++ environment variables before using this Makefile.)
$(info E.g. by running setting CHARM_ROOT pr run 'module load charm++')
$(error CHARM_ROOT is not set)
endif
CHARMINC=$(CHARM_ROOT)/include/
CHARMBIN=$(CHARM_ROOT)/bin
CXX=$(CHARMBIN)/charmc
OBJDIR=charm_obj
BINDIR=charm_bin

#OPTFLAGS=-g
OPTFLAGS= -g -O3 -I../CoEVP -I../CoEVP/CM/include -I../CoEVP/CM/src -I../flann/src/cpp
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

OBJS:=$(addprefix $(OBJDIR)/, CoarseScaleModel.o FineScaleModel.o NearestNeighborSearch.o Interpolate.o DBInterface.o input.o Main.o)
# ../CoEVP/LULESH/lulesh.o)
CXXFLAGS+=$(OPTFLAGS) $(BOOST_CFLAG) $(DEFS)
LDFLAGS= -L../CoEVP/CM/lib -lcm $(LINALG_LDFLAG) -lm -lrt 

LULESH_OBJ=../CoEVP/LULESH/lulesh.o

#target
NAME=$(BINDIR)/TabaSCo
default: all
all: $(SUBDIRS) $(OBJDIR) $(NAME)

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
$(OBJDIR)/%.o: $(SRCDIR)/%.cc
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

