
CXX=mpicxx
CXXFLAGS=-c -g -Wall -I../base

CXXFLAGS+=-O3

BASEOBJS = \
	../base/birth.o \
	../base/constants.o \
	../base/death.o \
	../base/global.o \
	../base/hierarchical.o \
	../base/hierarchicalprior.o \
	../base/ptexchange.o \
	../base/resample.o \
	../base/rng.o \
	../base/tdtwave2dexception.o \
	../base/tdtwave2dimage.o \
	../base/tdtwave2dutil.o \
	../base/value.o

BASELIB=../../TDTbase
LIBS = -L$(BASELIB)/hnk -lhnk \
	-L$(BASELIB)/tracking -ltracking \
	-L$(BASELIB)/log -llog \
	-L$(BASELIB)/wavetree -lwavetree \
	-L$(BASELIB)/sphericalwavelet -lsphericalwavelet \
	-L$(BASELIB)/oset -loset \
	-L$(BASELIB)/wavelet -lwavelet \
	$(shell gsl-config --libs) \
	-lgmp

TARGETOBJS = genericregression.o

TARGETS = mksyntheticregression \
	regressiontdtwave2d \
	regressiontdtwave2dpt

all : $(TARGETS)

genericregression.o : genericregression.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

../base/%.o : ../base/%.cpp
	make -C ../base  $*.o

mksyntheticregression : $(TARGETOBJS) $(BASEOBJS) ../base/mksyntheticobservations.o
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../base/mksyntheticobservations.o $(LIBS)

regressiontdtwave2d : $(TARGETOBJS) $(BASEOBJS) ../base/tdtwave2dinvert.o 
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../base/tdtwave2dinvert.o $(LIBS)

regressiontdtwave2dpt : $(TARGETOBJS) $(BASEOBJS) ../base/tdtwave2dinvert_pt.o 
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../base/tdtwave2dinvert_pt.o $(LIBS)

clean :
	rm -f $(TARGETS)
