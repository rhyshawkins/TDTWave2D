
CXX=mpicxx
CXXFLAGS=-c -g -Wall -I../base

CXXFLAGS+=-O3

BASEOBJS = \
	../base/tdtwave2dinvert \
	../base/generalvoronois2exception.o \
	../base/hierarchical_model.o \
	../base/pathutil.o \
	../base/prior.o \
	../base/rng.o \
	../base/global.o

LIBS = -L../libraries/log -llog \
	$(shell gsl-config --libs)

TARGETOBJS = genericregression.o

TARGETS = mksyntheticregression \
	regressiontdtwave2d \
	regressiontdtwave2dpt

all : $(TARGETS)

genericregression.o : genericregression.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

../base/%.o : ../base/%.cpp
	make -C ../base  $*.o

mksyntheticregression : $(TARGETOBJS) $(BASEOBJS) ../mksyntheticobservations.o
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../mksyntheticobservations.o $(LIBS)

regressiontdtwave2d : $(TARGETOBJS) $(BASEOBJS) ../base/tdtwave2dinvert.o 
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../base/tdtwave2dinvert .o $(LIBS)

regressiontdtwave2dpt : $(TARGETOBJS) $(BASEOBJS) ../base/tdtwave2dinvertpt.o 
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../base/tdtwave2dinvertpt.o $(LIBS)

clean :
	rm -f $(TARGETS)
