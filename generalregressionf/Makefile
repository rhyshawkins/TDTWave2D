
CXX=mpicxx
CXXFLAGS=-c -g -Wall -I../base

FCC=gfortran
FCCFLAGS=-c -g -Wall

#CXXFLAGS+=-O3

BASEOBJS = \
	../base/birth.o \
	../base/constants.o \
	../base/death.o \
	../base/global.o \
	../base/hierarchical.o \
	../base/hierarchicalmodel.o \
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
	-lgmp \
	-lgfortran

TARGETOBJS = genericregression.o

TARGETS = mksyntheticregression \
	regressiontdtwave2d \
	regressiontdtwave2dpt

POSTPROCESS = postprocess_khistory \
	postprocess_likelihood \
	postprocess_mean \
	postprocess_mean_mpi



all : $(TARGETS) $(POSTPROCESS)

genericregression.o : genericregression.f90
	$(FCC) $(FCCFLAGS) -o $@ $<

../base/%.o : ../base/%.cpp
	make -C ../base  $*.o

mksyntheticregression : $(TARGETOBJS) $(BASEOBJS) ../base/mksyntheticobservations.o
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../base/mksyntheticobservations.o $(LIBS)

regressiontdtwave2d : $(TARGETOBJS) $(BASEOBJS) ../base/tdtwave2dinvert.o 
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../base/tdtwave2dinvert.o $(LIBS)

regressiontdtwave2dpt : $(TARGETOBJS) $(BASEOBJS) ../base/tdtwave2dinvert_pt.o 
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../base/tdtwave2dinvert_pt.o $(LIBS)

postprocess_khistory : $(TARGETOBJS) $(BASEOBJS) ../base/postprocess_khistory.o 
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../base/postprocess_khistory.o $(LIBS)

postprocess_likelihood : $(TARGETOBJS) $(BASEOBJS) ../base/postprocess_likelihood.o 
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../base/postprocess_likelihood.o $(LIBS)

postprocess_mean : $(TARGETOBJS) $(BASEOBJS) ../base/postprocess_mean.o 
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../base/postprocess_mean.o $(LIBS)

postprocess_mean_mpi : $(TARGETOBJS) $(BASEOBJS) ../base/postprocess_mean_mpi.o 
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../base/postprocess_mean_mpi.o $(LIBS)

clean :
	rm -f $(TARGETOBJS) $(TARGETS)
