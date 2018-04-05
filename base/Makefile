
BASE=../../TDTbase

INCLUDES = \
	$(shell mpicc -showme:compile) \
	-I$(BASE)/log \
	-I$(BASE)/hnk \
	-I$(BASE)/tracking \
	-I$(BASE)/oset \
	-I$(BASE)/wavelet \
	-I$(BASE)/wavetree

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


CXX = mpicxx
CXXFLAGS = -c -g -Wall --std=c++11 $(INCLUDES)

#CXXFLAGS += -O3

OBJS = birth.o \
	constants.o \
	death.o \
	global.o \
	hierarchical.o \
	hierarchicalmodel.o \
	hierarchicalprior.o \
	ptexchange.o \
	resample.o \
	rng.o \
	tdtwave2dexception.o \
	tdtwave2dutil.o \
	value.o


all : $(OBJS) $(POSTPROCESS)


%.o : %.cpp
	$(CXX) $(CXXFLAGS) -o $*.o $*.cpp

clean :
	rm -f $(OBJS)

