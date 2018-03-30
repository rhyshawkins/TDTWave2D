
BASE=../../TDTbase

INCLUDES = \
	$(shell mpicc -showme:compile) \
	-I$(BASE)/log \
	-I$(BASE)/hnk \
	-I$(BASE)/tracking \
	-I$(BASE)/oset \
	-I$(BASE)/wavelet \
	-I$(BASE)/wavetree

EXTRA_LIBS = \
	-L$(BASE)/tracking -ltracking \
	-L$(BASE)/log -llog


CXX = g++
CXXFLAGS = -c -g -Wall --std=c++11 $(INCLUDES)

CXXFLAGS += -O3

OBJS = birth.o \
	constants.o \
	death.o \
	global.o \
	hierarchical.o \
	hierarchicalprior.o \
	ptexchange.o \
	resample.o \
	rng.o \
	tdtwave2dexception.o \
	tdtwave2dutil.o \
	value.o

all : $(OBJS)


%.o : %.cpp
	$(CXX) $(CXXFLAGS) -o $*.o $*.cpp
