
default : all

all :
	make -C base
	make -C generalregressioncpp
	make -C generaltomocpp
	make -C generalregressionf

BASESRCS = \
	base/birth.cpp \
	base/birth.hpp \
	base/constants.cpp \
	base/constants.hpp \
	base/death.cpp \
	base/death.hpp \
	base/genericinterface.hpp \
	global.cpp \
	global.hpp \
	hash.hpp \
	hierarchical.cpp \
	hierarchical.hpp \
	hierarchicalmodel.cpp \
	hierarchicalmodel.hpp \
	hierarchicalprior.cpp \
	hierarchicalprior.hpp \
	mksyntheticobservations.cpp \
	postprocess_khistory.cpp \
	postprocess_likelihood.cpp \
	postprocess_mean.cpp \
	postprocess_mean_mpi.cpp \
	ptexchange.cpp \
	ptexchange.hpp \
	resample.cpp \
	resample.hpp \
	rng.cpp \
	rng.hpp \
	tdtwave2dexception.cpp \
	tdtwave2dexception.hpp \
	tdtwave2dimage.cpp \
	tdtwave2dimage.hpp \
	tdtwave2dinvert.cpp \
	tdtwave2dinvert_pt.cpp \
	tdtwave2dutil.cpp \
	tdtwave2dutil.hpp \
	value.cpp \
	value.hpp

GENERALREGRESSIONCPPSRCS = \
	generalregressioncpp/Makefile \
	generalregressioncpp/genericregression.cpp \
	generalregressioncpp/example/Makefile \
	generalregressioncpp/example/priorproposal.txt

GENERALTOMOGRAPHYCPPSRCS = \
	generaltomocpp/Makefile \
	generaltomocpp/generictomography.cpp \
	generaltomocpp/linearweights.hpp \
	generaltomocpp/example/Makefile \
	generaltomocpp/example/priorproposal.txt

SRCS = $(BASESRCS)
