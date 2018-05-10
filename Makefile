
default : all

all :
	make -C base
	make -C generalregressioncpp
	make -C generaltomocpp
	make -C generalregressionf
	make -C generaltomof

BASESRCS = \
	base/birth.cpp \
	base/birth.hpp \
	base/constants.cpp \
	base/constants.hpp \
	base/death.cpp \
	base/death.hpp \
	base/genericinterface.hpp \
	base/global.cpp \
	base/global.hpp \
	base/hash.hpp \
	base/hierarchical.cpp \
	base/hierarchical.hpp \
	base/hierarchicalmodel.cpp \
	base/hierarchicalmodel.hpp \
	base/hierarchicalprior.cpp \
	base/hierarchicalprior.hpp \
	base/mksyntheticobservations.cpp \
	base/postprocess_khistory.cpp \
	base/postprocess_likelihood.cpp \
	base/postprocess_mean.cpp \
	base/postprocess_mean_mpi.cpp \
	base/ptexchange.cpp \
	base/ptexchange.hpp \
	base/resample.cpp \
	base/resample.hpp \
	base/rng.cpp \
	base/rng.hpp \
	base/tdtwave2dexception.cpp \
	base/tdtwave2dexception.hpp \
	base/tdtwave2dimage.cpp \
	base/tdtwave2dimage.hpp \
	base/tdtwave2dinvert.cpp \
	base/tdtwave2dinvert_pt.cpp \
	base/tdtwave2dutil.cpp \
	base/tdtwave2dutil.hpp \
	base/value.cpp \
	base/value.hpp

DOC = documentation\manual.tex \
	documentation\manual.pdf \
	documentation\bibliography.bib

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

GENERALREGRESSIONFSRCS = \
	generalregressionf/Makefile \
	generalregressionf/genericregression.f90 \
	generalregressionf/example/Makefile \
	generalregressionf/example/priorproposal.txt

GENERALTOMOGRAPHYFSRCS = \
	generaltomof/Makefile \
	generaltomof/generictomography.f90 \
	generaltomof/linearweights.f90 \
	generaltomof/example/Makefile \
	generaltomof/example/priorproposal.txt

SRCS = $(BASESRCS) \
	$(DOCUMENTATION) \
	$(GENERALREGRESSIONCPPSRCS) \
	$(GENERALTOMOGRAPHYCPPSRCS) \
	$(GENERALREGRESSIONFSRCS) \
	$(GENERALTOMOGRAPHYFSRCS)

documentation/manual.pdf : documentation/manual.tex
	cd documentation && pdflatex manual.tex && bibtex manual && pdflatex manual.tex && pdflatex manual.tex 

INSTALL = install
INSTALLFLAGS = -D
DIR = TDTWave2D
DATE = $(shell date +"%Y%m%d%H%M")
TGZ = $(DIR).tar.gz

dist :
	mkdir -p $(DIR)
	echo $(DATE) > $(DIR)/Version
	for f in $(SRCS); do \
	    $(INSTALL) $(INSTALLFLAGS) $$f $(DIR)/$$f ; \
	done
	tar -czf $(TGZ) $(DIR)/*
	rm -rf $(DIR)

