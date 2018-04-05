
TARGETS = datatemplate.txt \
	imagetemplateconstant.txt \
	imagetemplategaussian.txt \
	syntheticobs_constant.txt \
	syntheticobs_gaussian.txt \
	results_constant/ch.dat \
	results_gaussian/ch.dat

all : $(TARGETS)

datatemplate.txt :
	python2 generatetemplatepoints.py -N 100 -o $@

imagetemplateconstant.txt :
	python2 generatetemplateimage.py -W 32 -H 32 -o $@ -A 0.0

imagetemplategaussian.txt :
	python2 generatetemplateimage.py -W 32 -H 32 -o $@ -m Gaussian

syntheticobs_constant.txt : datatemplate.txt imagetemplateconstant.txt
	../mksyntheticregression -i imagetemplateconstant.txt -I datatemplate.txt \
	-o $@ -N 0.05

syntheticobs_gaussian.txt : datatemplate.txt imagetemplategaussian.txt
	../mksyntheticregression -i imagetemplategaussian.txt -I datatemplate.txt \
	-o $@ -N 0.05

results_constant/ch.dat : syntheticobs_constant.txt
	mkdir -p results_constant
	../regressiontdtwave2d -i syntheticobs_constant.txt \
	-M priorproposal.txt \
	-o results_constant/ \
	-x 5 -y 5 -t 100000 -w 4

results_constant/mean.txt : results_constant/ch.dat
	../postprocess_mean -w 4 -x 5 -y 5 -i $< -o $@ -t 5 -s 50000

results_gaussian/ch.dat : syntheticobs_gaussian.txt
	mkdir -p results_gaussian
	../regressiontdtwave2d -i syntheticobs_gaussian.txt \
	-M priorproposal.txt \
	-o results_gaussian/ \
	-v 10000 \
	-x 5 -y 5 -t 1000000 -w 4

results_gaussian/mean.txt : results_gaussian/ch.dat
	../postprocess_mean -w 4 -x 5 -y 5 -i $< -o $@ -t 5 -s 500000
