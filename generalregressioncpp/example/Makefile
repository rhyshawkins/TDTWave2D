
TARGETS = datatemplate.txt \
	imagetemplate.txt \
	syntheticobs.txt \
	results/ch.dat

all : $(TARGETS)

datatemplate.txt :
	python2 generatetemplatepoints.py -N 100 -o $@

imagetemplate.txt :
	python2 generatetemplateimage.py -W 32 -H 32 -o $@

syntheticobs.txt : datatemplate.txt imagetemplate.txt
	../mksyntheticregression -i imagetemplate.txt -I datatemplate.txt \
	-o $@ -N 1.0

results/ch.dat : syntheticobs.txt
	mkdir -p results
	../regressiontdtwave2d -i syntheticobs.txt \
	-M priorproposal.txt \
	-o results/ \
	-x 5 -y 5 -t 100000 -w 4

