
TARGETS = datatemplate.txt \
	imagetemplate.txt \
	syntheticobs.txt

all : $(TARGETS)

datatemplate.txt :
	python2 generatetemplatepoints.py -N 100 -o $@

imagetemplate.txt :
	python2 generatetemplateimage.py -W 32 -H 32 -o $@

syntheticobs.txt :
	../mksyntheticregression -i imagetemplate.txt -I datatemplate.txt \
	-o $@ -N 1.0
