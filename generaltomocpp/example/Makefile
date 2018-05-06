
TARGETS = datatemplate.txt \
	imagetemplateconstant.txt \
	imagetemplategaussian.txt \
	syntheticobs_constant.txt \
	syntheticobs_gaussian.txt \
	results_constant/ch.dat \
	results_gaussian/ch.dat

all : $(TARGETS)

datatemplate.txt :
	python2 ../../scripts/generatetemplatepaths.py -N 100 --samples 2 -o $@ \
	--minx -10 --maxx 10 \
	--miny -10 --maxy 10

imagetemplateconstant.txt :
	python2 ../../scripts/generatetemplateimage.py -W 32 -H 32 -o $@ -A 3.0

imagetemplategaussian.txt :
	python2 ../../scripts/generatetemplateimage.py -W 32 -H 32 -o $@ -m Gaussian -A 2.0 -O 2.5

syntheticobs_constant.txt : datatemplate.txt imagetemplateconstant.txt
	../mksynthetictomography -i imagetemplateconstant.txt -I datatemplate.txt \
	-o $@ -N 5.0

syntheticobs_gaussian.txt : datatemplate.txt imagetemplategaussian.txt
	../mksynthetictomography -i imagetemplategaussian.txt -I datatemplate.txt \
	-o $@ -N 5.0

#
# Serial Constant True model
#
results_constant/ch.dat : syntheticobs_constant.txt
	mkdir -p results_constant
	../tomographytdtwave2d -i syntheticobs_constant.txt \
	-M priorproposal.txt \
	-o results_constant/ \
	-x 5 -y 5 -t 100000 -w 4

results_constant/mean.txt : results_constant/ch.dat
	../postprocess_mean -w 4 -x 5 -y 5 -i $< -o $@ -t 5 -s 50000

#
# Serial Gaussian Hill True model
#
results_gaussian/ch.dat : syntheticobs_gaussian.txt
	mkdir -p results_gaussian
	../tomographytdtwave2d -i syntheticobs_gaussian.txt \
	-M priorproposal.txt \
	-o results_gaussian/ \
	-v 10000 \
	-x 5 -y 5 -t 1000000 -w 4

results_gaussian/mean.txt : results_gaussian/ch.dat
	../postprocess_mean -w 4 -x 5 -y 5 -i $< -o $@ -t 5 -s 500000

#
# Serial Hierarchical Prior Constant
#
results_constant_hp/ch.dat : syntheticobs_constant.txt
	mkdir -p results_constant_hp
	../tomographytdtwave2d -i syntheticobs_constant.txt \
	-M priorproposal.txt \
	-o results_constant_hp/ \
	-x 5 -y 5 -t 100000 -w 4 \
	-p 0.1

results_constant_hp/mean.txt : results_constant_hp/ch.dat
	../postprocess_mean -w 4 -x 5 -y 5 -i $< -o $@ -t 5 -s 50000

#
# Serial Hierarchical Prior Gaussian
#
results_gaussian_hp/ch.dat : syntheticobs_gaussian.txt
	mkdir -p results_gaussian_hp
	../tomographytdtwave2d -i syntheticobs_gaussian.txt \
	-M priorproposal.txt \
	-o results_gaussian_hp/ \
	-x 5 -y 5 -t 1000000 -w 4 \
	-p 0.1

results_gaussian_hp/mean.txt : results_gaussian_hp/ch.dat
	../postprocess_mean -w 4 -x 5 -y 5 -i $< -o $@ -t 5 -s 500000

#
# Serial Hierarchical Noise Constant
#
results_constant_hn/ch.dat : syntheticobs_constant.txt
	mkdir -p results_constant_hn
	../tomographytdtwave2d -i syntheticobs_constant.txt \
	-M priorproposal.txt \
	-o results_constant_hn/ \
	-x 5 -y 5 -t 100000 -w 4 \
	-L 0.01

results_constant_hn/mean.txt : results_constant_hn/ch.dat
	../postprocess_mean -w 4 -x 5 -y 5 -i $< -o $@ -t 5 -s 50000

#
# Serial Hierarchical Noise Gaussian
#
results_gaussian_hn/ch.dat : syntheticobs_gaussian.txt
	mkdir -p results_gaussian_hn
	../tomographytdtwave2d -i syntheticobs_gaussian.txt \
	-M priorproposal.txt \
	-o results_gaussian_hn/ \
	-x 5 -y 5 -t 1000000 -w 4 \
	-v 10000 \
	-L 0.01

results_gaussian_hn/mean.txt : results_gaussian_hn/ch.dat
	../postprocess_mean -w 4 -x 5 -y 5 -i $< -o $@ -t 5 -s 500000


#
# Parallel Constant
#
results_constant_mpi/ch.dat-000 : syntheticobs_constant.txt
	mkdir -p results_constant_mpi
	mpirun -np 4 ../tomographytdtwave2dpt -c 4 -T 1 \
	-i syntheticobs_constant.txt \
	-M priorproposal.txt \
	-o results_constant_mpi/ \
	-x 5 -y 5 -t 100000 -w 4 

results_constant_mpi/mean.txt : results_constant_mpi/ch.dat-000
	mpirun -np 4 ../postprocess_mean_mpi -i results_constant_mpi/ch.dat -o $@ \
	-w 4 -x 5 -y 5 -t 5 -s 50000

#
# Parallel Gaussian
#
results_gaussian_mpi/ch.dat-000 : syntheticobs_gaussian.txt
	mkdir -p results_gaussian_mpi
	mpirun -np 4 ../tomographytdtwave2dpt -c 4 -T 1 \
	-i syntheticobs_gaussian.txt \
	-M priorproposal.txt \
	-o results_gaussian_mpi/ \
	-x 5 -y 5 -t 1000000 -v 10000 -w 4 

results_gaussian_mpi/mean.txt : results_gaussian_mpi/ch.dat-000
	mpirun -np 4 ../postprocess_mean_mpi -i results_gaussian_mpi/ch.dat -o $@ \
	-w 4 -x 5 -y 5 -t 15 -s 1500000
