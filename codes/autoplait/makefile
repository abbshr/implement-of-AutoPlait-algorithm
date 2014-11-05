CFLAGS= -g
INCS=
CC=gcc
CODES = autoplait 
all :  $(CODES)
#------------------------------------------------------#
autoplait: autoplait.o plait.o cps.o tool.o nrutil.o hmm.h autoplait.h \
	viterbi.o hmmutils.o \
	baum.o forbackward.o kmeans.o segbox.o
	$(CC) -o autoplait \
	autoplait.o plait.o cps.o tool.o nrutil.o \
	viterbi.o hmmutils.o \
	baum.o forbackward.o kmeans.o segbox.o \
	-lm
#------------------------------------------------------#

#make plot: plots_all.m
#	matlab -r 'plots_all' 
#make demo: demo.sh
#	sh demo.sh

cleanall: 
	rm -rf core *.o *~ $(CODES) ./_out/dat_tmp/
clean :
	rm -f core *.o *~

