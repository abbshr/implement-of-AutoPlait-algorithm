#!/bin/sh
make cleanall
make
INPUTDIR="./_dat/"
OUTDIR="./_out/"

#----------------------#
echo "----------------------"
echo "mocap and googleTrend"
echo "----------------------"
outdir=$OUTDIR"dat_tmp"
dblist=$INPUTDIR"list"
n=3  # data size
d=4  # dimension
#----------------------#

mkdir $outdir
for (( i=1; i<=$n; i++ ))
do
  output=$outdir"/dat"$i"/" 
  mkdir $output
  input=$output"input"
  awk '{if(NR=='$i') print $0}' $dblist > $input
  ./autoplait $d $input $output 
done




