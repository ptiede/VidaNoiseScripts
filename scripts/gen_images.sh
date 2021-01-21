#!/bin/bash


if [ $# != 3 ];
then
  echo "ERROR: requires three arguments"
  echo "1. Chain file"
  echo "2. Number of samples"
  echo "3. Output directory"
  exit 1
fi


chain_file=$1
number_of_images=$2
outdir=$3

nlines=`wc -l $chain_file | awk '{print $1}'`
nlines=`expr $nlines - 2`
simg=`expr $nlines / 2`
echo "This chain file has $nlines"
dn=`expr $nlines / $number_of_images`
dn=`expr $dn / 2`
echo "Will sample $number_of_images images with a spacing of $dn after cutting first half of chain"

head -n 2 $chain_file > $outdir/subchain.dat
tail -n $simg $chain_file > $outdir/tmp.dat
awk -v NUM=$dn 'NR % NUM == 0' $outdir/tmp.dat >> $outdir/subchain.dat
rm $outdir/tmp.dat
nlines=`wc -l $outdir/subchain.dat | awk '{print $1}'`
nlines=`expr $nlines - 2`

for ((i=0; i<$nlines; i+=1)) do
  themispy_image_fit -c $outdir/subchain.dat -N 64 --limits 80 -i $i -o $outdir/themage_$i.fits -r
done

#Now create image list
readlink -f $outdir/*.fits > $outdir/imagelist

