#!/bin/bash

# pre_inter_dir='data/MotifPipeline/MeArray'
pre_inter_dir='data/MotifPipeline/val_mili_plusNONcpg/redo' ##if you fuck up, pull back down from caleb, parse with last tile in milipeed_tutorial.ipynb

prefiles=$(ls $pre_inter_dir/*)

mearray='data/MotifPipeline/ENCODE/methyl_array'
# mefiles=$(ls $mearray/*_MeArray.txt)

outdir='data/MotifPipeline/MeArrayintersect/'
rm -r -i -f $outdir
mkdir $outdir
counter=1 


for tfile in $prefiles
do
	sed -i -e $'s/\t\t/\tNaN\t/g' $tfile
	gtag=$(basename "$tfile" |cut -d _ -f1)
	bbname=$(basename "$tfile" .txt)
	eval "~/bedtools2/bin/bedtools intersect -wa -wb -a $tfile -b $mearray/$gtag"_"MeArrayHG38.txt" > $outdir$bbname 
	# eval "~/bedtools2/bin/bedtools intersect -v -a $tfile -b $mearray/$gtag"_"MeArray.txt" > tmp.txt
	# awk '{print $1,$2,$3,$4,$5,$6,$7,0,0,0,$4}' OFS='\t' tmp.txt >> $outdir$bbname 
done
