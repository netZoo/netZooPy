#!/bin/bash


# WGBSdir=$3
WGBSdir='data/MotifPipeline/test' ## chr start end tmp Me-value, multiple per cell line
WGBSfiles=$(ls $WGBSdir/*)

motifdir='data/MotifPipeline/hg38_refseq_100kb_tr_fromhg19/'  ## chr start stop pwm per gene name/
motiffiles=$(ls $motifdir*.txt)
aaa='.txt'
for tf in $WGBSfiles # targeted genes in ChIP cell lines
do
	out=$(basename $tf)
	outt="$(cut -d'_' -f2 <<<$out)"
if [ -e $motifdir$outt$aaa ]
then
	eval "~/bedtools2/bin/bedtools intersect -wa -wb -a $tf -b $motifdir$outt$aaa" >   temp.txt
	awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$13,$14,$18,$19}' OFS='\t' temp.txt > $WGBSdir/$out #$motifdir$out

	rm temp.txt
fi
done

