#!/bin/bash

WGBSdir='data/MotifPipeline/ENCODE/wgbsin' ## chr start end tmp Me-value, multiple per cell line
WGBSfiles=$(ls $WGBSdir/*.bed)
WGBSmeta=$WGBSdir/cellline_meta.txt
merged='A549both.txt'

# motifdir='data/MotifPipeline/MeArrayintersectFULLshuf/test/'  ## chr start stop pwm per gene name/
motifdir='data/MotifPipeline/MeArrayintersectFULLshuf/result/'  ## chr start stop pwm per gene name/
motiffiles=$(ls $motifdir*)

outdir='data/MotifPipeline/MeArrayintersectFULLshuf/final/'


for tf in $motiffiles # targeted genes in ChIP cell lines
do

	out=$(basename $tf)
	outt="$(cut -d'_' -f2 <<<$out)"
	celll="$(cut -d'_' -f1 <<<$out)"
	# awk -v pat=$celll '$2 ~ pat'  $WGBSmeta | cut -f 1 | sed 's/^/wgbsin\//' |sed 's/^/ENCODE\//' |sed 's/^/MotifPipeline\//' |sed 's/^/data\//' | sed 's/$/.bed/'> tempTF0.txt
	# mergeall=$(find $(< tempTF0.txt) | sort -nr | head -n 1 )


	# if [ -e $WGBSdir$celll ]
	# then
	# eval "~/bedtools2/bin/bedtools intersect -wa -wb -a $tf -b $mergeall" >   temp2.txt
	# eval "~/bedtools2/bin/bedtools intersect -v -a $tf -b $mergeall" >   temp.txt
	eval "~/bedtools2/bin/bedtools intersect -wa -wb -a $tf -b $WGBSdir/$merged" >   temp2.txt
	eval "~/bedtools2/bin/bedtools intersect -v -a $tf -b $WGBSdir/$merged" >   temp.txt
	
	awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,0}' OFS='\t' temp.txt >> temp2.txt 
	# awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$21}' OFS='\t' temp2.txt > $outdir/$out #$motifdir$out
	awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$15,$16,$17,$18}' OFS='\t' temp2.txt > $outdir/$out #$motifdir$out

	rm -rf temp.txt
	rm -rf temp2.txt
	rm -rf $tf
	rm -rf $WGBSdir/$merged
	# fi
done

