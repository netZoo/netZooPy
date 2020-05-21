#!/bin/bash

# pre_inter_dir='data/MotifPipeline/MeArray'
# pre_inter_dir='data/MotifPipeline/val_mili_plusNONcpg/redo' ##if you fuck up, pull back down from caleb, parse with last tile in milipeed_tutorial.ipynb
pre_inter_dir='data/MotifPipeline/MeArrayintersectFULLshuf/test'
prefiles=$(ls $pre_inter_dir/*)

mearray='data/MotifPipeline/ENCODE/methyl_array/crossQC'
# mefiles=$(ls $mearray/*_MeArray.txt)
# shufarray='data/MotifPipeline/ENCODE/methyl_array'

outdir='data/MotifPipeline/MeArrayintersectFULLshuf/result/'
# rm -r -i -f $outdir
# mkdir $outdir
counter=1 


for tfile in $prefiles
do
	sed -i '' -e $'s/\t\t/\tNaN\t/g' $tfile

	gtag=$(basename "$tfile" |cut -d _ -f1)
	bbname=$(basename "$tfile" .txt)
	


	cut -f4 $tfile|shuf > tmp3.txt
	cut -f5 $tfile|shuf >tmp4.txt
	paste $tfile tmp3.txt tmp4.txt | awk '{print $1,$2,$3,$4,$5,$6,$7,$14,$15}' OFS='\t' > tmp2.txt
	
	eval "~/bedtools2/bin/bedtools intersect -wa -wb -a tmp2.txt -b $mearray/$gtag"_"shufMeArrayHG38b.txt" > tmp33.txt
	eval "~/bedtools2/bin/bedtools intersect -v -a tmp2.txt -b $mearray/$gtag"_"shufMeArrayHG38b.txt" > tmp.txt

##fix 
	# cut -f1,2,3,6,8,12,13 $tfile> temppp.txt
	# cut -f4 temppp.txt|shuf > tmp3.txt
	# cut -f5 temppp.txt|shuf >tmp4.txt
	# paste temppp.txt tmp3.txt tmp4.txt | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS='\t' > tmp2.txt

	# eval "~/bedtools2/bin/bedtools intersect -wa -wb -a tmp2.txt -b $mearray/$gtag"_"shufMeArrayHG38QC.txt" > $outdir$bbname 
	# eval "~/bedtools2/bin/bedtools intersect -v -a tmp2.txt -b $mearray/$gtag"_"shufMeArrayHG38QC.txt" > tmp.txt
	
	awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,0,0,0,$4,$4}' OFS='\t' tmp.txt >> tmp33.txt 
	awk '{print $1,$2,$3,$4,$8,$5,$9,$13,$14,$6,$7}' OFS='\t' tmp33.txt > $outdir$bbname 

	rm -rf tmp2.txt
	rm -rf tmp3.txt
	rm -rf tmp33.txt
	rm -rf tmp4.txt
	rm -rf tmp.txt
	rm -rf temppp.txt
	rm -rf $tfile

done
