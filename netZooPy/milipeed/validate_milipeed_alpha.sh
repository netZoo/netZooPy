#!/bin/bash
###BOTH STEPS IN AUROC CALC btwn PWM + wgbsmeth + CHiP meth
## for all ChIP files associated with Motif which will me intersected with WGBS/Me array
#PARAMETER REGION
# validate_milipeed ChIPdir motifdir WGBSdir valoutdir {
# ChIPdir=$1 # %%% location of ChIP bed files
ChIPdir='/udd/redmo/data/MotifPipeline/remap/' ## chr start end gene-target per cell line name
# tfdb=$ChIPdir/meta2IDR.txt # metadata file including ChIP file name and gene
TFfiles=$(ls $ChIPdir*.bed)


# motifdir=$2
motifdir='/udd/redmo/data/MotifPipeline/hg38_refseq_100kb_tr_fromhg19/'  ## chr start stop pwm per gene name/
# motiffiles=$(ls $motifdir*.txt)

# WGBSdir=$3
WGBSdir='/udd/redmo/data/MotifPipeline/ENCODE/wgbsin' ## chr start end tmp Me-value, multiple per cell line
WGBSfiles=$(ls $WGBSdir/*.bed)
WGBSmeta=$WGBSdir/cellline_meta.txt


# outdir=$4
outdir='/udd/redmo/data/MotifPipeline/val_mili_plusNONcpg/'
# rm -r -i -f $outdir
# mkdir $outdir
counter=1 
# for wfile in $WGBSfiles
# do
# echo $wfile
for tfile in $TFfiles # cell line
do

TFs=$(eval "cat "$tfile" | cut -f4 | sort | uniq")

for tf in $TFs # targeted genes in ChIP cell lines
do
tftag="${tf}.txt"  #$(eval "echo "$tfile" | cut -d / -f7")

if [ -e $motifdir/$tftag ]
then
if (( $counter >250 ))
then
awk 'NR > 1 || $4 == $tf' $tfile > tffilee.txt ## cell-line and gene specific chip file

gtag=$(eval "echo "$tfile" | cut -d / -f7| cut -d . -f1")

###merge multiple wgbs for same cell line (replicates)
awk -v pat=$gtag '$2 ~ pat'  $WGBSmeta | cut -f 1 | sed 's/^/wgbsin\//' | sed 's/$/.bed/'> tempTF00.txt
# tempTF1=$(cat tempTF00.txt)
# cat `echo $tempTF1` | cut -f1,2,3,5| sort -V -k1,1 -k2,2n > mergeall.txt
# eval "~/../rekrg/Tools/bedtools2/bin/bedtools merge -i mergeall.txt" > tfile.txt #wgbs file

###use largest of replicates (faster)
mergeall=$(find $(< tempTF00.txt) | sort -nr | head -n 1 )

# cut -f1,2,3,5 $mergall > wfile.txt

echo "$counter : TF=$tf"
counter=$[$counter +1]
### restrict WGBS to only motif regions, return both PWM and %Me
eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a $motifdir$tftag -b $mergeall" > temp0aa.txt
cut -f1,2,3,4,17 temp0aa.txt >temp0bb.txt
# min=$(eval "cut -f4 temp0bb.txt | scale=4 | bc ")
# min=$(eval "cut -f4 temp0bb.txt | sort -n | head -1")
min=$(cat temp0bb.txt | cut -f4 |sort -n | head -1)
max=$(cat temp0bb.txt | cut -f4 |sort -n | tail -1)
awk '{print $1,$2,$3,$4,$5,($4-m)/(ma-m),1-$5/100,1-$5/100}' m="$min" ma="$max" OFS='\t' temp0bb.txt > temp0dd.txt ##standardize range of PWM and Methyl

eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -v -a $motifdir$tftag -b $mergeall" > temp0cc.txt
min=$(cat temp0cc.txt | cut -f4 |sort -n | head -1)
max=$(cat temp0cc.txt | cut -f4 |sort -n | tail -1)
awk '{print $1,$2,$3,$4,$5,($4-m)/(ma-m),($4-m)/(ma-m),1}' m="$min" ma="$max" OFS='\t' temp0cc.txt >> temp0dd.txt ##combine motif with and without CpG

# cat temp0cc.txt |awk 'BEGIN {FS="\t"}; {print $1 $2 $3 $4 $4}' > temp0dd.txt
### restrict that intersection above with hits on WGBS, and if no ChIP peak return zero
eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wao -a temp0dd.txt -b tffilee.txt" > $outdir${gtag}_${tf} ##compare entire motif with new methyaltion weights inserted

rm -i -f -r temp0aa.txt
rm -i -f -r temp0bb.txt
rm -i -f -r temp0cc.txt
rm -i -f -r temp0dd.txt
rm -i -f -r tempTF00.txt 
rm -i -f -r tffilee.txt
# rm -i -f -r sort0bed.txt
else
counter=$[$counter +1]
# break
fi
fi
done
done
# }
