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
motifdir='/udd/redmo/data/MotifPipeline/hg38_refseq_100kb_trFROMhg19/'  ## chr start stop pwm per gene name/
# motiffiles=$(ls $motifdir*.txt)

# WGBSdir=$3
WGBSdir='/udd/redmo/data/MotifPipeline/ENCODE/wgbsin' ## chr start end tmp Me-value, multiple per cell line
WGBSfiles=$(ls $WGBSdir/*.bed)
WGBSmeta=$WGBSdir/cellline_meta.txt


# outdir=$4
outdir='/udd/redmo/data/MotifPipeline/validate_milipeed2/'
rm -r -i -f $outdir
mkdir $outdir
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
awk 'NR > 1 || $4 == $tf' $tfile > tffile.txt ## cell-line and gene specific chip file

gtag=$(eval "echo "$tfile" | cut -d / -f7| cut -d . -f1")

###merge multiple wgbs for same cell line (replicates)
awk -v pat=$gtag '$2 ~ pat'  $WGBSmeta | cut -f 1 | sed 's/^/wgbsin\//' | sed 's/$/.bed/'> tempTF0.txt
# tempTF1=$(cat tempTF0.txt)
# cat `echo $tempTF1` | cut -f1,2,3,5| sort -V -k1,1 -k2,2n > mergeall.txt
# eval "~/../rekrg/Tools/bedtools2/bin/bedtools merge -i mergeall.txt" > tfile.txt #wgbs file

###use largest of replicates (faster)
mergeall=$(find $(< tempTF0.txt) | sort -nr | head -n 1 )

# cut -f1,2,3,5 $mergall > wfile.txt

echo "$counter : TF=$tf"
counter=$[$counter +1]
### restrict WGBS to only motif regions, return both PWM and %Me
eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a $motifdir$tftag -b $mergeall" > temp0a.txt
### restrict that intersection above with hits on WGBS, and if no ChIP peak return zero
eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wao -a temp0a.txt -b tffile.txt" > $outdir$gtag$tf

rm -i -f -r temp0a.txt
rm -i -f -r tempTF0.txt 
rm -i -f -r tffile.txt
rm -i -f -r sort0bed.txt
# 

fi
done
done
# }
