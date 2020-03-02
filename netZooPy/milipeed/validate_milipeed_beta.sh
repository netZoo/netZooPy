#!/bin/bash
###BOTH STEPS IN AUROC CALC btwn PWM + wgbsmeth + CHiP meth
## for all ChIP files associated with Motif which will me intersected with WGBS/Me array
#PARAMETER REGION
validate_milipeed TFdir motifdir bsfile valoutdir {
TFdir=$1 # %%% location of ChIP bed files
tfdb=$TFdir/meta2IDR.txt # metadata file including ChIP file name and gene
motifdir=$2
motiffiles=$(ls $motifdir*.txt)
TFfiles=$(ls $TFdir/*.bed)
bsfile=$3
outdir=$4
rm -r -i -f $outdir
mkdir $outdir

# while read mfile #
for mfile in $motiffiles
do

if [ -e $mfile ]
then

mtag=$(eval "echo "$mfile" | cut -d / -f3|cut -d . -f1")
bstag=$(eval "echo "$bsfile" | cut -d . -f1| cut -d / -f2")

awk -v pat=$mtag '$2 ~ pat'  $tfdb | cut -f 1 | sed 's/^/'TFdir'\//'| sed 's/$/.bed/'> tempTF0.txt

if [ -s "tempTF0.txt" ]
then

tempTF1=$(cat tempTF0.txt)
cat `echo $tempTF1` | cut -f1,2,3,5| sort -V -k1,1 -k2,2n > tffile.txt
echo "$count = wgbs=$bstag gene=$gtag TF=$mtag"
sort -k1,1 -k2,2n $mfile > sort_motif.txt

### restrict WGBS to only motif regions, return both PWM and %Me
eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a sort_motif.txt -b wgbsin/ENCFF005TID.txt " > temp0a.txt
### restrict that intersection above with hits on WGBS, and if no ChIP peak return zero
eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wao -a temp0a.txt -b tffile.txt" > $outdir$mtag$ttag 

rm -i -f -r temp0a.txt
rm -i -f -r tempTF0.txt 
rm -i -f -r tffile.txt
rm -i -f -r sort0bed.txt

fi
fi
done
}
