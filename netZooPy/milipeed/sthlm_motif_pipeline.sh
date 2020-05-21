#make single motif file
# cat data/MotifPipeline/MeArrayintersectFULLshuf/final/* >> all.txt 
# cut -f1,2,3 all.txt|uniq >all123.txt  

##combine wgbs duplicates
# paste data/MotifPipeline/ENCODE/wgbsin/ENCFF003JVR.txt data/MotifPipeline/ENCODE/wgbsin/ENCFF005TID.txt |awk '{print $1,$2,$3,$4,$5,$9,$10}' OFS='\t' > data/MotifPipeline/ENCODE/wgbsin/A549both.txt
# ~/bedtools2/bin/bedtools intersect -wa -wb -a data/MotifPipeline/ENCODE/wgbsin/A549both.txt -b data/MotifPipeline/ENCODE/methyl_array/crossQC/A-549_shufMeArrayHG38b.txt data/MotifPipeline/remap/A-549.bed -names me chip > A549_overlap_NOmotif.txt
# ~/bedtools2/bin/bedtools intersect -wao -a data/MotifPipeline/ENCODE/wgbsin/A549both.txt -b data/MotifPipeline/ENCODE/methyl_array/crossQC/A-549_shufMeArrayHG38b.txt data/MotifPipeline/remap/A-549.bed >> A549_overlap_NOmotif.txt

# paste data/MotifPipeline/ENCODE/wgbsin/ENCFF003JVR.bed data/MotifPipeline/ENCODE/wgbsin/ENCFF005TID.bed |awk '{print $1,$2,$3,$10,$11}' OFS='\t' > data/MotifPipeline/ENCODE/wgbsin/A549both.txt
# paste data/MotifPipeline/ENCODE/wgbsin/ENCFF804NTQ.bed data/MotifPipeline/ENCODE/wgbsin/ENCFF696OLO.bed |awk '{print $1,$2,$3,$10,$11}' OFS='\t' > data/MotifPipeline/ENCODE/wgbsin/HeLaboth.txt
# paste data/MotifPipeline/ENCODE/wgbsin/ENCFF847OWL.bed data/MotifPipeline/ENCODE/wgbsin/ENCFF390OZB.bed |awk '{print $1,$2,$3,$10,$11}' OFS='\t' > data/MotifPipeline/ENCODE/wgbsin/HepG2both.txt
# paste data/MotifPipeline/ENCODE/wgbsin/ENCFF279HCL.bed data/MotifPipeline/ENCODE/wgbsin/ENCFF835NTC.bed |awk '{print $1,$2,$3,$10,$11}' OFS='\t' > data/MotifPipeline/ENCODE/wgbsin/GM12878both.txt
# paste data/MotifPipeline/ENCODE/wgbsin/ENCFF867JRG.bed data/MotifPipeline/ENCODE/wgbsin/ENCFF721JMB.bed |awk '{print $1,$2,$3,$10,$11}' OFS='\t' > data/MotifPipeline/ENCODE/wgbsin/K562both.txt
# paste data/MotifPipeline/ENCODE/wgbsin/ENCFF940XWW.bed data/MotifPipeline/ENCODE/wgbsin/ENCFF179VKR.bed |awk '{print $1,$2,$3,$10,$11}' OFS='\t' > data/MotifPipeline/ENCODE/wgbsin/SKNSHboth.txt


# eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a data/MotifPipeline/ENCODE/wgbsin/A549both.txt -b data/MotifPipeline/ENCODE/methyl_array/crossQC/A-549_shufMeArrayHG38b.txt"|cut -f1,2,3,4,5,10 > A549_overlap_NOmotif_a.txt
# eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wao -a A549_overlap_NOmotif_a.txt -b data/MotifPipeline/remap/A-549.bed "|cut -f1,2,3,4,5,6,10,11 > A549_overlap_NOmotif_b.txt
# eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -v -a A549_overlap_NOmotif_a.txt -b data/MotifPipeline/remap/A-549.bed ">> A549_overlap_NOmotif_b.txt

###***this second way works better to keep columns similar***###
eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a data/MotifPipeline/ENCODE/wgbsin/A549both.txt -b data/MotifPipeline/ENCODE/methyl_array/crossQC/A-549_shufMeArrayHG38b.txt"|cut -f1,2,3,4,5,9 > A549_overlap_NOmotif_a.txt
eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wao -a A549_overlap_NOmotif_a.txt -b data/MotifPipeline/remap/A-549.bed "|cut -f1,2,3,4,5,6,10,11 > A549_overlap_NOmotif_b.txt

# eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a data/MotifPipeline/ENCODE/wgbsin/HeLaboth.txt -b data/MotifPipeline/ENCODE/methyl_array/crossQC/HeLa-S3_shufMeArrayHG38b.txt data/MotifPipeline/remap/HeLa-S3.bed -names me chip" > HeLa_overlap_NOmotif.txt
# eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wao -a data/MotifPipeline/ENCODE/wgbsin/HeLaboth.txt -b data/MotifPipeline/ENCODE/methyl_array/crossQC/HeLa-S3_shufMeArrayHG38b.txt data/MotifPipeline/remap/HeLa-S3.bed -names wg other" >> HeLa_overlap_NOmotif.txt
eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a data/MotifPipeline/ENCODE/wgbsin/HeLaboth.txt -b data/MotifPipeline/ENCODE/methyl_array/crossQC/HeLa-S3_shufMeArrayHG38b.txt"|cut -f1,2,3,4,5,9 > HeLa_overlap_NOmotif_a.txt
eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wao -a HeLa_overlap_NOmotif_a.txt -b data/MotifPipeline/remap/HeLa-S3.bed "|cut -f1,2,3,4,5,6,10,11 > HeLa_overlap_NOmotif_b.txt

# eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a data/MotifPipeline/ENCODE/wgbsin/HepG2both.txt -b data/MotifPipeline/ENCODE/methyl_array/crossQC/Hep-G2_shufMeArrayHG38b.txt data/MotifPipeline/remap/Hep-G2.bed -names me chip" > HepG2_overlap_NOmotif.txt
# eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wao -a data/MotifPipeline/ENCODE/wgbsin/HepG2both.txt -b data/MotifPipeline/ENCODE/methyl_array/crossQC/Hep-G2_shufMeArrayHG38b.txt data/MotifPipeline/remap/Hep-G2.bed -names wg other" >> HepG2_overlap_NOmotif.txt
eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a data/MotifPipeline/ENCODE/wgbsin/HepG2both.txt -b data/MotifPipeline/ENCODE/methyl_array/crossQC/Hep-G2_shufMeArrayHG38b.txt"|cut -f1,2,3,4,5,9 > HepG2_overlap_NOmotif_a.txt
eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wao -a HepG2_overlap_NOmotif_a.txt -b data/MotifPipeline/remap/Hep-G2.bed "|cut -f1,2,3,4,5,6,10,11 > HepG2_overlap_NOmotif_b.txt

# eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a data/MotifPipeline/ENCODE/wgbsin/GM12878both.txt -b data/MotifPipeline/ENCODE/methyl_array/crossQC/GM12878_shufMeArrayHG38b.txt data/MotifPipeline/remap/GM12878.bed -names me chip" > GM12878_overlap_NOmotif.txt
# eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wao -a data/MotifPipeline/ENCODE/wgbsin/GM12878both.txt -b data/MotifPipeline/ENCODE/methyl_array/crossQC/GM12878_shufMeArrayHG38b.txt data/MotifPipeline/remap/GM12878.bed -names wg other" >> GM12878_overlap_NOmotif.txt
eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a data/MotifPipeline/ENCODE/wgbsin/GM12878both.txt -b data/MotifPipeline/ENCODE/methyl_array/crossQC/GM12878_shufMeArrayHG38b.txt"|cut -f1,2,3,4,5,9 > GM12878_overlap_NOmotif_a.txt
eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wao -a GM12878_overlap_NOmotif_a.txt -b data/MotifPipeline/remap/GM12878.bed "|cut -f1,2,3,4,5,6,10,11 > GM12878_overlap_NOmotif_b.txt

# eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a data/MotifPipeline/ENCODE/wgbsin/K562both.txt -b data/MotifPipeline/ENCODE/methyl_array/crossQC/K-562_shufMeArrayHG38b.txt data/MotifPipeline/remap/K-562.bed -names me chip" > K562_overlap_NOmotif.txt
# eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wao -a data/MotifPipeline/ENCODE/wgbsin/K562both.txt -b data/MotifPipeline/ENCODE/methyl_array/crossQC/K-562_shufMeArrayHG38b.txt data/MotifPipeline/remap/K-562.bed -names wg other" >> K562_overlap_NOmotif.txt
eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a data/MotifPipeline/ENCODE/wgbsin/K562both.txt -b data/MotifPipeline/ENCODE/methyl_array/crossQC/K-562_shufMeArrayHG38b.txt"|cut -f1,2,3,4,5,9 > K562_overlap_NOmotif_a.txt
eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wao -a K562_overlap_NOmotif_a.txt -b data/MotifPipeline/remap/K-562.bed "|cut -f1,2,3,4,5,6,10,11 > K562_overlap_NOmotif_b.txt

# eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a data/MotifPipeline/ENCODE/wgbsin/SKNSHboth.txt -b data/MotifPipeline/ENCODE/methyl_array/crossQC/SK-N-SH_shufMeArrayHG38b.txt data/MotifPipeline/remap/SK-N-SH.bed -names me chip" > SKNSH_overlap_NOmotif.txt
# eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wao -a data/MotifPipeline/ENCODE/wgbsin/SKNSHboth.txt -b data/MotifPipeline/ENCODE/methyl_array/crossQC/SK-N-SH_shufMeArrayHG38b.txt data/MotifPipeline/remap/SK-N-SH.bed -names wg other" >> SKNSH_overlap_NOmotif.txt
eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a data/MotifPipeline/ENCODE/wgbsin/SKNSHboth.txt -b data/MotifPipeline/ENCODE/methyl_array/crossQC/SK-N-SH_shufMeArrayHG38b.txt"|cut -f1,2,3,4,5,9 > SKNSH_overlap_NOmotif_a.txt
eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wao -a SKNSH_overlap_NOmotif_a.txt -b data/MotifPipeline/remap/SK-N-SH.bed "|cut -f1,2,3,4,5,6,10,11 > SKNSH_overlap_NOmotif_b.txt


tfdb='data/MotifPipeline/ENCODE/Homo_sapiens_motifinfo.txt'
motifdir='../rekrg/MotifScans/MotifScans/hg38_bed/'
motifs=$(ls $motifdir*)
outdir='data/MotifPipeline/sthlm_motif'

for TF in $motifs
do
tf=$(eval "echo "$TF" | cut -d / -f6| cut -d _ -f1")

gene=$(awk -v pat=$tf '$1 ~ pat' $tfdb | cut -f 2)

eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a $TF -b ~/A549_overlap_NOmotif_b.txt"|cut -f1,2,3,5,8,9,10,11,12,13,14 > $outdir/A549_$gene
eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a $TF -b ~/HeLa_overlap_NOmotif_b.txt"|cut -f1,2,3,5,8,9,10,11,12,13,14 > $outdir/HeLa_$gene
eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a $TF -b ~/HepG2_overlap_NOmotif_b.txt" |cut -f1,2,3,5,8,9,10,11,12,13,14> $outdir/HepG2_$gene
eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a $TF -b ~/GM12878_overlap_NOmotif_b.txt"|cut -f1,2,3,5,8,9,10,11,12,13,14 > $outdir/GM12878_$gene
eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a $TF -b ~/K562_overlap_NOmotif_b.txt" |cut -f1,2,3,5,8,9,10,11,12,13,14> $outdir/K562_$gene
# eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a $TF -b ~/SKNSH_overlap_NOmotif_b.txt" |cut -f1,2,3,5,8,9,10,11,12,13,14> $outdir/SKNSH_$gene
find . -size 0 -delete

done


# cut -f1,2,3,8 A549_overlap_NOmotif.txt > 1238.txt

# XX='me'
# awk ' ($5 == $XX ) ' A549_overlap.txt  | cut -f1,2,3,4,5 |uniq >> meA549_overlap.txt
# XX='chip'
# awk ' ($5 == $XX ) ' A549_overlap.txt |cut -f1,2,3,4,5 |uniq >> chipA549_overlap.txt
# XX='wg'
# awk ' ($5 == $XX ) ' A549_overlap.txt |cut -f1,2,3,4,5 |uniq >> wgA549_overlap.txt
