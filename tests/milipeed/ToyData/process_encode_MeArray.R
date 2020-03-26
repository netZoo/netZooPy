library(minfi,IlluminaHumanMethylationEPICmanifest,IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
rgSet<-minfi::read.metharray("data/MotifPipeline/ENCODE/methyl_array/SKNSH/idat")
getManifest(rgSet)
MSet.illumina <- preprocessNoob(rgSet)
# MSet <- preprocessRaw(rgSet) 
# qc <- getQC(MSet.illumina)
# plotQC(qc)
RSet <- ratioConvert(MSet.illumina, what = "both", keepCN = TRUE)
# RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
beta <- getBeta(RSet)
GRset <- mapToGenome(RSet)
annotation <- getAnnotation(GRset)
names(annotation)
gr <- granges(GRset)
write.table(GRset@assays@data$Beta,'data/MotifPipeline/ENCODE/methyl_array/SKNSH.txt',sep ='\t',quote = FALSE,dec = '.',col.names = FALSE)
# write.table(GRset@rowRanges@ranges@start,'data/MotifPipeline/ENCODE/methyl_array/start.txt',sep ='\t',quote = FALSE,dec = '.',row.names=FALSE,col.names = FALSE)
write.table(annotation@listData$AddressA,'data/MotifPipeline/ENCODE/methyl_array/startA.txt',sep ='\t',quote = FALSE,dec = '.',row.names=FALSE,col.names = FALSE)

# write.table(annotation@listData$chr,'data/MotifPipeline/ENCODE/methyl_array/chr.txt',sep ='\t',quote = FALSE,dec = '.',row.names=FALSE,col.names = FALSE)


##bash to format as bed file cd to data/MotifPipeline/ENCODE/methyl_array/
paste chr.txt start.txt A-549.txt | awk '{print $1,$2,$2+1,$4}' OFS='\t'> A-549_MeArray.txt
paste chr.txt start.txt GM12878.txt | awk '{print $1,$2,$2+1,$4}' OFS='\t'> GM12878_MeArray.txt
paste chr.txt start.txt HeLa-S3.txt | awk '{print $1,$2,$2+1,$4}' OFS='\t'> HeLa-S3_MeArray.txt
paste chr.txt start.txt Hep-G2.txt | awk '{print $1,$2,$2+1,$4}' OFS='\t'> Hep-G2_MeArray.txt
paste chr.txt start.txt K-562.txt | awk '{print $1,$2,$2+1,$4}' OFS='\t'> K-562_MeArray.txt
paste chr.txt start.txt SK-N-SH.txt | awk '{print $1,$2,$2+1,$4}' OFS='\t'> SK-N-SH_MeArray.txt


##map from 19 to 38 cd to data
# ./liftOver MotifPipeline/ENCODE/methyl_array/A-549_MeArray.txt hg19ToHg38.over.chain MotifPipeline/ENCODE/methyl_array/A-549_MeArrayHG38.txt unmapped
# ./liftOver MotifPipeline/ENCODE/methyl_array/GM12878_MeArray.txt hg19ToHg38.over.chain MotifPipeline/ENCODE/methyl_array/GM12878_MeArrayHG38.txt unmapped
# ./liftOver MotifPipeline/ENCODE/methyl_array/HeLa-S3_MeArray.txt hg19ToHg38.over.chain MotifPipeline/ENCODE/methyl_array/HeLa-S3_MeArrayHG38.txt unmapped
# ./liftOver MotifPipeline/ENCODE/methyl_array/Hep-G2_MeArray.txt hg19ToHg38.over.chain MotifPipeline/ENCODE/methyl_array/Hep-G2_MeArrayHG38.txt unmapped
# ./liftOver MotifPipeline/ENCODE/methyl_array/K-562_MeArray.txt hg19ToHg38.over.chain MotifPipeline/ENCODE/methyl_array/K-562_MeArrayHG38.txt unmapped
# ./liftOver MotifPipeline/ENCODE/methyl_array/SK-N-SH_MeArray.txt hg19ToHg38.over.chain MotifPipeline/ENCODE/methyl_array/SK-N-SH_MeArrayHG38.txt unmapped

##build shuffled versions
cut -f4 A-549_MeArrayHG38.txt|cut -f4 |shuf> shuf_A549.txt
paste A-549_MeArrayHG38.txt shuf_A549.txt | awk '{print $1,$2,$3,$4,$5}' OFS='\t' > A-549_shufMeArrayHG38.txt
cut -f4 GM12878_MeArrayHG38.txt|cut -f4 |shuf> shuf_GM12878.txt
paste GM12878_MeArrayHG38.txt shuf_GM12878.txt | awk '{print $1,$2,$3,$4,$5}' OFS='\t' > GM12878_shufMeArrayHG38.txt
cut -f4 HeLa-S3_MeArrayHG38.txt|cut -f4 |shuf> shuf_HeLa-S3.txt
paste HeLa-S3_MeArrayHG38.txt shuf_HeLa-S3.txt | awk '{print $1,$2,$3,$4,$5}' OFS='\t' > HeLa-S3_shufMeArrayHG38.txt
cut -f4 Hep-G2_MeArrayHG38.txt|cut -f4 |shuf> shuf_Hep-G2.txt
paste Hep-G2_MeArrayHG38.txt shuf_Hep-G2.txt | awk '{print $1,$2,$3,$4,$5}' OFS='\t' > Hep-G2_shufMeArrayHG38.txt
cut -f4 K-562_MeArrayHG38.txt|cut -f4 |shuf> shuf_K-562.txt
paste K-562_MeArrayHG38.txt shuf_K-562.txt | awk '{print $1,$2,$3,$4,$5}' OFS='\t' > K-562_shufMeArrayHG38.txt
cut -f4 SK-N-SH_MeArrayHG38.txt|cut -f4 |shuf> shuf_SK-N-SH.txt
paste SK-N-SH_MeArrayHG38.txt shuf_SK-N-SH.txt | awk '{print $1,$2,$3,$4,$5}' OFS='\t' > SK-N-SH_shufMeArrayHG38.txt


