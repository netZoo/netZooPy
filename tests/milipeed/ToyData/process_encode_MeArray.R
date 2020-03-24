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
# write.table(annotation@listData$chr,'data/MotifPipeline/ENCODE/methyl_array/chr.txt',sep ='\t',quote = FALSE,dec = '.',row.names=FALSE,col.names = FALSE)


##bash to format as bed file
# paste chr.txt start.txt A-549.txt | awk '{print $1,$2,$2+1,$4}' OFS='\t'> A-549_MeArray.txt
# paste chr.txt start.txt GM12878.txt | awk '{print $1,$2,$2+1,$4}' OFS='\t'> GM12878_MeArray.txt
# paste chr.txt start.txt HeLa-S3.txt | awk '{print $1,$2,$2+1,$4}' OFS='\t'> HeLa-S3_MeArray.txt
# paste chr.txt start.txt Hep-G2.txt | awk '{print $1,$2,$2+1,$4}' OFS='\t'> Hep-G2_MeArray.txt
# paste chr.txt start.txt K-562.txt | awk '{print $1,$2,$2+1,$4}' OFS='\t'> K-562_MeArray.txt
# paste chr.txt start.txt SK-N-SH.txt | awk '{print $1,$2,$2+1,$4}' OFS='\t'> SK-N-SH_MeArray.txt