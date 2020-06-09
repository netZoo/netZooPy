# library(minfi,IlluminaHumanMethylationEPICmanifest)#,IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
### for H1 450k
# BiocManager::install("minfiData")

# BiocManager::install("IlluminaHumanMethylation450kmanifest")
# BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library(IlluminaHumanMethylation450kmanifest)#,IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

RGSet<-minfi::read.metharray.exp("data/MotifPipeline/ENCODE/methyl_array/idats/450k/")
# phenoData <- pData(RGSet)
manifest <- getManifest(RGSet)
head(getProbeInfo(manifest))
MSet <- preprocessRaw(RGSet) 
MSet
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
RSet
beta <- getBeta(RSet)
GRset <- mapToGenome(RSet)
GRset

##QC
RGSet <- preprocessQuantile(RGSet)#, fixOutliers = TRUE,removeBadSamples = TRUE, badSampleCutoff = 10.5,quantileNormalize = TRUE, stratified = TRUE, mergeManifest = FALSE, sex = NULL)
# MSet.swan <- preprocessQuantile(RGSet)
# RSet2 <- ratioConvert(MSet.swan)

snps <- getSnpInfo(RGSet)
head(snps,10)
RGSet <- addSnpInfo(GRset)
##filter out snps
RGSet <- dropLociWithSnps(RGSet, snps=c("SBE","CpG"), maf=0)

###filter out cross-reactive probes
devtools::install_github("markgene/maxprobes")
library(maxprobes)
RGSet<-maxprobes::dropXreactiveLoci(RGSet)
write.table(beta[,1],'data/MotifPipeline/ENCODE/methyl_array/H1a.txt',sep ='\t',quote = FALSE,dec = '.',col.names = FALSE)


##else comment above out
beta <- getBeta(RGSet)
M <- getM(GRset)
CN <- getCN(GRset)
sampleNames <- sampleNames(GRset)
probeNames <- featureNames(GRset)
pheno <- pData(GRset)
gr <- granges(GRset)
annotation <- getAnnotation(RGSet)
qc <- getQC(MSet)
head(qc)
plotQC(qc)
densityPlot(beta, sampGroups = GRset@colData@rownames,main = 'preQC_beta',legend = None)
densityBeanPlot(beta, sampGroups = GRset@colData@rownames,main = 'preQC_beta')
controlStripPlot(RGSet, controls="BISULFITE CONVERSION II")
qcReport(RGSet, pdf= "qcReport.pdf",sampNames = MSet@colData@rownames,sampGroups = MSet@colData@rownames)

# write.table(beta[,1],'data/MotifPipeline/ENCODE/methyl_array/A-549b.txt',sep ='\t',quote = FALSE,dec = '.',col.names = FALSE)
# write.table(beta[,2],'data/MotifPipeline/ENCODE/methyl_array/GM12878b.txt',sep ='\t',quote = FALSE,dec = '.',col.names = FALSE)
# write.table(beta[,3],'data/MotifPipeline/ENCODE/methyl_array/HeLa-S3b.txt',sep ='\t',quote = FALSE,dec = '.',col.names = FALSE)
# write.table(beta[,4],'data/MotifPipeline/ENCODE/methyl_array/Hep-G2b.txt',sep ='\t',quote = FALSE,dec = '.',col.names = FALSE)
# write.table(beta[,5],'data/MotifPipeline/ENCODE/methyl_array/K562b.txt',sep ='\t',quote = FALSE,dec = '.',col.names = FALSE)
# write.table(beta[,6],'data/MotifPipeline/ENCODE/methyl_array/SKNSHb.txt',sep ='\t',quote = FALSE,dec = '.',col.names = FALSE)
write.table(RGSet@rowRanges@ranges@start,'data/MotifPipeline/ENCODE/methyl_array/startc.txt',sep ='\t',quote = FALSE,dec = '.',row.names=FALSE,col.names = FALSE)
write.table(annotation@listData$chr,'data/MotifPipeline/ENCODE/methyl_array/chrc.txt',sep ='\t',quote = FALSE,dec = '.',row.names=FALSE,col.names = FALSE)
write.table(annotation$Relation_to_Island,'data/MotifPipeline/ENCODE/methyl_array/R2islandc.txt',sep ='\t',quote = FALSE,dec = '.',row.names=FALSE,col.names = FALSE)



# ##bash to format as bed file cd to data/MotifPipeline/ENCODE/methyl_array/
# paste chrb.txt startb.txt A-549b.txt | awk '{print $1,$2,$2+1,$4}' OFS='\t'> A-549_MeArrayb.txt
# paste chrb.txt startb.txt GM12878b.txt | awk '{print $1,$2,$2+1,$4}' OFS='\t'> GM12878_MeArrayb.txt
# paste chrb.txt startb.txt HeLa-S3b.txt | awk '{print $1,$2,$2+1,$4}' OFS='\t'> HeLa-S3_MeArrayb.txt
# paste chrb.txt startb.txt Hep-G2b.txt | awk '{print $1,$2,$2+1,$4}' OFS='\t'> Hep-G2_MeArrayb.txt
# paste chrb.txt startb.txt K562b.txt | awk '{print $1,$2,$2+1,$4}' OFS='\t'> K-562_MeArrayb.txt
# paste chrb.txt startb.txt SKNSHb.txt | awk '{print $1,$2,$2+1,$4}' OFS='\t'> SK-N-SH_MeArrayb.txt
paste chrc.txt startc.txt H1a.txt R2islandc.txt| awk '{print $1,$2,$2+1,$4,$5}' OFS='\t'> H1_MeArrayb.txt


# 
# ##map from 19 to 38 cd to data
# ./liftOver MotifPipeline/ENCODE/methyl_array/A-549_MeArrayb.txt hg19ToHg38.over.chain MotifPipeline/ENCODE/methyl_array/A-549_MeArrayHG38b.txt unmapped
# ./liftOver MotifPipeline/ENCODE/methyl_array/GM12878_MeArrayb.txt hg19ToHg38.over.chain MotifPipeline/ENCODE/methyl_array/GM12878_MeArrayHG38b.txt unmapped
# ./liftOver MotifPipeline/ENCODE/methyl_array/HeLa-S3_MeArrayb.txt hg19ToHg38.over.chain MotifPipeline/ENCODE/methyl_array/HeLa-S3_MeArrayHG38b.txt unmapped
# ./liftOver MotifPipeline/ENCODE/methyl_array/Hep-G2_MeArrayb.txt hg19ToHg38.over.chain MotifPipeline/ENCODE/methyl_array/Hep-G2_MeArrayHG38b.txt unmapped
# ./liftOver MotifPipeline/ENCODE/methyl_array/K-562_MeArrayb.txt hg19ToHg38.over.chain MotifPipeline/ENCODE/methyl_array/K-562_MeArrayHG38b.txt unmapped
# ./liftOver MotifPipeline/ENCODE/methyl_array/SK-N-SH_MeArrayb.txt hg19ToHg38.over.chain MotifPipeline/ENCODE/methyl_array/SK-N-SH_MeArrayHG38b.txt unmapped
# 

./data/liftOver MotifPipeline/ENCODE/methyl_array/H1_MeArrayb.txt hg19ToHg38.over.chain MotifPipeline/ENCODE/methyl_array/H1_MeArrayHG38b.txt unmapped

# ##build shuffled versions
# cut -f4 A-549_MeArrayHG38b.txt|cut -f4 |shuf> shuf_A549.txt
# paste A-549_MeArrayHG38b.txt shuf_A549.txt | awk '{print $1,$2,$3,$4,$5}' OFS='\t' > A-549_shufMeArrayHG38b.txt
# cut -f4 GM12878_MeArrayHG38b.txt|cut -f4 |shuf> shuf_GM12878.txt
# paste GM12878_MeArrayHG38b.txt shuf_GM12878.txt | awk '{print $1,$2,$3,$4,$5}' OFS='\t' > GM12878_shufMeArrayHG38b.txt
# cut -f4 HeLa-S3_MeArrayHG38b.txt|cut -f4 |shuf> shuf_HeLa-S3.txt
# paste HeLa-S3_MeArrayHG38b.txt shuf_HeLa-S3.txt | awk '{print $1,$2,$3,$4,$5}' OFS='\t' > HeLa-S3_shufMeArrayHG38b.txt
# cut -f4 Hep-G2_MeArrayHG38b.txt|cut -f4 |shuf> shuf_Hep-G2.txt
# paste Hep-G2_MeArrayHG38b.txt shuf_Hep-G2.txt | awk '{print $1,$2,$3,$4,$5}' OFS='\t' > Hep-G2_shufMeArrayHG38b.txt
# cut -f4 K-562_MeArrayHG38b.txt|cut -f4 |shuf> shuf_K-562.txt
# paste K-562_MeArrayHG38b.txt shuf_K-562.txt | awk '{print $1,$2,$3,$4,$5}' OFS='\t' > K-562_shufMeArrayHG38b.txt
# cut -f4 SK-N-SH_MeArrayHG38b.txt|cut -f4 |shuf> shuf_SK-N-SH.txt
# paste SK-N-SH_MeArrayHG38b.txt shuf_SK-N-SH.txt | awk '{print $1,$2,$3,$4,$5}' OFS='\t' > SK-N-SH_shufMeArrayHG38b.txt
# 
cut -f4 H1_MeArrayHG38b.txt|cut -f4 |shuf> shuf_H1.txt
paste H1_MeArrayHG38b.txt shuf_H1.txt | awk '{print $1,$2,$3,$4,$5}' OFS='\t' > H1_shufMeArrayHG38b.txt
# 
#
#
###VENN DIAGRAM
library(VennDiagram)
library(readr)
wgbs=read.csv('data/MotifPipeline/ENCODE/wgbsin/A549both.txt',sep='\t',header=FALSE)#, colClasses=c(NA, NA, NA,"NULL","NULL","NULL",)))#,usecols=[0,1,2],sep='\t',names=['chr','start','end'])
chip=read.csv('data/MotifPipeline/remap/A-549.bed',sep='\t',header=FALSE)#,usecols=[0,1,2],sep='\t',names=['chr','start','end'])
me=read.csv('data/MotifPipeline/ENCODE/methyl_array/crossQC/A-549_shufMeArrayHG38b.txt',sep='\t',header=FALSE)#,usecols=[0,1,2],sep='\t',names=['chr','start','end'])
me$V4<-NULL
me$V5<-NULL
chip$V4<-NULL

M1<-paste(me$V1,me$V2,me$V3,sep="_")
C1<-paste(chip$V1,chip$V2,chip$V3,sep="_")
W1<-paste(wgbs$V1,wgbs$V2,wgbs$V3,sep="_")
# sample four-set Venn Diagram
# Mo1<-read.csv('data/MotifPipeline/glassHG38.txt',sep='\t',header=FALSE)#,usecols=[0,1,2],sep='\t',names=['chr','start','end'])

venn.plot <- venn.diagram(
  x = list(
    '    A. WGBS' = W1,
    # motif = Mo1,
    'C. ChIP ' = C1,
    'B. meArray' = M1
  ),
  filename = "interX_Venn3.tiff",
  col = "transparent",
  fill = c("blue", "pink", "orange"),
  alpha = 0.50,
  # label.col = c("orange", "white", "darkorchid4", "white", 
  #               "white", "white", "white", "white", "darkblue", "white", 
  #               "white", "white", "white", "darkgreen", "white"),
  cex = 2,
  fontfamily = "serif",
  fontface = "bold",
  # cat.col = c("blue", "darkpink", "orange"),#, "darkorchid4"),
  cat.cex = 2,
  # cat.pos = 0,
  # cat.dist = 0.07,
  cat.fontfamily = "serif",
  # rotation.degree = 270,
  # margin = 0.2
);




#
#
#
# par(mfrow=c(6,2))
MSet.swan <- preprocessSWAN(RGSet)
RSet2 <- ratioConvert(MSet.swan)
beta_SWANQC <- getBeta(RSet2)
hist(beta_SWANQC)
A=densityPlot(MSet.swan, sampGroups = MSet.swan@colData@rownames,main = 'preQC_beta',legend = None)
B<-densityBeanPlot(MSet.swan, sampGroups = MSet.swan@colData@rownames,main = 'preQC_beta')
C=densityPlot(beta_SWANQC, sampGroups = MSet.swan@colData@rownames,main = 'SWANQC_beta',legend = None)
D=densityBeanPlot(beta_SWANQC, sampGroups = MSet.swan@colData@rownames,main = 'SWANQC_beta')
# 
# 
# 
# GRset.quantile <- preprocessQuantile(RGSet, fixOutliers = TRUE,removeBadSamples = TRUE, badSampleCutoff = 10.5,quantileNormalize = TRUE, stratified = TRUE, mergeManifest = FALSE, sex = NULL)
# # RSet2 <- ratioConvert(GRset.quantile)
# beta_quantileQC <- getBeta(GRset.quantile)
# hist(beta_quantileQC)
# 
# E=densityPlot(beta_quantileQC, sampGroups = MSet.swan@colData@rownames,main = 'quantilQC_beta',legend = None)
# F=densityBeanPlot(beta_quantileQC, sampGroups = MSet.swan@colData@rownames,main = 'quantilQC_beta')
# 
# densityPlot(beta_quantileQC, sampGroups = MSet@colData@rownames)
# densityBeanPlot(beta_quantileQC, sampGroups = MSet@colData@rownames)
# controlStripPlot(GRset.quantile, controls="BISULFITE CONVERSION II")
# qcReport(RGSet, pdf= "qcReport.pdf",sampNames = MSet@colData@rownames,sampGroups = MSet@colData@rownames)
# 
# ###swan == quantile ~~ funnorm
# ### illumina and noob are terrible -- off center mean
# 
# 
# GRset.funnorm <- preprocessFunnorm(RGSet)
# # RSet3 <- ratioConvert(MSet.noob)
# beta_funnormQC <- getBeta(GRset.funnorm)
# hist(beta_funnormQC)
# G=densityPlot(beta_funnormQC, sampGroups = MSet.swan@colData@rownames,main = 'beta_funnormQC',legend = None)
# H=densityBeanPlot(beta_funnormQC, sampGroups = MSet.swan@colData@rownames,main = 'beta_funnormQC')
# 
# 
# MSet.illumina <- preprocessIllumina(RGSet, bg.correct = TRUE,normalize = "controls")
# RSet <- ratioConvert(MSet.illumina)
# beta_illuminaQC <- getBeta(RSet)
# hist(beta_illuminaQC)
# # densityPlot(MSet.illumina, sampGroups = MSet.illumina@colData@rownames,main = 'pre-illuminaQC_beta')
# # densityBeanPlot(MSet.illumina, sampGroups = MSet.illumina@colData@rownames,main = 'pre-illuminaQC_beta')
# I<-densityPlot(beta_illuminaQC, sampGroups = MSet.swan@colData@rownames,main = 'beta_illuminaQC',legend = None)
# J<-densityBeanPlot(beta_illuminaQC, sampGroups = MSet.swan@colData@rownames,main = 'beta_illuminaQC')
# 
# 
# MSet.noob <- preprocessNoob(RGSet)
# RSet3 <- ratioConvert(MSet.noob)
# beta_noobQC <- getBeta(RSet3)
# hist(beta_noobQC)
# # densityPlot(MSet.noob, sampGroups = MSet.noob@colData@rownames,main = 'pre-noobQC_beta')
# # densityBeanPlot(MSet.noob, sampGroups = MSet.noob@colData@rownames,main = 'pre-noobQC_beta')
# K=densityPlot(beta_noobQC, sampGroups = MSet.swan@colData@rownames,main = 'beta_noobQC',legend = None)
# L=densityBeanPlot(beta_noobQC, sampGroups = MSet.swan@colData@rownames,main = 'beta_noobQC')
# 
# # ggarrange(A,C,E,G,I,K,#B,C,D,E,F,G,H,I,J,K,L, # labels = c("A", "B", "C"),ncol = 2, nrow = 6)
# 
# # write.table(beta_SWANQC[,1],'data/MotifPipeline/ENCODE/methyl_array/A-549_QC.txt',sep ='\t',quote = FALSE,dec = '.',col.names = FALSE)
# # write.table(beta_SWANQC[,2],'data/MotifPipeline/ENCODE/methyl_array/GM12878_QC.txt',sep ='\t',quote = FALSE,dec = '.',col.names = FALSE)
# # write.table(beta_SWANQC[,3],'data/MotifPipeline/ENCODE/methyl_array/HeLa-S3_QC.txt',sep ='\t',quote = FALSE,dec = '.',col.names = FALSE)
# # write.table(beta_SWANQC[,4],'data/MotifPipeline/ENCODE/methyl_array/Hep-G2_QC.txt',sep ='\t',quote = FALSE,dec = '.',col.names = FALSE)
# # write.table(beta_SWANQC[,5],'data/MotifPipeline/ENCODE/methyl_array/K562_QC.txt',sep ='\t',quote = FALSE,dec = '.',col.names = FALSE)
# # write.table(beta_SWANQC[,6],'data/MotifPipeline/ENCODE/methyl_array/SKNSH_QC.txt',sep ='\t',quote = FALSE,dec = '.',col.names = FALSE)
# # 
# # 
# # # getManifest(rgSet)
# # # MSet.noobs <- preprocessNoob(rgSet)
# # # # MSet <- preprocessRaw(rgSet) 
# # # # qc <- getQC(MSet.illumina)
# # # # plotQC(qc)
# # # RSet <- ratioConvert(MSet.noobs)
# # # # RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
# # # beta <- getBeta(RSet)
# # # GRset <- mapToGenome(RSet)
# # # annotation <- getAnnotation(GRset)
# # # names(annotation)
# # # gr <- granges(GRset)
# # # write.table(GRset@assays@data$Beta,'data/MotifPipeline/ENCODE/methyl_array/SKNSH.txt',sep ='\t',quote = FALSE,dec = '.',col.names = FALSE)
# # # # write.table(GRset@rowRanges@ranges@start,'data/MotifPipeline/ENCODE/methyl_array/start.txt',sep ='\t',quote = FALSE,dec = '.',row.names=FALSE,col.names = FALSE)
# # # # write.table(annotation@listData$AddressA,'data/MotifPipeline/ENCODE/methyl_array/startA.txt',sep ='\t',quote = FALSE,dec = '.',row.names=FALSE,col.names = FALSE)
# # # 
# # # # write.table(annotation@listData$chr,'data/MotifPipeline/ENCODE/methyl_array/chr.txt',sep ='\t',quote = FALSE,dec = '.',row.names=FALSE,col.names = FALSE)
# # 
# # 
# # ##bash to format as bed file cd to data/MotifPipeline/ENCODE/methyl_array/
# # paste chr.txt start.txt A-549_QC.txt | awk '{print $1,$2,$2+1,$4}' OFS='\t'> A-549_MeArrayQC.txt
# # paste chr.txt start.txt GM12878_QC.txt | awk '{print $1,$2,$2+1,$4}' OFS='\t'> GM12878_MeArrayQC.txt
# # paste chr.txt start.txt HeLa-S3_QC.txt | awk '{print $1,$2,$2+1,$4}' OFS='\t'> HeLa-S3_MeArrayQC.txt
# # paste chr.txt start.txt Hep-G2_QC.txt | awk '{print $1,$2,$2+1,$4}' OFS='\t'> Hep-G2_MeArrayQC.txt
# # paste chr.txt start.txt K562_QC.txt | awk '{print $1,$2,$2+1,$4}' OFS='\t'> K-562_MeArrayQC.txt
# # paste chr.txt start.txt SKNSH_QC.txt | awk '{print $1,$2,$2+1,$4}' OFS='\t'> SK-N-SH_MeArrayQC.txt
# # 
# # 
# # ##map from 19 to 38 cd to data
# # ./liftOver MotifPipeline/ENCODE/methyl_array/A-549_MeArrayQC.txt hg19ToHg38.over.chain MotifPipeline/ENCODE/methyl_array/A-549_MeArrayHG38QC.txt unmapped
# # ./liftOver MotifPipeline/ENCODE/methyl_array/GM12878_MeArrayQC.txt hg19ToHg38.over.chain MotifPipeline/ENCODE/methyl_array/GM12878_MeArrayHG38QC.txt unmapped
# # ./liftOver MotifPipeline/ENCODE/methyl_array/HeLa-S3_MeArrayQC.txt hg19ToHg38.over.chain MotifPipeline/ENCODE/methyl_array/HeLa-S3_MeArrayHG38QC.txt unmapped
# # ./liftOver MotifPipeline/ENCODE/methyl_array/Hep-G2_MeArrayQC.txt hg19ToHg38.over.chain MotifPipeline/ENCODE/methyl_array/Hep-G2_MeArrayHG38QC.txt unmapped
# # ./liftOver MotifPipeline/ENCODE/methyl_array/K-562_MeArrayQC.txt hg19ToHg38.over.chain MotifPipeline/ENCODE/methyl_array/K-562_MeArrayHG38QC.txt unmapped
# # ./liftOver MotifPipeline/ENCODE/methyl_array/SK-N-SH_MeArrayQC.txt hg19ToHg38.over.chain MotifPipeline/ENCODE/methyl_array/SK-N-SH_MeArrayHG38QC.txt unmapped
# # 
# # ##build shuffled versions
# # cut -f4 A-549_MeArrayHG38QC.txt|cut -f4 |shuf> shuf_A549.txt
# # paste A-549_MeArrayHG38QC.txt shuf_A549.txt | awk '{print $1,$2,$3,$4,$5}' OFS='\t' > A-549_shufMeArrayHG38QC.txt
# # cut -f4 GM12878_MeArrayHG38QC.txt|cut -f4 |shuf> shuf_GM12878.txt
# # paste GM12878_MeArrayHG38QC.txt shuf_GM12878.txt | awk '{print $1,$2,$3,$4,$5}' OFS='\t' > GM12878_shufMeArrayHG38QC.txt
# # cut -f4 HeLa-S3_MeArrayHG38QC.txt|cut -f4 |shuf> shuf_HeLa-S3.txt
# # paste HeLa-S3_MeArrayHG38QC.txt shuf_HeLa-S3.txt | awk '{print $1,$2,$3,$4,$5}' OFS='\t' > HeLa-S3_shufMeArrayHG38QC.txt
# # cut -f4 Hep-G2_MeArrayHG38QC.txt|cut -f4 |shuf> shuf_Hep-G2.txt
# # paste Hep-G2_MeArrayHG38QC.txt shuf_Hep-G2.txt | awk '{print $1,$2,$3,$4,$5}' OFS='\t' > Hep-G2_shufMeArrayHG38QC.txt
# # cut -f4 K-562_MeArrayHG38QC.txt|cut -f4 |shuf> shuf_K-562.txt
# # paste K-562_MeArrayHG38QC.txt shuf_K-562.txt | awk '{print $1,$2,$3,$4,$5}' OFS='\t' > K-562_shufMeArrayHG38QC.txt
# # cut -f4 SK-N-SH_MeArrayHG38QC.txt|cut -f4 |shuf> shuf_SK-N-SH.txt
# # paste SK-N-SH_MeArrayHG38QC.txt shuf_SK-N-SH.txt | awk '{print $1,$2,$3,$4,$5}' OFS='\t' > SK-N-SH_shufMeArrayHG38QC.txt
# # 
# # 
