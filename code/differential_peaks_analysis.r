#######################################################################################
##  Script to predict the differential peaks. 										  #
##  Here is the script to analyze h4k16ac data                                        #
#######################################################################################

# load library
library(GenomicRanges)
library(rtracklayer)
library(IRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Gviz)
library(biomaRt)
library(Rsubread)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


###############################################################################################
# readin the peaks of macs2 analysis
inp=import.bed(file.path("output_h3k27ac", "all.bed") ) 

# read in the roadmap peaks
inp1=import.bed(file.path("roadmap", "BI.Brain_Angular_Gyrus.H3K27ac.149.bed") )  
inp2=import.bed(file.path("roadmap", "BI.Brain_Anterior_Caudate.H3K27ac.149.bed") ) 
inp3=import.bed(file.path("roadmap", "BI.Brain_Cingulate_Gyrus.H3K27ac.149.bed") ) 
inp4=import.bed(file.path("roadmap", "BI.Brain_Hippocampus_Middle.H3K27ac.150.bed") ) 
inp5=import.bed(file.path("roadmap", "BI.Brain_Inferior_Temporal_Lobe.H3K27ac.149.bed") ) 
inp6=import.bed(file.path("roadmap", "BI.Brain_Mid_Frontal_Lobe.H3K27ac.149.bed") ) 
inp7=import.bed(file.path("roadmap", "BI.Brain_Substantia_Nigra.H3K27ac.149.bed") )

# validate the peaks identified by macs2
cutoff=0.6
lec=100
o1=findOverlaps(inp,inp1,minoverlap=lec)
overlaps <- pintersect(inp[queryHits(o1)], inp1[subjectHits(o1)])
percentOverlap <- width(overlaps) / width(inp1[subjectHits(o1)])
hits1 <- o1[percentOverlap > cutoff]
c1=hits1@from
o2=findOverlaps(inp,inp2,minoverlap=lec)
overlaps <- pintersect(inp[queryHits(o2)], inp2[subjectHits(o2)])
percentOverlap <- width(overlaps) / width(inp2[subjectHits(o2)])
hits2 <- o2[percentOverlap > cutoff]
c2=hits2@from
o3=findOverlaps(inp,inp3,minoverlap=lec)
overlaps <- pintersect(inp[queryHits(o3)], inp3[subjectHits(o3)])
percentOverlap <- width(overlaps) / width(inp3[subjectHits(o3)])
hits3 <- o3[percentOverlap > cutoff]
c3=hits3@from
o4=findOverlaps(inp,inp4,minoverlap=lec)
overlaps <- pintersect(inp[queryHits(o4)], inp4[subjectHits(o4)])
percentOverlap <- width(overlaps) / width(inp4[subjectHits(o4)])
hits4 <- o4[percentOverlap > cutoff]
c4=hits4@from
o5=findOverlaps(inp,inp5,minoverlap=lec)
overlaps <- pintersect(inp[queryHits(o5)], inp5[subjectHits(o5)])
percentOverlap <- width(overlaps) / width(inp5[subjectHits(o5)])
hits5 <- o5[percentOverlap > cutoff]
c5=hits5@from
o6=findOverlaps(inp,inp6,minoverlap=lec)
overlaps <- pintersect(inp[queryHits(o6)], inp6[subjectHits(o6)])
percentOverlap <- width(overlaps) / width(inp6[subjectHits(o6)])
hits6 <- o6[percentOverlap > cutoff]
c6=hits6@from
o7=findOverlaps(inp,inp7,minoverlap=lec)
overlaps <- pintersect(inp[queryHits(o7)], inp7[subjectHits(o7)])
percentOverlap <- width(overlaps) / width(inp7[subjectHits(o7)])
hits7 <- o7[percentOverlap > cutoff]
c7=hits7@from
olp=unique(c(c1,c2,c3,c4,c5,c6,c7))
new.inp=inp[olp,]

# transform peaks into data.frame format
ann.ext=as.data.frame(new.inp) 
ann.ext=cbind(GeneID=paste("peak",1:nrow(ann.ext),sep=""),ann.ext)
names(ann.ext) <-c("GeneID","Chr","Start","End","Width","Strand")

##################################################################################################

#read in sample annotation
an=read.table("sample.csv",sep=",",header=T)
input.files=sapply(as.vector(an$Run), function(x) file.path("h3k27ac", paste(x,".filtered.sorted.nodup.bam",sep="")))
cc <- featureCounts(input.files,annot.ext=ann.ext,isPairedEnd=F,allowMultiOverlap=T,primaryOnly=T,nthreads=6)
mx<-cc$counts
temp=colnames(mx)
temp=sapply(temp, function(x) strsplit(x,"\\.")[[1]][7])
colnames(mx) <-temp
tag=apply(mx, 1, function(x) length(x[x==0]))/ncol(mx) < 0.1
mx=mx[tag,]
age2<-cut(as.vector(an$age), breaks=5, ordered_result = T,labels=c("a","b","c","d","e"))
np2<-cut(as.vector(an$np), breaks=5, ordered_result = T,labels=c("a","b","c","d","e"))
np2[is.na(np2)]<-"c"
sex=an$sex
row.names(an)<-as.vector(an$Run)

#############################################################

# 1st round differential peak analysis
library(edgeR)
group=factor(as.vector(an$type),levels=c("Control","AD"))
ADcountList<-DGEList(counts=mx[,as.vector(an$Run)], group=group)
design<-model.matrix(~ age2 + np2 + sex +group) 
ADcountList_n <- calcNormFactors(ADcountList)
ADcountList_nd<-estimateDisp(ADcountList_n, design)
ee=cpm(ADcountList_nd)
fit <- glmQLFit(ADcountList_nd, design)
qlf <- glmQLFTest(fit,"groupAD")
res=subset(as.data.frame(topTags(qlf, 100000)),PValue < 0.01  ) # result of the first-round differential peaks analysis 
dif=row.names(res) # differential peaks
dr0=as.vector(res$logFC) # fold changes


# next round of differential peak analysis
# in this process, the differential peaks were removed for a new round analysis
new.mx=mx
pp=0;
while(length(dif)!=0){ # to remove the effects of regulation loss by multiple round of DE analysis
  pp=pp+1
  print(paste("Run round ",pp," ...",sep=""))
  print(length(dif))
  new.mx=new.mx[! row.names(new.mx) %in% dif,as.vector(sub.an$Run)] 
  ADcountList2<-DGEList(counts=new.mx, group=group)
  ADcountList_n2 <- calcNormFactors(ADcountList2)
  ADcountList_nd2 <-estimateDisp(ADcountList_n2, design)
  fit2 <- glmQLFit(ADcountList_nd2, design)
  qlf2 <- glmQLFTest(fit2,"groupAged, diseased")
  #qlf2$table$QValue <- p.adjust(qlf2$table$PValue, method="BH")
  res2=subset(as.data.frame(topTags(qlf2, 100000)),PValue <0.05) # to minimize the effects of DPs
  dif=row.names(res2)
  drr=as.vector(res2$logFC)
  #print(quantile(drr))
}

# update the parameters 
nf=ADcountList_n2$samples$norm.factors 
ADcountList_n$samples$norm.factors=nf
ADcountList_n$samples$lib.size=ADcountList_n2$samples$lib.size
ADcountList_nd<-estimateDisp(ADcountList_n, design)
fit <- glmQLFit(ADcountList_nd, design)
qlf <- glmQLFTest(fit,"groupAged, diseased")
res=subset(as.data.frame(topTags(qlf, 100000)),PValue < 0.05 & FDR < 0.05 )
dif=row.names(res)
dr=as.vector(res$logFC)
ee=cpm(ADcountList_nd, log=T) # normalized expression matrix 




# the ratio of regulation loss against regulation gain
library("RColorBrewer")
tb=table(dr < 0)
names(tb) <-c("AD","Normal")
pdf("bar_peak_h4k16ac.pdf",height=4,width=4)
barplot(tb,col=brewer.pal(n = 3, name = "RdBu")[c(1,3)])
legend("topleft",paste("Ratio=",round(tb[1]/tb[2],digits=2),": 1"))
dev.off()

# annoate the differential peaks
used=row.names(subset(subset(as.data.frame(topTags(qlf, 100000)),PValue <0.01 ), logFC < 0))
used=as.numeric(sapply(used, function(x) sub("peak","",x)))
peak=inp[used,]
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(clusterProfiler)
library(org.Hs.eg.db)

peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
pdf("peak_ann_h4k16ac.pdf", height=5,width=6)
plotAnnoPie(peakAnno)
dev.off()


