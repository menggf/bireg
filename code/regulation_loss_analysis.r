## load the R packages
library(parallel)
library(clusterProfiler)
library(devtools)
library(edgeR)
library(DESeq2)
library(MASS)

## read the AD data
da<-read.table("AMP-AD_MSBB_MSSM_IlluminaHiSeq2500_raw_counts_September_2016.txt")
clinical=read.table("AMP-AD_MSBB_MSSM_clinicalInfo_RNAseq_Dec2015_release.csv",header=T,sep="\t");
covariate=read.table("MSBB_RNAseq_covariates.csv",sep=",",header=T)

## data processing for normalized expression data
tag=apply(da,1, function(x) if(length(x[x<=10]) > 0.1* length(x)) FALSE else TRUE)
cc=as.matrix(da[tag,])
y1=DGEList(counts=cc)
y1 <- calcNormFactors(y1)
cpm(y1, log=T)->preexp
pas=colnames(preexp)

clinical=clinical[!is.na(clinical$final),]  #find the subjects with both RNA-seq and clinical 
row.names(clinical)<-paste("X",as.vector(clinical$final),sep="")
covariate=unique(subset(covariate, fileType=="bam" & batch !=""))
row.names(covariate)<-paste("X", as.vector(covariate$barcode),sep="")
used.pas=pas[(pas %in% row.names(clinical)) & (pas %in% row.names(covariate))]

preexp=preexp[,used.pas];
clinical=cbind(clinical[used.pas,], covariate[used.pas, ])
age=as.vector(clinical$AOD);
age[age=="90+"]=95;
age=as.numeric(age);
clinical$age=age
clinical$NP.1=as.factor(clinical$NP.1)

## differential expression expression
#diff=lapply(c("BM10","BM22","BM36","BM44"), function(tis){
#  testids=row.names(sub.clin)
#  dds <- DESeqDataSetFromMatrix(countData = cc[,testids],
#                              colData = sub.clin,
#                              design= ~ NP.1 + RACE + sex + pmi + batch )
#   dds <- DESeq(dds)
#   res <- results(dds, contrast=c("NP.1",2,1))  
#   res
#})

## remove the covariates of the expression data
exp=t(apply(preexp, 1, function(x){
	ip=clinical
	ip$ee=x
	fit = lm(ee ~ tissue + PMI + RACE + age + batch + RIN, data=ip)
	return(resid(fit))
}))


# readin TF genes
tfs=as.vector(as.matrix(read.table("tf_ensembl.txt")))
tfloc=which(row.names(exp) %in% tfs)

# bi-clustering analysis
cutoff=0.8
min.pas=50
install_github("menggf/bireg")
bi.results=bca(exp, tfloc, min.pas=min.pas, min.ges=50, max.ges=1000, cutoff=cutoff, cores=50)
res=transform(bi.results)

# find the missed or weakened regulators 
patients=colnames(exp)
genes=row.names(exp)
output=mclapply(tfloc, function(tf){ # check the TFs one by one
	pas<-res[[as.character(tf)]][["pas"]] 
	nges<-res[[as.character(tf)]][["nges"]]
	brk=which.max(nges) # the point reaching to the maximum number of regulated genes
	select.pas=patients[-1*pas[1:brk]] # the patients with regulation
	bk=length(patients)-which(nges > 20)[1] # check if enough genes were regulated
	rr1=apply(exp[, select.pas], 1, function(x) cor(x, exp[tf, select.pas], method="spearman")) 
	select.ges=genes[rr > cutoff]
  	select.ges=select.ges[select.ges != genes[tf] ]	    
	rr2=abs(sapply(select.ges, function(x) cor(exp[tf, !patients %in%select.pas], exp[x, !patients %in%select.pas], 
			method="spearman")))
	type="NR"
	if(quantile(rr2)[4] > 0.6 | quantile(rr2)[2] > 0.5) 
		type="DR" 
	else if(quantile(r2)[4] < 0.3) 
		type="MR" 
	else if(quantile(r2)[4] >=0.3 & quantile(r2)[4] <= 0.6) 
		type="WR";
	select.ges=genes[abs(rr) > 0.85 ]
	select.ges=select.ges[select.ges != genes[tf] ] #find the regulated genes
	cl1=clin[select.pas,] # the clinical for subjects with TF rgulation
	cl2=clin[patients[!patients %in% select.pas], ] # the clinical for subjects without TF rgulation
	cl=rbind(cbind(cl1,label=1),cbind(cl2, label = 2 )) # label them different group
	tb.tissue=table(cl[,c(5, 15)])
	p.tissue=chisq.test(tb.tissue)$p.value
	p.cdr1=vector(); #p-value
	p.plaq1=vector(); 
	p.bbscore1=vector(); 
	p.age1=vector(); 
	p.np11=vector();
	p.sex1=vector(); 
	f.cdr1=vector(); # fold changes
	f.plaq1=vector(); 
	f.bbscore1=vector(); 
	f.age1=vector(); 
	f.np11=vector();
	f.sex1=vector(); 
	for(tiss in unique(as.vector(clin$tissue))){ # evaluation in each brain regions
		sub.cl1=subset(cl1, tissue==tiss)
		sub.cl2=subset(cl2, tissue==tiss)
		sub.cl=subset(cl, tissue==tiss)
		if(dim(sub.cl1)[1] <= 30 | dim(sub.cl2)[1] <=30){ #set the p-value to 1 if there is not missed or weakened regulation
			p.cdr1=append(p.cdr1, 1)
			p.plaq1=append(p.plaq1, 1)
			p.bbscore1=append(p.bbscore1, 1)
			p.age1=append(p.age1, 1 )
			p.np11=append(p.np11, 1)
			p.sex=append(p.age, 1)
			
			f.cdr1=append(f.cdr1, 0)
			f.plaq1=append(f.plaq1, 0)
			f.bbscore1=append(f.bbscore1,0)
			f.age1=append(f.age1, 0 )
			f.np11=append(f.np11, 0)
		} else {
			p.cdr1=append(p.cdr1, ks.test(sub.cl1$CDR,sub.cl2$CDR)$p.value)
			p.plaq1=append(p.plaq1, ks.test(sub.cl1$PlaqueMean, sub.cl2$PlaqueMean)$p.value)
			p.bbscore1=append(p.bbscore1, ks.test(sub.cl1$bbscore, sub.cl2$bbscore)$p.value)
			p.age1=append(p.age1, ks.test(sub.cl1$age, sub.cl2$age)$p.value )
			p.np11=append(p.np11, ks.test(sub.cl1$NP.1, sub.cl2$NP.1)$p.value)
			tb.sex=table(sub.cl[,c(6, 11)])
			p.sex=append(p.np11, chisq.test(tb.sex)$p.value)
			
			f.cdr1=append(f.cdr1, mean(sub.cl1$CDR,na.rm = T)/mean(sub.cl2$CDR,,na.rm = T))
			f.plaq1=append(f.plaq1, mean(sub.cl1$PlaqueMean,,na.rm = T)/mean(sub.cl2$PlaqueMean,,na.rm = T))
			f.bbscore1=append(f.bbscore1,mean(sub.cl1$bbscore,,na.rm = T)/mean(sub.cl2$bbscore,,na.rm = T))
			f.age1=append(f.age1, mean(sub.cl1$age, ,na.rm = T)/mean(sub.cl2$age,,na.rm = T) )
			f.np11=append(f.np11, mean(sub.cl1$NP.1, ,na.rm = T)/mean(sub.cl2$NP.1, ,na.rm = T))
		}
	}
	# evaluation using all brain regions
	cl1=unique(cl1[,c(-2,-3,-4, -5)])
	cl2=unique(cl2[,c(-2,-3,-4, -5)])
	cl=rbind(cbind(cl1,label=1),cbind(cl2, label=2))
	p.cdr=ks.test(cl1$CDR, cl2$CDR)$p.value
	f.cdr=mean(cl1$CDR)/mean(cl2$CDR)
	p.plaq=ks.test(cl1$PlaqueMean, cl2$PlaqueMean)$p.value
	f.plaq=mean(cl1$PlaqueMean)/mean(cl2$PlaqueMean)
	p.bbscore=ks.test(cl1$bbscore, cl2$bbscore)$p.value
	f.bbscore=mean(cl1$bbscore,na.rm = T)/mean(cl2$bbscore,na.rm = T)
	p.age=ks.test(cl1$age, cl2$age)$p.value 
	p.np1=ks.test(cl1$NP.1, cl2$NP.1)$p.value
	f.np1=mean(cl1$NP.1)/mean(cl2$NP.1)
	tb.sex=table(cl[,c(6, 11)])
	p.sex=chisq.test(tb.sex)$p.value
	tb.race=table(cl[,c(3, 11)])
	p.race=chisq.test(tb.race)$p.value
	return(c(tf, type=type,length(select.ges),length(select.pas),bk, p.cdr,f.cdr, p.plaq,f.plaq, p.bbscore, 
		f.bbscore, p.age, p.np1, f.np1, p.tissue, p.sex, p.race, p.cdr1[1],f.cdr1[1], p.plaq1[1],f.plaq1[1], 
		p.bbscore1[1] ,f.bbscore1[1] ,p.age1[1],f.age1[1], p.np11[1],f.np11[1], p.cdr1[2],f.cdr1[2], 
		p.plaq1[2],f.plaq1[2], p.bbscore1[2],f.bbscore1[2] ,p.age1[2],f.age1[2], p.np11[2],f.np11[2], p.cdr1[3],
		f.cdr1[3], p.plaq1[3],f.plaq1[3],p.bbscore1[3],f.bbscore1[3] ,p.age1[3],f.age1[3], p.np11[3],f.np11[3], 
		p.cdr1[4],f.cdr1[4], p.plaq1[4],f.plaq1[4], p.bbscore1[4],f.bbscore1[4] ,p.age1[4],f.age1[4], p.np11[4], 
		f.np11[4]))
}, mc.cores=40)
names(output)<-as.character(tfloc)

report=as.data.frame(t(sapply(1:length(tfloc), function(i) temp[[i]])))
names(report)<-c("TF","type","no.ges", "no.pas", "biggest","cdr","f.cdr", "plaq", "f.plaq", "bbscore","f.bbscore",
	"age", "NP1","f.np1","tissue","sex", "race", "cdr.bm10", "f.cdr.bm10","plaq.bm10","f.plaq.bm10", "bbscore.bm10",
	"f.bbscore.bm10", "age.bm10","f.age.bm10", "NP1.bm10","f.NP1.bm10", "cdr.bm22","f.cdr.bm22", "plaq.bm22", 
	"f.plaq.bm22","bbscore.bm22","f.bbscore.bm22", "age.bm22","f.age.bm22", "NP1.bm22","f.NP1.bm22",  "cdr.bm36",
	"f.cdr.bm36", "plaq.bm36","f.plaq.bm36", "bbscore.bm36","f.bbscore.bm36", "age.bm36","f.age.bm36","NP1.bm36",
	"f.NP1.bm36", "cdr.bm44", "f.cdr.bm44", "plaq.bm44","f.plaq.bm44", "bbscore.bm44","f.bbscore.bm44", "age.bm44",
	"f.age.bm44", "NP1.bm44","f.NP1.bm44")

ids=as.character(report$TF)
mytype=tf.type[ids]
pre.tf=genes[as.vector(report$TF)]
tf.symbol=as.vector(as.matrix(ann[pre.tf, 2])) # tranfrom the tf location to symbols
report=cbind(symbol=tf.symbol, diff=diffs[ids], report)
cf=0.001 # cutoff for AD-related regulation loss
sub.report=subset(report,  cdr.bm10 <cf | plaq.bm10 <cf |bbscore.bm10 <cf |cdr.bm22 <cf | plaq.bm22 <cf | bbscore.bm22 < cf 
	|cdr.bm36 <cf | plaq.bm36 <cf |bbscore.bm36 <cf | cdr.bm44 <cf | plaq.bm44 <cf |bbscore.bm44 <cf  ) # the results in Table S2


# regulation loss burden
for(tis in c("BM10", "BM22","BM36","BM44")){
	pas.bm=patients[patients %in% row.names(subset(clin, tissue==tis))]
	load(paste("report_",tis, ".rda", sep="")) # the sub report for each tissue
	load(paste("output_",tis, ".rda", sep="")) # bicluster results
	tp=as.vector(sub.report$type)
	used.ones=as.vector(sub.report$TF)
	cc.reg.pa=sapply(used.ones, function(tf){ #transform the the patient's regulation status into 1 or 0
		op=rep(0, length(pas.bm))
    	nms=as.vector(as.matrix(ann[genes[tf], 2]))
    	pas<-res[[as.character(tf)]][["pas"]]
    	nges<-res[[as.character(tf)]][["nges"]]
    	brk=which.max(nges)
		bk=which(nges > 5 )[1]
    	select.pas=pas.bm[-1 * pas[1:brk]]
		non.select.pas=pas.bm[!pas.bm %in% select.pas]
    	if(is.na(bk) | bk < 10){
			return(op)
		}
		op[pas.bm%in% select.pas]=1
		return(op);
    })
    row.names(cc.reg.pa) <- pas.bm
    colnames(cc.reg.pa) <- as.vector(as.matrix(ann[genes[used.ones], 2]))
    cof=sapply(tp, function(x) if(x == "MR" ) 0.7 else 1)
    rlb.bm=rowSums(sapply(1 : ncol(cc.reg.pa), function(i) cc.reg.pa[, i] * cof[i]))
}
	