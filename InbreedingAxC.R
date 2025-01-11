#ijob -A berglandlab -c20 -p standard --mem=40G
#module load gcc/11.4.0  openmpi/4.1.4 R/4.3.1; R

#Determine alleles that are depressed in different regions?
#Putative regions that code for inbreeding depression?


library(data.table)
library(foreach)
library(SeqArray)
library(doMC)
  registerDoMC(20)
### working dir
  setwd("/project/berglandlab/Robert")

### open
  genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/MapJune2020_ann.seq.gds")



### filtered SNP set
  snp.dt <- fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/finalsetsnpset01pulex_table_20200623")

### sample metadata
  samps <- fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")

### get heterozygous sites in C parents
  seqSetFilter(genofile, sample.id=samps[AxCF1Hybrid=="1"]$clone, variant.id=snp.dt$variant.ids)

unique(samps$SC)
### get genotype tables

scsubids <- seqGetData(genofile, "sample.id")

samps2 <- as.data.table(scsubids)

setnames(samps2, "scsubids", "clone")


  seqSetFilter(genofile, sample.id=samps[AxCF1Hybrid=="1"]$clone, variant.id=snp.dt$variant.ids)

unique(samps$SC)
### get genotype tables
  genos <- seqGetData(genofile, "$dosage")
  clone <- seqGetData(genofile, "sample.id")

  # genos.tab <- apply(genos, 2, function(x) mean(x==1, na.rm=T))
  # genos.tab0 <- apply(genos, 2, function(x) mean(x==0, na.rm=T))
  # genos.tab2 <- apply(genos, 2, function(x) mean(x==2, na.rm=T))

  table(round(genos.tab, 1))

### heterozygous sites
  seqResetFilter(genofile)
  seqSetFilter(genofile, variant.id=snp.dt$variant.ids, sample.id=samps[AxCF1Hybrid=="1"]$clone)
  genos <- seqGetData(genofile, "$dosage")
  table(genos)


  # seqResetFilter(genofile)
  # seqSetFilter(genofile, variant.id=snp.dt[which(genos.tab>.9)]$variant.ids, sample.id=samps[SC=="A"]$clone)
  # genos <- seqGetData(genofile, "$dosage")
  # Csites <- table(genos)
  # save(Csites, file="/project/berglandlab/Robert/Csites.Rdata")



  selfedC_hetsites <- data.table(genotype=expand.grid(genos)$Var1,
             clone=rep(seqGetData(genofile, "sample.id"), times=dim(genos)[2]),
             variant.id=rep(seqGetData(genofile, "variant.id"), each=dim(genos)[1]))

selfedC_hetsites_byclone <- selfedC_hetsites[, .N, by=list(clone, genotype)]


selfedC_hetsites2 <- na.omit(selfedC_hetsites2)

o <- foreach(snp.i=snp.dt$variant.ids[1:100], .errorhandling="remove")%do%{
    # snp.i <- snp.dt[freq>.05 & freq<.95 ]$variantId[1:100]
  seqSetFilter(genofile, variant.id=snp.dt$variant.ids, sample.id=samps[AxCF1Hybrid=="1"]$clone)

    ### iterate through comparisons
      foreach(set.i=unique(samps2$clone), .combine="rbind", .errorhandling="remove")%do%{
        #set.i <- April17_2018_D8_Male1
        #tmp.tmp <- snp.tmp[Lifestage%in%sets[set==set.i]$V2]

# selfedC_hetsites2 <- subset(selfedC_hetsites, variant.id ==76 | variant.id == 3719894)

          t1<- glm(genotype ~ clone, data=selfedC_hetsites2, family=poisson())
          t0 <- glm(genotype~1, data=selfedC_hetsites2, family=poisson())

          modOut <- anova(t1, t0, test="Chisq")

          ret <- data.table(variant.id=snp.i,
                     p.value=c(modOut$Pr[2]))

      }
     return(ret)
    #snp.tmp[,male:=Lifestage=="Male"]

  }
  o <- rbindlist(o)








  genos.tab <- apply(genos, 2, function(x) mean(x==1, na.rm=T))

  table(round(genos.tab, 1))

# ### heterozygous sites
#   seqResetFilter(genofile)
#   seqSetFilter(genofile, variant.id=snp.dt[which(genos.tab>.9)]$variant.ids, sample.id=samps[SC=="C"]$clone)
#   genos <- seqGetData(genofile, "$dosage")
#   table(genos)


  # seqResetFilter(genofile)
  # seqSetFilter(genofile, variant.id=snp.dt[which(genos.tab>.9)]$variant.ids, sample.id=samps[SC=="A"]$clone)
  # genos <- seqGetData(genofile, "$dosage")
  # Csites <- table(genos)
  # save(Csites, file="/project/berglandlab/Robert/Csites.Rdata")



#   selfedC_hetsites <- data.table(genotype=expand.grid(genos)$Var1,
#              clone=rep(seqGetData(genofile, "sample.id"), times=dim(genos)[2]),
#              variant.id=rep(seqGetData(genofile, "variant.id"), each=dim(genos)[1]))

# selfedC_hetsites_byclone <- selfedC_hetsites[, .N, by=list(clone, genotype)]

#   save(selfedA_hetsites_byclone, file="/project/berglandlab/Robert/selfedA_hetsites_byclone.Rdata")


#   out <- merge(selfedC_hetsites, snp.dt, by.x="variant.id", by.y="variant.ids")


### save
  save(out, file="/project/berglandlab/Robert/selfed_c.Rdata")


### dl
#  scp aob2x@rivanna.hpc.virginia.edu:~/selfed_c.Rdata ~/.


### load
  library(ggplot2)
  library(data.table)
  library(ggtern)

  load("/project/berglandlab/Robert/selfed_c.Rdata")

  out.ag <- out[,list(fAA=mean(genotype==2, na.rm=T), fAa=mean(genotype==1, na.rm=T), faa=mean(genotype==0, na.rm=T), N=sum(!is.na(genotype))),
                 list(variant.id)]


  table(out.ag$faa==0)

   ggtern(data=out.ag[N>20], aes(fAA,fAa,faa)) + geom_point()

   base + geom_point() + geom_jitter()

### another approach

  out.ag2 <- out[,list(fAA=mean(genotype==2, na.rm=T), fAa=mean(genotype==1, na.rm=T), faa=mean(genotype==0, na.rm=T), N=sum(!is.na(genotype))),
                 list(clone)]

  out.ag2.long <- melt(out.ag2[,-"N",with=F], id.vars="clone")


  ggplot(data=out.ag2.long, aes(y=value, x=variable, color=variable)) + geom_boxplot()

  prop.table(table(out$genotype))










library(gdsfmt)
library(SNPRelate)
library(data.table)
library(ggplot2)
library(foreach)
library(lattice)
library(tidyr)
library(SeqArray)
library(dplyr)
library(tidyverse)
### working dir

### open

#genofile <- seqOpen("MapJune2020_ann.seq.gds")

### filtered SNP set

load("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/finalsetsnpset01_20200623.Rdata")

### sample metadata

samps <- fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")

AxCF1 <- subset(samps, AxCF1Hybrid == 1)

table(AxCF1$SC)

sampleId <- seqGetData(genofile, "sample.id")
scsubids <- AxCF1$clone

sc <- AxCF1$SC
seqSetFilter(genofile, sample.id=scsubids)

scsubids <- seqGetData(genofile, "sample.id")


### set some global parameters
maf <- 0.001
missing.rate <- 0.15
threads <- 10

ibs <- snpgdsIBS(genofile, snp.id=finalsetsnpset01, sample.id=scsubids, num.thread=20, maf=maf,
                 missing.rate=0.15, autosome.only = FALSE)
#Working space: 76 samples, 29,811 SNVs

### a bit of re-formating of the ibs matrix
ibs.mat <- ibs$ibs
rownames(ibs.mat) <- ibs$sample.id
colnames(ibs.mat) <- ibs$sample.id

ibs.matdt <- as.data.table(ibs.mat)
ibs.matdt$cloneA <- scsubids
ibs.long<- melt(ibs.matdt, measure.vars=scsubids, variable.name="cloneB", value.name="IBS")

### remake ibs.long
ibs.matdt <- as.data.table(ibs.mat)
ibs.matdt$cloneA <- scsubids
ibs.long<- melt(ibs.matdt, measure.vars=scsubids, variable.name="cloneB", value.name="IBS")
ibs.long <- na.omit(ibs.long)
# 





### first, need to tack in SC identities
AxCF1

sctmp <- AxCF1[,c("clone","SC")]

setnames(ibs.long, "cloneA", "clone")
setkey(ibs.long, "clone")
setkey(sctmp, "clone")
ibs.long <- merge(ibs.long, sctmp)
setnames(ibs.long, "clone", "cloneA")
setnames(ibs.long, "SC", "SC.A")

setnames(ibs.long, "cloneB", "clone")
setkey(ibs.long, "clone")
setkey(sctmp, "clone")
ibs.long <- merge(ibs.long, sctmp)
setnames(ibs.long, "clone", "cloneB")
setnames(ibs.long, "SC", "SC.B")
# 
# ibs.long <- ibs.long[Species.x=="pulex" & Species.y=="pulex"]
# ibs.long <- ibs.long[population.x!="Dramp" & population.y!="Dramp" & population.x!="DLily" &
#                        population.y!="DLily" & population.x!="DOil" & population.y!="DOil" & population.x!="DMud" &
#                        population.y!="DMud"]

### re-rank SCs based on size in pond
# 
ibs.long <- ibs.long[,c("cloneA", "cloneB", "SC.A", "SC.B", "IBS"), with=F]


ibs.long.ag <- ibs.long[,list(clone.n = length(IBS)), list(SC.A, cloneA) ]
ibs.long.ag$clone.n<- ifelse(ibs.long.ag$SC.A=="OO", 1, ibs.long.ag$clone.n)
setkey(ibs.long.ag, SC.A, clone.n)


# ibs.long.ag <- ibs.long[,list(pond.n = length(IBS)), list(pondA, SC.A) ]
# ibs.long.ag$pond.n<- ifelse(ibs.long.ag$SC.A=="OO", 1, ibs.long.ag$pond.n)
# setkey(ibs.long.ag, pondA, pond.n)
ibs.long.ag[,clone.sc.rank := ibs.long.ag[,list(clone.sc.rank = rank(-clone.n, ties="random")), list(SC.A)]$clone.sc.rank]
ibs.long.ag[,clone.sc.rank := letters[clone.sc.rank]]

### be lazy and write a loop
ibs.long.2 <- foreach(i=1:dim(ibs.long.ag)[1], .combine="rbind")%do%{
  temp <- ibs.long[cloneA==ibs.long.ag$cloneA[i] & SC.A==ibs.long.ag$SC.A[i]]
  temp[,sc.a:=ibs.long.ag$clone.sc.rank[i]]
  temp
}

ibs.long.3 <- foreach(i=1:dim(ibs.long.ag)[1], .combine="rbind")%do%{
  temp <- ibs.long.2[cloneB==ibs.long.ag$cloneA[i] & SC.B==ibs.long.ag$SC.A[i]]
  temp[,sc.b:=ibs.long.ag$clone.sc.rank[i]]
  temp
}

ibs.long <- ibs.long.3

setkey(ibs.long, cloneA, cloneB, sc.a, sc.b)

### next, generate [s]uper[c]lone[i]ds for individual 'A' and 'B'
# 
ibs.long[,scid.a := paste(sc.a, sprintf("%03d", as.numeric(as.factor(cloneA))), sep=".")]
ibs.long[,scid.b := paste(sc.b, sprintf("%03d", as.numeric(as.factor(cloneB))), sep=".")]

### group on pond
# # 
# ibs.long[,pondA := factor(pondA, levels=c("D10", "DCat", "D8", "DBunk"))]
# ibs.long[,pondB := factor(pondB, levels=c("D10", "DCat", "D8", "DBunk"))]
# 
# ibs.long[,scid.a := paste(LETTERS[as.numeric(SC.A)], scid.a, sep=".")]
# ibs.long[,scid.b := paste(LETTERS[as.numeric(SC.B)], scid.b, sep=".")]
# 
# ibs.long[,scid.a := as.factor(scid.a)]
# ibs.long[,scid.b := as.factor(scid.b)]

### tack in buffer cloneIds for graphical purposes
# ibs.long <- rbind(ibs.long,
#                   data.table(scid.a=c(paste(unique(paste(LETTERS[as.numeric(ibs.long[pondA%in%c("D10", "D8", "DBunk",
#                                                                                                 "DCat")]$pondA)], min(ibs.long$SC.A), sep=".")), c("000"), sep="."),
#                                       paste(unique(paste(LETTERS[as.numeric(ibs.long[pondA%in%c("D10", "D8", "DBunk",
#                                                                                                 "DCat")]$pondA)], max(ibs.long$SC.A), sep=".")), c("999"), sep=".")),
#                              scid.b=c(paste(unique(paste(LETTERS[as.numeric(ibs.long[pondB%in%c("D10", "D8", "DBunk",
#                                                                                                 "DCat")]$pondB)], min(ibs.long$SC.B), sep=".")), c("000"), sep="."),
#                                       paste(unique(paste(LETTERS[as.numeric(ibs.long[pondB%in%c("D10", "D8", "DBunk",
#                                                                                                 "DCat")]$pondB)], max(ibs.long$SC.B), sep=".")), c("999"), sep="."))),
#                   fill=T)
# ibs.long[,scid.a := factor(scid.a, levels=sort(unique(as.character(scid.a))))]
# ibs.long[,scid.b := factor(scid.b, levels=sort(unique(as.character(scid.b))))]
# 
# ### make lower triangle poofy-de-poof
# ibs.long[,dist.noTri := IBS]
# ibs.long[as.numeric(scid.a)>as.numeric(scid.b), dist.noTri:=NA]
# 
# ### make pond bounding boxes
# ibs.long.ag <- data.table(scid.a.min=paste(unique(paste(LETTERS[as.numeric(ibs.long[pondA%in%c("D10", "D8", "DBunk", "DCat")]$pondA)], min(ibs.long$SC.A, na.rm=T), sep=".")), c("000"), sep="."),
#                           scid.a.max=paste(unique(paste(LETTERS[as.numeric(ibs.long[pondA%in%c("D10", "D8", "DBunk", "DCat")]$pondA)], max(ibs.long$SC.A, na.rm=T), sep=".")), c("999"), sep="."),
#                           scid.b.min=paste(unique(paste(LETTERS[as.numeric(ibs.long[pondB%in%c("D10", "D8", "DBunk", "DCat")]$pondB)], min(ibs.long$SC.B, na.rm=T), sep=".")), c("000"), sep="."),
#                           scid.b.max=paste(unique(paste(LETTERS[as.numeric(ibs.long[pondB%in%c("D10", "D8", "DBunk", "DCat")]$pondB)], max(ibs.long$SC.B, na.rm=T), sep=".")), c("999"), sep="."))
# 
# ibs.long.ag[,scid.a.min := as.numeric(factor(scid.a.min, levels=sort(unique(as.character(ibs.long$scid.a)))))]
# ibs.long.ag[,scid.a.max := as.numeric(factor(scid.a.max, levels=sort(unique(as.character(ibs.long$scid.a)))))]
# ibs.long.ag[,scid.b.min := as.numeric(factor(scid.b.min, levels=sort(unique(as.character(ibs.long$scid.b)))))]
# ibs.long.ag[,scid.b.max := as.numeric(factor(scid.b.max, levels=sort(unique(as.character(ibs.long$scid.b)))))]


library(ggplot2)
library(viridis)

load("ibs.long.Rdata")

### plot it
h.just <- .25
v.just <- .25
l.size <- 1.5


ibs.long2111 <- ibs.long %>% group_by(cloneA) %>% ungroup()

ibs.long2111 <- ibs.long %>% group_by(cloneA) %>% mutate(ibs_mean = mean(IBS)) %>% 
  arrange(ibs_mean) %>% ungroup()

ibs.long2111 <- as.data.frame(ibs.long2111)




#Next, if two numbes are above a certain amount (0.98?) group together with a different column
#filter out selfed C's, see if there are any chromosomes/regions that are deprecated


corrmatrix <- ggplot(data=subset(ibs.long, IBS >= .50), aes(cloneA, cloneB, fill=IBS)) +
  geom_raster() +
  scale_fill_viridis(option="D")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
corrmatrix

ggplot(data=ibs.long, aes(x=IBS)) +
  geom_histogram(bins=100)+
  ylab("Number of pairwise comb.")+ 
  xlab("IBS") + theme_bw()


corrmatrix <- ggplot(data=ibs.long, aes(scid.a, scid.b, fill=IBS)) +
  geom_raster() +
  scale_fill_viridis(option="D")
corrmatrix


ggsave(corrmatrix, file="corrmatrix_RP.pdf")














### libraries
  library(SeqArray)
  library(data.table)

### working dir
### filtered SNP set


genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/MapJune2020_ann.seq.gds")

genofile
##seqResetFilter(genofile)

### filtered SNP set
  snp.dt <- fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/finalsetsnpset01pulex_table_20200623")

### sample metadata
  samps <- fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")


### load sampleMeta data

# samps[,sampleId:=paste("Sample", SampleID, sep="")]

### get sample names
sampleId <- scsubids

AxCF1

### sample metadata
### get heterozygous sites in C parents
  seqSetFilter(genofile, sample.id=samps[SC=="AxCF1"]$clone, variant.id=snp.dt$variant.ids)

unique(samps$SC)
### get genotype tables
  genos <- seqGetData(genofile, "$dosage")

  genos.tab <- apply(genos, 2, function(x) mean(x==1, na.rm=T))

  table(round(genos.tab, 1))

### heterozygous sites
  seqResetFilter(genofile)
  seqSetFilter(genofile, variant.id=snp.dt[which(genos.tab>.9)]$variant.ids, sample.id=samps[SC=="C"]$clone)
  genos <- seqGetData(genofile, "$dosage")
  table(genos)


  seqResetFilter(genofile)
  seqSetFilter(genofile, variant.id=snp.dt[which(genos.tab>.9)]$variant.ids, sample.id=samps[SC=="A"]$clone)
  genos <- seqGetData(genofile, "$dosage")
  Csites <- table(genos)
  save(Csites, file="/project/berglandlab/Robert/Csites.Rdata")



  selfedC_hetsites <- data.table(genotype=expand.grid(genos)$Var1,
             clone=rep(seqGetData(genofile, "sample.id"), times=dim(genos)[2]),
             variant.id=rep(seqGetData(genofile, "variant.id"), each=dim(genos)[1]))

selfedC_hetsites_byclone <- selfedC_hetsites[, .N, by=list(clone, genotype)]

  save(selfedA_hetsites_byclone, file="/project/berglandlab/Robert/selfedA_hetsites_byclone.Rdata")


  out <- merge(selfedC_hetsites, snp.dt, by.x="variant.id", by.y="variant.ids")


### save
  save(out, file="~/selfed_c.Rdata")


### dl
#  scp aob2x@rivanna.hpc.virginia.edu:~/selfed_c.Rdata ~/.


### load
  library(ggplot2)
  library(data.table)
  library(ggtern)

  load("~/selfed_c.Rdata")

  out.ag <- out[,list(fAA=mean(genotype==2, na.rm=T), fAa=mean(genotype==1, na.rm=T), faa=mean(genotype==0, na.rm=T), N=sum(!is.na(genotype))),
                 list(variant.id)]


  table(out.ag$faa==0)

   ggtern(data=out.ag[N>20], aes(fAA,fAa,faa)) + geom_point()

   base + geom_point() + geom_jitter()

### another approach

  out.ag2 <- out[,list(fAA=mean(genotype==2, na.rm=T), fAa=mean(genotype==1, na.rm=T), faa=mean(genotype==0, na.rm=T), N=sum(!is.na(genotype))),
                 list(clone)]

  out.ag2.long <- melt(out.ag2[,-"N",with=F], id.vars="clone")


  ggplot(data=out.ag2.long, aes(y=value, x=variable, color=variable)) + geom_boxplot()

  prop.table(table(out$genotype))

















# ### libraries
# library(data.table)
# library(foreach)
# library(SeqArray)
# library(doMC)
#   registerDoMC(20)

# ### genofile

#   genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/MapJune2020_ann.seq.gds")

# genofile
# ##seqResetFilter(genofile)

# ### load sampleMeta data
# #samps <- fread("/project/berglandlab/Robert/DorsetPooledSequencing2018_2019/AllPooledVCFsIncludingDcatDbunk/DorsetpoolsALL_metadata.csv")

# # samps[,sampleId:=paste("Sample", SampleID, sep="")]

# ### get sample names
# sampleId <- scsubids

# ### iterate through each sample to get missing rate & coverage
# sampleStats <- foreach(sample.i=scsubids, .combine="rbind") %dopar%{
#   ## sample.i = sampleId[1]
#   ### subsample to sample.i
#   seqSetFilter(genofile, sample.id=sample.i)
  
#   ### get missing data rate
#   missingRate <- mean(is.na(seqGetData(genofile, "$dosage")))
  
#   ### coverage
#   depth <- seqGetData(genofile, "annotation/format/DP")
  
#   ### output
#   data.table(sampleId=sample.i,
#              missingRate=missingRate,
#              meanDepth=mean(depth, na.rm=T),
#              medianDepth=median(depth, na.rm=T),
#              upper_quantile=quantile(depth, .975, na.rm=T),
#              lower_quantile=quantile(depth, .025, na.rm=T))
# }

# sampleStats

# samps <- as.data.table(scsubids)
# setnames(samps, "scsubids", "sampleId")

# samps <- merge(samps, sampleStats, by="sampleId")

# ### create a snp.dt object
# ### get basic info
# snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
#                      pos=seqGetData(genofile, "position"),
#                      nAlleles=seqGetData(genofile, "$num_allele"),
#                      variantId=seqGetData(genofile, "variant.id"))


# ### iterate through sites to identify outliers
# snpStats <- foreach(samp.i=sampleStats$sampleId)%dopar%{
#   # samp.i=sampleStats$sampleId[1]
#   message(samp.i)
#   seqSetFilter(genofile, sample.id=samp.i, variant.id=snp.dt[nAlleles==2]$variantId)
  



#   #breaks after this point, difference between pooled seq and regular?

#   ### coverage
#   snp.tmp <- data.table(variantId=seqGetData(genofile, "variant.id"),
#                         dp=seqGetData(genofile,"annotation/format/DP")[1,],
#                         ad=seqGetData(genofile,"annotation/format/AD")[1,],
#                         position=seqGetData(genofile,"position"),
#                         sampleId=samp.i)
  
  
#   ### add in quantiles
#   snp.tmp[,use:=(dp>5 & dp<=quantile(dp, .95, na.rm=T) & !is.na(dp))]
#   snp.tmp[,freq:=ad/dp]
  
# }


# snpStats <- rbindlist(snpStats)




# snpStats_chr <- merge(snp.dt, snpStats, by="variantId")


# #write.csv(snpStats_chr, file="/project/berglandlab/Robert/DorsetPooledSequencing2018_2019/snpStats_chr.csv")


# #rerun with USE = TRUE in frequency ,  freq=mean(freq[use=t]

# snpStats.ag <- snpStats[,list(nUse=sum(use), allUse=all(use), freq=mean(freq, na.rm=T)), variantId]
# table(snpStats.ag$allUse)

# snp.dt <- merge(snp.dt, snpStats.ag, by="variantId")

# snp.dt[allUse==T & freq>.05 & freq<.95]


# ### save
# ### save(snp.dt, file="/project/berglandlab/Robert/DorsetPooledSequencing2018_2019/DorsetPooled_2018_2019_snpdt.csv")

# write.csv(snp.dt,"/project/berglandlab/Robert/DorsetPooledSequencing2018_2019/AllPooledVCFsIncludingDcatDbunk/DorsetpoolsALL_snpdt.csv")
# write.csv(snpStats,"/project/berglandlab/Robert/DorsetPooledSequencing2018_2019/AllPooledVCFsIncludingDcatDbunk/DorsetpoolsALL_snpStats.csv")


# usetrue <- subset(snp.dt, use=="TRUE")

# snp.dt