# install.packages("data.table")
# install.packages("caTools")
# install.packages("gridExtra")

library(data.table)
library(ggplot2)
library(caTools)
library(gridExtra)

source("R/auxiliary_functions.R")

# import data
bulkStandard <- fread("raw-data/1101.tsv") # 1ug of bulk
bulkDP <- fread("raw-data/1110.tsv") # 120pg of bulk
sampleDP <- fread("raw-data/1112.tsv") # 10-20 cells on HiSeq

# process data
bulkStandard <- summarise.data(bulkStandard)
bulkDP <- summarise.data(bulkDP)
sampleDP <- summarise.data(sampleDP)

### sanity check data
# no SNPs should have more than one base with numbers in
bulkDP[,zero.allele:=sum(a==0,c==0,t==0,g==0),by=c("chr","pos")]
bulkDP[zero.allele!=3 & zero.allele!=4]

# how many snps with zero reads in each sample
nrow(bulkStandard[ref==0 & a==0 & c==0 & t==0 & g==0])
nrow(bulkDP[ref==0 & a==0 & c==0 & t==0 & g==0])
nrow(sampleDP[ref==0 & a==0 & c==0 & t==0 & g==0])

# overview distribution plots
hist(bulkStandard$total,breaks=1000,xlim=c(0,100))
hist(bulkDP$total,breaks=1000,xlim=c(0,100))
hist(sampleDP$total,breaks=1000,xlim=c(0,100))

hist(bulkStandard$bin.median,breaks=1000,xlim=c(0,100))
hist(bulkDP$bin.median,breaks=1000,xlim=c(0,100))

hist(bulkStandard$bin.mean,breaks=1000,xlim=c(0,100))
hist(bulkDP$bin.mean,breaks=1000,xlim=c(0,100))

hist(bulkStandard$BAF,breaks=1000,ylim=c(0,70000))
hist(bulkDP$BAF,breaks=1000,ylim=c(0,70000))
hist(sampleDP$BAF,breaks=1000,ylim=c(0,70000))


# plot binned median over points for WGS and DigiPico
grid.arrange(grobs=list(plot.genomic.loadings(bulkStandard,"Bulk tumour, Standard WGS"),
                        plot.genomic.loadings(bulkDP,"Bulk tumour, DigiPico Sequencing")),
             layout_matrix=rbind(c(1),c(2)))


####
#### changepoint
####
# install.packages("changepoint")
library("changepoint")

cptST <- cpt.meanvar(bulkStandard[chr==6]$total,test.stat="Normal",method="BinSeg",Q=30,minseglen=3000)
cptDP <- cpt.meanvar(bulkDP[chr==6]$total,test.stat="Poisson",method="BinSeg",Q=30,minseglen=3000)

# PELT versions
# cptST <- cpt.meanvar(bulkStandard[chr==6]$total,test.stat="Normal",method="PELT",penalty="Manual",pen.value=10000)
# cptDP <- cpt.meanvar(bulkDP[chr==6]$total,test.stat="Poisson",method="PELT",penalty="Manual",pen.value=8000)

# pen.value=c(2*log(length(bulkDP[chr==6]$total)),100*log(length(bulkDP[chr==6]$total))),penalty="CROPS"


# convert segment snp number to genomic coordinates
BPs.bulkDP <- bulkDP[chr==6][cpts(cptDP)]$pos
BPs.bulkStandard <- bulkStandard[chr==6][cpts(cptST)]$pos

# calculate start and end points of segments in genomic co-ordinates
segments.bulkDP <- data.frame(start=c(1,BPs.bulkDP),end=c(BPs.bulkDP-1,max(bulkDP[chr==6]$pos)))
segments.bulkStandard <- data.frame(start=c(1,BPs.bulkStandard),end=c(BPs.bulkStandard-1,max(bulkDP[chr==6]$pos)))

# annotate SNPs with segment number
bulkDP$segment.no <- NULL
for (i in 1:nrow(segments.bulkDP)){
  bulkDP[chr==6 & pos>segments.bulkDP[i,]$start & pos<segments.bulkDP[i,]$end,segment.no:=i]
}

# bulkDP[chr==6,.N,by=segment.no]

for (i in 1:nrow(segments.bulkStandard)){
  bulkStandard[chr==6 & pos>segments.bulkStandard[i,]$start & pos<segments.bulkStandard[i,]$end,segment.no:=i]
}

# calculate mode total read count per segment
bulkDP$segment.mode <- NULL
bulkStandard$segment.mode <- NULL

bulkDP[chr==6 & !is.na(segment.no),segment.mode:=Mode(total),by=c("chr","segment.no")]
bulkStandard[chr==6 & !is.na(segment.no),segment.mode:=Mode(total),by=c("chr","segment.no")]

# one NA inbetween each segment so vertical lines arn't plotted
bulkDP[chr==6 & is.na(segment.mode)]
bulkStandard[chr==6 & is.na(segment.mode)]

#####

# common mu's from mclust on whole genome (+0.66 &1.33 at start)
initial.parameters.realistic <- data.table(
  w=c(0.5),
  mu=c(2,2.61,4.81,7.72,10.44,13.00,17.17,23.84,38.32,50,60,70,300,0.66,1.33),
  sigma2=c(8^2),
  component.type=c(rep("common",9),rep("specific",5)))


### fit actual data to model
chr6.input <- bulkDP[chr==6 & segment.no %in% c(1:31),.(vals=total,seg=as.character(segment.no))]
chr6.output <- fit.model(chr6.input,rho=0.5,initial.parameters.realistic)
plot(chr6.output$likelihood[2:length(chr6.output$likelihood)])
chr6.input$comp <- "A"

grid.arrange(grobs=list(plot.components(chr6.input[seg %in% c(5,15,22,31)],chr6.output,segment.subset = c("5","15","22","31")),
                        plot.components(chr6.input[seg %in% c(5,15,22,31)],chr6.output,segment.subset = c("5","15","22","31"),type="density"),
                        plot.components(chr6.input[seg %in% c(5,15,22,31)],chr6.output,segment.subset = c("5","15","22","31"),type="QQ")),
             layout_matrix=rbind(c(1),c(2),c(3)))

######

# compare rho to mode per segment and segment modes of normal vs digipico sequencing
grid.arrange(grobs=list(plot.rho(bulkDP[chr==6],parse.output.parameters(chr6.output)),
                        plot.breakpoints(bulkDP[chr==6],BPs.bulkDP,"BinSeg - DP",Mode(bulkDP[chr==6]$total)),
                        plot.breakpoints(bulkStandard[chr==6],BPs.bulkStandard,"BinSeg - standard",Mode(bulkStandard[chr==6]$total))),
             layout_matrix=rbind(c(1),c(2),c(3)))

# save breakpoints and model parameters in human + machine readable files
sink("output-human.txt")
chr6.output$breakpoints <- BPs.bulkDP
print(chr6.output)
sink()
saveRDS(chr6.output,"output.rds")

#### Check distribution of per segment read count

grid.arrange(grobs=list(plot.histogram.stat(0,18568208,"Normal"),
                        plot.histogram.stat(64622787,133557441,"Normal"),
                        plot.histogram.stat(165193244,1165193244,"Deleted"),
                        plot.histogram.stat(47653485,57239589,"Amplified")),
             layout_matrix=rbind(c(1,2),c(3,4)))

# QQ Plots
#install.packages("car")
library(car)
qqp(bulkDP[chr==6 & pos>lower.bound & pos<upper.bound]$total, "nbinom",size=2.67,mu=16.2)
qqp(bulkDP[chr==6 & pos>lower.bound & pos<upper.bound]$total, "pois",lambda=16.2)


