# install.packages("data.table")
# install.packages("caTools")
# install.packages("gridExtra")

library(data.table)
library(ggplot2)
library(caTools)
library(gridExtra)

summarise.data <- function(data){
  data[,non.ref:=sum(a,c,t,g),by=c("chr","pos")]
  data[,total:=sum(ref,a,c,t,g),by=c("chr","pos")]
  data[(non.ref+ref)>10,BAF.thresh:=non.ref/(non.ref+ref),by=c("chr","pos")]
  data[,BAF:=non.ref/(non.ref+ref),by=c("chr","pos")] # B allele freq.
  data[,running.mean:=runmean(total, k=500, alg="fast", endrule="keep", align = "center")]
  data[,running.median:=runmed(total, k=501)] 
  number.of.bins=20000 # 1000
  bins <- rep(1:number.of.bins, each = ceiling(nrow(data)/number.of.bins)) # variable genomic width, equal data point bins
  bins <- bins[1:nrow(data)] # crop bins down to fit
  data[,bin:=bins]
  data[,bin.median:=median(as.numeric(total)),by=bin]
  data[,bin.mean:=mean(as.numeric(total)),by=bin]
}

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


# set chromosome lengths for calculating plotting position
# http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/data/ GRCh38.p7
chromosome.lengths <- data.table(chr=factor(c(1:23)),chr_length=c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468,156040895)+5000000)
chromosome.lengths[,genomic_offset:=cumsum(as.numeric(chr_length))-chr_length]
setkey(chromosome.lengths,chr)

plot.genomic.loadings <- function(temp,title){
  
  # Correct / create chromosome positions / genomic locations
  temp$chr <- factor(temp$chr, levels = c(1:23))
  
  setkey(temp,chr)
  temp <- chromosome.lengths[temp]
  temp[,genomic_position:=genomic_offset+pos]
  
  #if(autozoom==TRUE){
    temp <- temp[chr=="6"]
  #}
    
    P <- ggplot(temp, aes(genomic_position,total)) +
    ylim(0,100) +
      geom_point(size=0.5,alpha=0.05,aes(color=chr)) +
      scale_colour_manual(values =rep_len(c("black", "skyblue4"),23)) +
      geom_line(aes(genomic_position,bin.median),colour="red")+
    ggtitle(title) +
    ylab("Total Read Count") +
    xlab("Genomic Coordinate")


#   BAF plot
#   P <- ggplot(temp, aes(genomic_position,BAF)) +
#     geom_point(size=0.5,alpha=0.01,aes(color=chr)) +
#     scale_colour_manual(values =rep_len(c("black", "skyblue4"),23))
}

summary(bulkStandard)

# plot binned median over points for WGS and DigiPico
grid.arrange(grobs=list(plot.genomic.loadings(bulkStandard,"Bulk tumour, Standard WGS"),
                        plot.genomic.loadings(bulkDP,"Bulk tumour, DigiPico Sequencing")),
             layout_matrix=rbind(c(1),c(2)))

####
#### DNAcopy
####
# source("https://bioconductor.org/biocLite.R")
# biocLite("DNAcopy")
library(DNAcopy)
CNA.object <- CNA( genomdat = log2(bulkStandard$total/29), chrom = bulkStandard$chr, maploc = bulkStandard$pos, data.type = 'logratio')
CNA.smoothed <- smooth.CNA(CNA.object)
subset.CNA.object <- subset(CNA.smoothed,chromlist=6)
segs <- segment(subset.CNA.object, verbose=0, min.width=5)

plot(segs, plot.type="p")
plot(segs, plot.type="w")


####
#### changepoint
####
# install.packages("changepoint")
library("changepoint")

cptST <- cpt.meanvar(bulkStandard[chr==6]$total,test.stat="Normal",method="BinSeg",Q=20,minseglen=1500)
cptDP <- cpt.meanvar(bulkDP[chr==6]$total,test.stat="Poisson",method="BinSeg",Q=20,minseglen=1500)

# PELT versions
# cptST <- cpt.meanvar(bulkStandard[chr==6]$total,test.stat="Normal",method="PELT",penalty="Manual",pen.value=10000)
# cptDP <- cpt.meanvar(bulkDP[chr==6]$total,test.stat="Poisson",method="PELT",penalty="Manual",pen.value=8000)

# pen.value=c(2*log(length(bulkDP[chr==6]$total)),100*log(length(bulkDP[chr==6]$total))),penalty="CROPS"

# main plotting function
plot.breakpoints <- function(data,breakpoints,title="",average){
 P <- ggplot(data[chr==6], aes(pos,total))+
    ylim(0,100) +
    geom_point(size=0.5,alpha=0.05) +
    ggtitle(title) +
    ylab("Total Read Count") +
    xlab("Genomic Coordinate") +
    geom_line(aes(pos,segment.mode),colour="red",size=1.5) +
    geom_hline(yintercept=average,colour="blue",size=4,alpha=0.4)
 
  #geom_vline(xintercept=breakpoints,colour="red")
  #geom_line(aes(pos,bin.median),colour="blue",size=0.1) +
}

# plot rho over geonomic segments by mixed model output
plot.rho <- function(data,model.output){
  model.output$segment.no <- as.integer(output.chr6$segment)
  setkey(model.output,segment.no)
  setkey(data,segment.no)
  data <- unique(model.output[,.(rho,segment.no)])[data[,.(pos,segment.no)]]
  
  P <- ggplot(data, aes(pos,rho))+
    ylim(0,1) +
  ylab("Rho") +
  xlab("Genomic Coordinate") +
  geom_line(colour="red",size=1.5)
}

# convert segment snp number to genomic coordinates
BPs.bulkDP <- bulkDP[chr==6][cpts(cptDP)]$pos
BPs.bulkStandard <- bulkStandard[chr==6][cpts(cptST)]$pos

# calculate start and end points of segments in genomic co-ordinates
segments.bulkDP <- data.frame(start=c(1,BPs.bulkDP),end=c(BPs.bulkDP-1,max(bulkDP[chr==6]$pos)))
segments.bulkStandard <- data.frame(start=c(1,BPs.bulkStandard),end=c(BPs.bulkStandard-1,max(bulkDP[chr==6]$pos)))

# annotate SNPs with segment number
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

### fit actual data to model
chr6.input <- bulkDP[chr==6 & segment.no %in% c(1:21),.(vals=total,seg=as.character(segment.no))]
chr6.output <- fit.model(chr6.input,rho=0.5,initial.parameters.realistic)
plot(chr6.output$likelihood[2:length(chr6.output$likelihood)])
chr6.input$comp <- "A"
plot.components(chr6.input[seg==16],parse.output.parameters(chr6.output)[segment=="16"])


######


# compare segment modes of normal vs digipico sequencing
grid.arrange(grobs=list(plot.breakpoints(bulkStandard,BPs.bulkStandard,"BinSeg - standard",Mode(bulkStandard[chr==6]$total)),
                        plot.breakpoints(bulkDP,BPs.bulkDP,"BinSeg - DP",Mode(bulkDP[chr==6]$total))),
             layout_matrix=rbind(c(1),c(2)))

# compare rho to mode per segment
grid.arrange(grobs=list(plot.rho(bulkDP[chr==6],chr6.output),
                        plot.breakpoints(bulkDP,BPs.bulkDP,"BinSeg - DP",Mode(bulkDP[chr==6]$total))),
             layout_matrix=rbind(c(1),c(2)))

#### Check distribution of per segment read count

# x.wei<-rweibull(n=200,shape=2.1,scale=1.1) ## sampling from a Weibull distribution with parameters shape=2.1 and scale=1.1
# x.teo<-rweibull(n=200,shape=2, scale=1) ## theorical quantiles from a Weibull population with known paramters shape=2 e scale=1
# qqplot(x.teo,x.wei,main="QQ-plot distr. Weibull") ## QQ-plot
# abline(0,1) ## a 45-degree reference line is plotted

#install.packages("fitdistrplus")
#install.packages("gamlss")
# library(fitdistrplus)    # fits distributions using maximum likelihood
# library(gamlss)          # defines pdf, cdf of ZIP
# 
# plot.fitdist <- function(lower.bound,upper.bound){
# # FIT DISTRIBUTION (mu = mean of poisson, sigma = P(X = 0)
# fit_zip = fitdist(bulkDP[chr==6 & pos>lower.bound & pos<upper.bound]$total, 'nbinom', start = list(mu = 2, size = 0.5),lower = c(0, 0),upper=c(50,10))
# #fit_zip = fitdist(bulkDP[chr==6 & pos>0 & pos<18568208]$total, 'ZIP', start = list(mu = 2, sigma = 0.5),lower = c(0, 0),upper=c(50,1))
# #fit_zip = fitdist(bulkDP[chr==6 & pos>0 & pos<18568208]$total, 'pois')
# # VISUALIZE TEST AND COMPUTE GOODNESS OF FIT    
# plot(fit_zip)
# }
# plot.fitdist(0,18568208) # normal

# fit segments to negative binomial distributions
library(MASS)
plot.histogram.stat <- function(lower.bound,upper.bound,title){
paramsNB <- as.list(MASS::fitdistr(bulkDP[chr==6 & pos>lower.bound & pos<upper.bound]$total, "negative binomial"))

#paramsP <- as.list(MASS::fitdistr(bulkDP[chr==6 & pos>lower.bound & pos<upper.bound]$total, "Poisson"))$estimate
#paramsG <- as.list(MASS::fitdistr(bulkDP[chr==6 & pos>lower.bound & pos<upper.bound]$total, "gamma"))

P <- ggplot(bulkDP[chr==6 & pos>lower.bound & pos<upper.bound,.(total)], aes(total)) +
  geom_histogram(aes(y=..density..),binwidth=1) +
  xlim(0,80) +
  geom_point(aes(y=dnbinom(total,size=paramsNB$estimate['size'],mu=paramsNB$estimate['mu'])), colour="red",shape=18) +
  annotate("text", x = 60, y=0.02, label = paste("size:",signif(paramsNB$estimate['size'],3),"\n Mu: ",signif(paramsNB$estimate['mu'],3))) +
  ggtitle(title) +
  ylab("Density") +
  xlab("Total Read Count")

  # geom_point(aes(y=dpois(total,paramsP)), colour="red")
  # geom_point(aes(y=dpois(total,paramsG$estimate['shape'],paramsG$estimate['rate'])), colour="blue")
}

grid.arrange(grobs=list(plot.histogram.stat(0,18568208,"Normal"),
                        plot.histogram.stat(64622787,133557441,"Normal"),
                        plot.histogram.stat(165193244,1165193244,"Deleted"),
                        plot.histogram.stat(47653485,57239589,"Amplified")),
             layout_matrix=rbind(c(1,2),c(3,4)))

#install.packages("car")
library(car)
qqp(bulkDP[chr==6 & pos>lower.bound & pos<upper.bound]$total, "nbinom",size=2.67,mu=16.2)
qqp(bulkDP[chr==6 & pos>lower.bound & pos<upper.bound]$total, "pois",lambda=16.2)

#####

## Segment with negative binomail underlying distribution

install.packages("Segmentor3IsBack")
library(Segmentor3IsBack)
N=2000 
y=c(rnbinom(N,prob=0.3,size=0.15),rnbinom(2*N,prob=0.1,size=0.15),
    rnbinom(N,prob=0.6,size=0.15),compress=FALSE)
res2=Segmentor(bulkDP[chr==6]$total, model=3,Kmax=25) # 3=negative binomial

BPs.bulkDP.segmentor3 <- bulkDP[chr==6][res2@breaks[20,]]$pos

print(plot.breakpoints(bulkDP,BPs.bulkDP.segmentor3,"Segmentor3 - DP",Mode(bulkDP[chr==6]$total)))

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

####
#### ecp
####

#install.packages("ecp")
library("ecp")
data("ACGH")
acghData = ACGH$data
# temp <- e.agglo(acghData, member = 1:nrow(acghData), alpha = 2, penalty = function(cp){0})

mem = rep(seq(1,50),times=rep(180000/50,50)) # /10
# temp <- e.agglo(X=as.matrix(bulkStandard[chr==6]$total),mem=mem[1:length(bulkStandard[chr==6]$total)],penalty=function(cp) 0)

# too slow

####
#### Breakout Detection (Twitter)
####

# install.packages("devtools")
# devtools::install_github("twitter/BreakoutDetection")
library(devtools)
library(BreakoutDetection)
plot(unique(bulkDP[chr==6,.(bin.median,bin)])$bin.median)

bulkStandard[,bin.position:=median(pos),by=bin]
bulkDP[,bin.position:=median(pos),by=bin]

# find breakpoints using bins, return actual co-ordinate
breakout.detection <- function(data){
  BDout = breakout(unique(data[chr==6,.(bin.median,bin)])$bin.median, min.size=10, method='multi', beta=.001, degree=1, plot=TRUE)
  BDout.bins <- unique(data[chr==6,.(bin.median,bin)])[BDout$loc]$bin
  return(unique(data[bin %in% BDout.bins]$bin.position))
}

BD.standard <- breakout.detection(bulkStandard)
BD.DP <- breakout.detection(bulkDP)

grid.arrange(grobs=list(plot.breakpoints(bulkStandard,BD.standard,"Bulk tumour, Standard WGS"),
                        plot.breakpoints(bulkDP,BD.DP,"Bulk tumour, DigiPico")),
             layout_matrix=rbind(c(1),c(2)))
