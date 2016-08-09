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

#####


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

## Segment with negative binomail underlying distribution

install.packages("Segmentor3IsBack")
library(Segmentor3IsBack)
N=2000 
y=c(rnbinom(N,prob=0.3,size=0.15),rnbinom(2*N,prob=0.1,size=0.15),
    rnbinom(N,prob=0.6,size=0.15),compress=FALSE)
res2=Segmentor(bulkDP[chr==6]$total, model=3,Kmax=25) # 3=negative binomial

BPs.bulkDP.segmentor3 <- bulkDP[chr==6][res2@breaks[20,]]$pos

print(plot.breakpoints(bulkDP,BPs.bulkDP.segmentor3,"Segmentor3 - DP",Mode(bulkDP[chr==6]$total)))

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