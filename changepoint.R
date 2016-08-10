# install.packages("data.table")
# install.packages("caTools")
# install.packages("gridExtra")

library(data.table)
library(ggplot2)
library(caTools)
library(gridExtra)

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

summarise.data <- function(data){
  data[,non.ref:=sum(a,c,t,g),by=c("chr","pos")]
  data[,total:=sum(ref,a,c,t,g),by=c("chr","pos")]
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


# plot binned median over points for WGS and DigiPico
grid.arrange(grobs=list(plot.genomic.loadings(bulkStandard,"Bulk tumour, Standard WGS"),
                        plot.genomic.loadings(bulkDP,"Bulk tumour, DigiPico Sequencing")),
             layout_matrix=rbind(c(1),c(2)))


####
#### changepoint
####
# install.packages("changepoint")
library("changepoint")

cptST <- cpt.meanvar(bulkStandard[chr==6]$total,test.stat="Normal",method="BinSeg",Q=20,minseglen=1500)
cptDP <- cpt.meanvar(bulkDP[chr==6]$total,test.stat="Poisson",method="BinSeg",Q=30,minseglen=3000)

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
  model.output$segment.no <- as.integer(model.output$segment)
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

# common mu's from mclust on whole genome (+2 at start)
initial.parameters.realistic <- data.table(
  w=c(0.5),
  mu=c(2,2.61,4.81,7.72,10.44,13.00,17.17,23.84,38.32,50,60,70,300,1),
  sigma2=c(8^2),
  component.type=c(rep("common",9),rep("specific",5)))


### fit actual data to model
chr6.input <- bulkDP[chr==6 & segment.no %in% c(1:31),.(vals=total,seg=as.character(segment.no))]
chr6.output <- fit.model(chr6.input,rho=0.5,initial.parameters.realistic)
plot(chr6.output$likelihood[2:length(chr6.output$likelihood)])
chr6.input$comp <- "A"
plot.components(chr6.input[seg==23],chr6.output,segment.subset = "23")
plot.components(chr6.input[seg==22],chr6.output,segment.subset = "22")
plot.components(chr6.input[seg==15],chr6.output,segment.subset = "15")
plot.components(chr6.input[seg==31],chr6.output,segment.subset = "31")
plot.components(chr6.input[seg==5],chr6.output,segment.subset = "5")

######


# compare segment modes of normal vs digipico sequencing
grid.arrange(grobs=list(plot.breakpoints(bulkStandard,BPs.bulkStandard,"BinSeg - standard",Mode(bulkStandard[chr==6]$total)),
                        plot.breakpoints(bulkDP,BPs.bulkDP,"BinSeg - DP",Mode(bulkDP[chr==6]$total))),
             layout_matrix=rbind(c(1),c(2)))

# compare rho to mode per segment
grid.arrange(grobs=list(plot.rho(bulkDP[chr==6],parse.output.parameters(chr6.output,segment.subset = NULL)),
                        plot.breakpoints(bulkDP,BPs.bulkDP,"BinSeg - DP",Mode(bulkDP[chr==6]$total))),
             layout_matrix=rbind(c(1),c(2)))

#### Check distribution of per segment read count

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

# QQ Plots
#install.packages("car")
library(car)
qqp(bulkDP[chr==6 & pos>lower.bound & pos<upper.bound]$total, "nbinom",size=2.67,mu=16.2)
qqp(bulkDP[chr==6 & pos>lower.bound & pos<upper.bound]$total, "pois",lambda=16.2)


