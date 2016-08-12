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

# main plotting function
plot.breakpoints <- function(data,breakpoints,title="",average){
 P <- ggplot(data, aes(pos,total))+
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

