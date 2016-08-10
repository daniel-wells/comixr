###### Mixture Modeling

library(data.table)
# install.packages("mclust")
library(mclust)

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# plot histogram of individual bins
plot.histogram <- function(data,lower.bound,upper.bound,title){
  hist(data[chr==6 & pos>lower.bound & pos<upper.bound]$total,breaks=1000,
       xlim=c(0,60),
       main=title,
       xlab="Total Read Count")
  abline(v=mean(bulkDP[chr==6 & pos>lower.bound & pos<upper.bound]$total),col="blue")
  abline(v=median(bulkDP[chr==6 & pos>lower.bound & pos<upper.bound]$total),col="green")
  abline(v=Mode(bulkDP[chr==6 & pos>lower.bound & pos<upper.bound]$total),col="red")
}

layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
plot.histogram(bulkDP,0,18568208,"Normal CN, (bulk digiPico)") # normal
plot.histogram(bulkDP,64622787,133557441,"Normal CN, (bulk digiPico)") # normal
plot.histogram(bulkDP,165193244,1165193244,"Deleted reigon, (bulk digiPico)") # del
plot.histogram(bulkDP,47653485,57239589,"Amplified reigon, (bulk digiPico)") # amp

# reset layout
par(mfrow=c(1,1))
dev.off()


plot.gmm <- function(data,lower.bound,upper.bound,title){
  mixture.model <- Mclust(bulkDP[chr==6 & pos>lower.bound & pos<upper.bound]$total)
  xseq <- seq(0,80,0.1)
  
  m <- matrix(, nrow = length(mixture.model$parameters$mean), ncol = length(xseq))
  for (i in 1:length(mixture.model$parameters$mean)){
    m[i,] <- mixture.model$parameters$pro[i]*dnorm(xseq, mixture.model$parameters$mean[i], mixture.model$parameters$variance$sigmasq[i]^0.5)
  }
  
  matplot(xseq,t(m), type="l",ylim=c(0,max(colSums(m)+0.005)),col="blue",main=paste0(title," Gs:",mixture.model$G),xlab="Total Read Count",ylab="Density")
  lines(xseq,colSums(m))
  lines(as.numeric(table(bulkDP[chr==6 & pos>lower.bound & pos<upper.bound]$total)/sum(table(bulkDP[chr==6 & pos>lower.bound & pos<upper.bound]$total))),col="red")
  legend(40,0.04,c("Components","Overall MM Desnity","Actual Data Density"),lty=c(1,1,1),lwd=c(2.5,2.5,2.5),col=c("blue","black","red"))
}

plot.gmm(bulkDP,0,18568208,"Normal CN, (bulk digiPico)")

plot.gmm(bulkDP,64622787,133557441,"Normal CN, (bulk digiPico)")

plot.gmm(bulkDP,165193244,1165193244,"Deleted reigon, (bulk digiPico)")

plot.gmm(bulkDP,47653485,57239589,"Amplified reigon, (bulk digiPico)")

# fit GMM to data as a whole
normal.model2 <- Mclust(bulkDP$total)
