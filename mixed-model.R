library("ggplot2")
library("dplyr")
library("reshape2")
library("data.table")

options(scipen = 999) # don't print scientific notation numbers
set.seed(1)

# Simulate test data, 2 segments, each with 4 components, 2 shared and 2 unique

# Common components
comp1.vals <- data.table(comp = "A", vals = rnorm(2000, mean = 5, sd = 0.5),seg="seg1")
comp2.vals <- data.table(comp = "A", vals = rnorm(2000, mean = 5, sd = 0.5),seg="seg2")

comp3.vals <- data.table(comp = "B", vals = rnorm(1500, mean = 9, sd = 0.5),seg="seg1")
comp4.vals <- data.table(comp = "B", vals = rnorm(1500, mean = 9, sd = 0.5),seg="seg2")

# Unique components
comp5.vals <- data.table(comp = "C", vals = rnorm(3000, mean = 1, sd = 0.5),seg="seg1")
comp6.vals <- data.table(comp = "D", vals = rnorm(1000, mean = 3, sd = 0.5),seg="seg2")

comp7.vals <- data.table(comp = "E", vals = rnorm(1500, mean = 10, sd = 0.5),seg="seg1")
comp8.vals <- data.table(comp = "F", vals = rnorm(1500, mean = 12, sd = 0.5),seg="seg2")

vals.df <- bind_rows(comp1.vals,comp2.vals,comp3.vals,comp4.vals,comp5.vals,comp6.vals,comp7.vals,comp8.vals)
vals.df$source <- "Known"

# Overall histogram per segment
ggplot(vals.df, aes(vals)) +
  geom_density() +
  facet_wrap(~seg,nrow = 2)

# Histogram broken down by true components
ggplot(vals.df, aes(vals, colour = comp)) +
  geom_density() +
  facet_wrap(~seg,nrow = 2)

# including weightings
ggplot(vals.df, aes(vals, colour = comp)) +
  geom_freqpoly(binwidth=0.1) +
  facet_wrap(~seg,nrow = 2)

# Generic Plot function
plot.components <- function(){
  temp <- data.table(vals=numeric(),comp=character(),seg=character())
  
  for (segment in segment.names){
    indexes <- segment.indicies[[segment]]
    
    w.temp <- c((1-rho[segment])*com.param$w,rho[segment]*unique(w.k[indexes,]))
    mu.temp <- c(com.param$mu,unique(mu.k[indexes,]))
    sigma2.temp <- c(com.param$sigma2,unique(sigma2.k[indexes,]))
    
    for (component in 1:length(w.temp)){
      temp <- rbind(temp,data.table(vals=rnorm(7000*w.temp[component], mean = mu.temp[component], sd = sigma2.temp[component]^0.5),comp=component,seg=segment))
    }
    
  }
  
  temp$source <- "Inferred"
  
  temp <- rbind(temp,vals.df)
  
  ggplot(temp, aes(vals, group = comp,colour=source)) +
    geom_freqpoly(binwidth=0.1) +
    ggtitle(paste("iteration:",iter.count)) +
    facet_wrap(~seg,nrow = 2)
}

plot.components()

### INITIALISATION

# get index ranges of each segment
segment.indicies <- vals.df[,.(index = list(.I)),by=seg]
segment.names <- segment.indicies$seg
segment.indicies <- segment.indicies$index
names(segment.indicies) <- segment.names

yi.vals <- vals.df$vals

# for each segment
rho <- c(0.5,0.5)
names(rho) <- segment.names

# set common component parameters, one per common component
com.param <- data.table(w=c(0.5,0.5),
                        mu=c(3.5,10),
                        sigma2=c(0.6^2,0.6^2))

# initialise segment specific parameter holding matricies
mu.k <- sigma2.k <- w.k <- matrix(nrow=length(yi.vals),ncol=2)

w.k[,1] <- 0.5
w.k[,2] <- 0.5

mu.k[,1] <- 2 # actually 1
mu.k[,2] <- 12 # actually 3

sigma2.k[,1] <- 0.5^2 # actually 0.5
sigma2.k[,2] <- 0.5^2 # actually 0.5

psi_i <- rep(0.5,length(yi.vals))

iter.count <- 0

plot.components()

#library(profvis)
#library(microbenchmark)
# profvis({ #CODE HERE })
#microbenchmark(CODE HERE)

# UPDATES
for (round in 1:40){
#### E Updates

# for each common component, cacluate prob for each data point, outputs i x c matrix
common <- apply(com.param, 1, function(params){ params['w'] * dnorm(yi.vals, params['mu'], params['sigma2']^0.5)})

# output i x k matrix
specific <- w.k * dnorm(yi.vals, mu.k, sigma2.k^0.5)

for (segment in segment.names){
  indexes <- segment.indicies[[segment]]
  psi_i[indexes] <- ( (1-rho[segment]) * rowSums(common[indexes,]) ) / (( (1-rho[segment]) * rowSums(common[indexes,]) ) + rho[segment] * rowSums(specific[indexes,]) )
}

phi_ic <- common/rowSums(common)

nu_ik <- specific/rowSums(specific)

#### M Updates

# returns c parameters
topsum <- colSums(phi_ic * psi_i)

com.param$w <- topsum / sum(psi_i)

com.param$mu <- colSums(phi_ic * psi_i * yi.vals) / topsum

com.param$sigma2 <- colSums(phi_ic * psi_i * outer(yi.vals,com.param$mu, FUN="-")^2) / topsum

for (segment in segment.names){
  # for each segment j, there are k weightings
  indexes <- segment.indicies[[segment]]
  
  rho[segment] <- 1 - sum(psi_i[indexes])/length(indexes)
  
  temp.sum <- (1-psi_i[indexes]) * nu_ik[indexes,]
  topsum <- colSums(temp.sum) # for each component...
  
  w.k[indexes,] <- rep(topsum / sum(1-psi_i[indexes]), each=length(indexes))

  mu.k[indexes,] <- rep(colSums(temp.sum * yi.vals[indexes]) /  topsum, each=length(indexes))
  
  sigma2.k[indexes,] <- rep(colSums(temp.sum * (yi.vals[indexes] - mu.k[indexes,])^2) /  topsum, each=length(indexes))

}

iter.count <- iter.count + 1

plot.components()

} # EM repetition loop
