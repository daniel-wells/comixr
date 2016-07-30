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
ggplot(vals.df, aes(vals, colour = comp)) +
  geom_density() +
  facet_wrap(~seg,nrow = 2) +
  ggtitle(paste("iteration:",iter.count)) +
  stat_function(fun = dnorm, colour = "black",aes(linetype="common"), args = list(mean = com.param[1,]$mu, sd = com.param[1,]$sigma^0.5)) +
  stat_function(fun = dnorm, colour = "black",aes(linetype="common"), args = list(mean = com.param[2,]$mu, sd = com.param[2,]$sigma^0.5)) +
  stat_function(fun = dnorm, colour = "black",aes(linetype="segment\n specific"), args = list(mean = unique(mu.k[segment.indicies[['seg1']],1]), sd = unique(sigma2.k[segment.indicies[['seg1']],1])^0.5)) +
    stat_function(fun = dnorm, colour = "black",aes(linetype="segment\n specific"), args = list(mean = unique(mu.k[segment.indicies[['seg1']],2]), sd = unique(sigma2.k[segment.indicies[['seg1']],2])^0.5)) +
    stat_function(fun = dnorm, colour = "black",aes(linetype="segment\n specific"), args = list(mean = unique(mu.k[segment.indicies[['seg2']],1]), sd = unique(sigma2.k[segment.indicies[['seg1']],1])^0.5)) +
    stat_function(fun = dnorm, colour = "black",aes(linetype="segment\n specific"), args = list(mean = unique(mu.k[segment.indicies[['seg2']],2]), sd = unique(sigma2.k[segment.indicies[['seg2']],2])^0.5))
}

plot.components()

### INITIALISATION

# get index ranges of each segment
segment.indicies <- vals.df[,.(index = list(.I)),by=seg]
segment.names <- segment.indicies$seg
segment.indicies <- segment.indicies$index
names(segment.indicies) <- segment.names

# for each segment
rho <- c(0.5,0.5)
names(rho) <- segment.names

# set common component parameters, one per common component
com.param <- data.table(w=c(0.5,0.5),
                        mu=c(3.5,10),
                        sigma2=c(0.6^2,0.6^2))

# initialise segment specific parameter holding matricies
mu.k <- sigma2.k <- w.k <- matrix(nrow=nrow(vals.df),ncol=2)

w.k[,1] <- 0.5
w.k[,2] <- 0.5

mu.k[,1] <- 2 # actually 1
mu.k[,2] <- 12 # actually 3

sigma2.k[,1] <- 0.5^2 # actually 0.5
sigma2.k[,2] <- 0.5^2 # actually 0.5

vals.df$psi_i <- 0.5

iter.count <- 0

plot.components()

#library(profvis)
#library(microbenchmark)
# profvis({ CODE HERE })
#microbenchmark(CODE HERE)

# UPDATES

#### E Updates

# for each common component, cacluate prob for each data point, outputs i x c matrix
common <- apply(com.param, 1, function(params){ params['w'] * dnorm(vals.df$vals, params['mu'], params['sigma2']^0.5)})

# output i x k matrix
specific <- w.k * dnorm(vals.df$vals, mu.k, sigma2.k^0.5)

for (segment in segment.names){
  indexes <- segment.indicies[[segment]]
  vals.df[seg==segment]$psi_i <- ( (1-rho[segment]) * rowSums(common[indexes,]) ) / (( (1-rho[segment]) * rowSums(common[indexes,]) ) + rho[segment] * rowSums(specific[indexes,]) )
}

phi_ic <- common/rowSums(common)

nu_ik <- specific/rowSums(specific)

#### M Updates

# returns c parameters
topsum <- colSums(phi_ic * vals.df$psi_i)

com.param$w <- topsum / sum(vals.df$psi_i)

com.param$mu <- colSums(phi_ic * vals.df$psi_i * vals.df$vals) / topsum

com.param$sigma2 <- colSums(phi_ic * vals.df$psi_i * outer(vals.df$vals,com.param$mu, FUN="-")^2) / topsum

for (segment in segment.names){
  # for each segment j, there are k weightings
  indexes <- segment.indicies[[segment]]
  
  rho[segment] <- 1 - sum(vals.df[seg==segment]$psi_i)/length(indexes)
  
  temp.sum <- (1-vals.df[seg==segment]$psi_i) * nu_ik[indexes,]
  topsum <- colSums(temp.sum) # for each component...
  
  w.k[indexes,] <- rep(topsum / sum(1-vals.df[seg==segment]$psi_i), each=length(indexes))

  mu.k[indexes,] <- rep(colSums(temp.sum * vals.df[seg==segment]$vals) /  topsum, each=length(indexes))
  
  sigma2.k[indexes,] <- rep(colSums(temp.sum * (vals.df[seg==segment]$vals - mu.k[indexes,])^2) /  topsum, each=length(indexes))

}

iter.count <- iter.count + 1

plot.components()

