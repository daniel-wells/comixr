library("ggplot2")
library("dplyr")
library("reshape2")
library("data.table")

options(scipen = 999) # don't print scientific notation numbers
set.seed(1)

# Simulate data, two segments, each with two components, a shared and a unique

comp1.vals <- data.table(comp = "A", vals = rnorm(1000, mean = 1, sd = 0.5),seg="seg1")
comp2.vals <- data.table(comp = "B", vals = rnorm(2000, mean = 2, sd = 0.5),seg="seg1")

comp3.vals <- data.table(comp = "B", vals = rnorm(2000, mean = 2, sd = 0.5),seg="seg2")
comp4.vals <- data.table(comp = "C", vals = rnorm(1000, mean = 3, sd = 0.5),seg="seg2")

comp5.vals <- data.table(comp = "D", vals = rnorm(1500, mean = 5, sd = 0.5),seg="seg1")
comp6.vals <- data.table(comp = "D", vals = rnorm(1500, mean = 5, sd = 0.5),seg="seg2")

vals.df <- bind_rows(comp1.vals, comp2.vals, comp3.vals, comp4.vals,comp5.vals,comp6.vals)

# hist per segment
ggplot(vals.df, aes(vals)) +
  geom_density() +
  facet_wrap(~seg,nrow = 2)

# broken down by true components
ggplot(vals.df, aes(vals, colour = comp)) +
  geom_density() +
  facet_wrap(~seg,nrow = 2)

ggplot(vals.df, aes(vals, colour = comp)) +
  geom_freqpoly(binwidth=0.1) +
  facet_wrap(~seg,nrow = 2)

# broken down by true components
plot.components <- function(){
ggplot(vals.df, aes(vals, colour = comp)) +
  geom_density() +
  facet_wrap(~seg,nrow = 2) +
  ggtitle(paste("iteration:",iter.count)) +
  stat_function(fun = dnorm, colour = "black",aes(linetype="common"), args = list(mean = com.param[1,]$mu, sd = com.param[1,]$sigma^0.5)) +
  stat_function(fun = dnorm, colour = "black",aes(linetype="common"), args = list(mean = com.param[2,]$mu, sd = com.param[2,]$sigma^0.5)) +
  stat_function(fun = dnorm, colour = "black",aes(linetype="segment\n specific"), args = list(mean = unique(vals.df[seg=="seg1"]$k1mu), sd = unique(vals.df[seg=="seg1"]$k1sigma^0.5))) +
  stat_function(fun = dnorm, colour = "black",aes(linetype="segment\n specific"), args = list(mean = unique(vals.df[seg=="seg2"]$k1mu), sd = unique(vals.df[seg=="seg2"]$k1sigma^0.5)))
}

plot.components()

### INITIALISATION

rho <- 0.5

# set common component parameters
com.param <- data.table(w=c(0.5,0.5),
                          mu=c(1.5,3.5),
                          sigma=c(0.6^2,0.6^2))

vals.df[seg=="seg1",k1w:=1] # something wrong here
vals.df[seg=="seg1",k1mu:=0.5] # actually 1
vals.df[seg=="seg1",k1sigma:=0.3^2] # actually 0.5

vals.df[seg=="seg2",k1w:=1]
vals.df[seg=="seg2",k1mu:=4] # actually 3
vals.df[seg=="seg2",k1sigma:=0.6^2] # actually 0.5

iter.count <- 0

plot.components()

#library(profvis)
#library(microbenchmark)
# profvis({ CODE HERE })
#microbenchmark(CODE HERE)

# UPDATES
# for the case of c=1, k=1 (per segment)

# for each common component, cacluate prob for each data point, outputs i x c matrix
common <- apply(com.param, 1, function(params){ params['w'] * dnorm(vals.df$vals, params['mu'], params['sigma']^0.5)})

specific <- vals.df$k1w * dnorm(vals.df$vals, vals.df$k1mu, vals.df$k1sigma^0.5)

vals.df$psi_i <- ( (1-rho) * rowSums(common) ) / (( (1-rho) * rowSums(common) ) + rho*specific )

phi_ic <- common/rowSums(common)

vals.df$v_ik <- 1 # specific/specific

####
  
rho <- 1 - sum(vals.df$psi_i)/nrow(vals.df)

# returns c parameters
topsum <- colSums(phi_ic * vals.df$psi_i)

com.param$w <- topsum / sum(vals.df$psi_i)

com.param$mu <- colSums(phi_ic * vals.df$psi_i * vals.df$vals) / topsum

# mu_0c <- sum(vals.df$phi_ic * vals.df$psi_i * vals.df$vals) / topsum

com.param$sigma <- colSums(phi_ic * vals.df$psi_i * outer(vals.df$vals,com.param$mu, FUN="-")^2) / topsum

# sigma_0c <- sum(vals.df$phi_ic * vals.df$psi_i * (vals.df$vals-mu_0c)^2) / topsum

for (segment in c("seg1","seg2")){
  # for each segment j, there are k weightings
  topsum <- sum((1-vals.df[seg==segment]$psi_i) * vals.df[seg==segment]$v_ik)
  
  vals.df[seg==segment,k1w:= topsum / sum(1-vals.df[seg==segment]$psi_i)]
  
  temp.sum <- 1-vals.df[seg==segment]$psi_i * vals.df[seg==segment]$v_ik
  
  vals.df[seg==segment,k1mu:=sum(temp.sum * vals.df[seg==segment]$vals) /  topsum]
  
  vals.df[seg==segment,k1sigma:=sum(temp.sum * (vals.df[seg==segment]$vals - vals.df[seg==segment]$k1mu)^2) /  topsum]
}

iter.count <- iter.count + 1

plot.components()

