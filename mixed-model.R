library("ggplot2")
#library("dplyr")
#library("reshape2")
library("data.table")

# parse learned/output parameters
parse.output.parameters <- function(output){

  output.parameters <- data.table(rho=numeric(),w=numeric(),rho_w=numeric(),mu=numeric(),sigma2=numeric(),component.type=character(),segment=character(),iteration=numeric())
  
  for (segment in 1:nrow(output$w.specific)){
    
    output.parameters.temp <- data.frame(
      rho = as.numeric(output$rho[segment]),
      w = c(output$w.common,output$w.specific[segment,]),
      rho_w = c((1-output$rho[segment])*output$w.common,output$rho[segment]*output$w.specific[segment,]),
      mu = c(output$mu.common,output$mu.specific[segment,]),
      sigma2 = c(output$sigma2.common,output$sigma2.specific[segment,]),
      component.type = c(rep("common",length(output$mu.common)),rep("specific",length(output$mu.specific[segment,]))),
      segment = segment,
      iteration = output$iteration
    )
    
    output.parameters <- rbind(output.parameters,output.parameters.temp)
  }
  return(output.parameters)
}

# Generic Plot function
plot.components <- function(vals.df,output.parameters){
  
  temp <- data.table(vals=numeric(),comp=character(),seg=character(),component.type=character(),source=character())
  
  for (component in 1:nrow(output.parameters)){
    temp <- rbind(temp,data.table(vals=c(rnorm(70000*output.parameters[component]$rho_w, 
                                             mean = output.parameters[component]$mu, 
                                             sd = output.parameters[component]$sigma2^0.5),0),
                                            comp=component,
                                            seg=output.parameters[component]$segment,
                                            component.type=output.parameters[component]$component.type,
                                            source="Inferred"))
  }
  
  vals.df$component.type <- "Original Data"
  vals.df$source <- "Original Data"
  
  temp <- rbind(temp,vals.df)
  
  # ggplot(temp, aes(vals, group = comp,colour=component.type)) +
  #   geom_freqpoly(binwidth=0.1) +
  #   ggtitle(paste("iteration:",iter.count)) +
  #   scale_colour_manual(values = c("blue","red","black")) +
  #   facet_wrap(~seg,nrow = 2)

  ggplot(temp, aes(vals, colour=component.type)) +
    geom_freqpoly(binwidth=0.1) +
    ggtitle(paste("iteration:",unique(output.parameters$iteration))) +
    scale_colour_manual(values = c("blue","red","black")) +
    facet_grid(source~seg,scales="free_y")

  
  # ggplot(temp, aes(vals, colour=component.type)) +
  #   geom_density() +
  #   ggtitle(paste("iteration:",iter.count)) +
  #   scale_colour_manual(values = c("blue","red","black")) +
  #   facet_wrap(~seg,nrow = 2)
}

# input required: 
# 1) data.frame of read counts and segment name for each position
# 2) initial rho
# 3) data.frame of initial parameters for each component, both common and specific
# fit.model(data=vals.df,rho=0.5,input.parameters=initial.parameters,init.max=40)
# returns a data frame of output parameters

fit.model <- function(vals.df,rho.input,input.parameters,init.max=40){

### INITIALISATION

## split input read count data frame into 2: 
# 1) a list of indexes specifying which read count is in which segment
# 2) a vector of read counts

# get index ranges of each segment
segment.indicies <- vals.df[,.(index = list(.I)),by=seg]
segment.names <- segment.indicies$seg
segment.indicies <- segment.indicies$index
names(segment.indicies) <- segment.names
print(paste(length(segment.indicies),"segments"))

read.count <- vals.df$vals

## parse input parameters

# rho parameter for each segment
rho <- rep(rho.input,length(segment.indicies))
stopifnot(length(segment.indicies)==length(rho))
names(rho) <- segment.names

# set common component parameters, one per common component
com.param <- input.parameters[component.type=="common"]
com.param$component.type <- NULL

n.specific.components <- nrow(input.parameters[component.type=="specific"])
print(paste(n.specific.components,"segment specific components"))
print(paste(nrow(input.parameters[component.type=="specific"]),"common components"))


# initialise segment specific parameter holding matricies
mu.k <- sigma2.k <- w.k <- matrix(nrow=length(read.count),ncol=n.specific.components)

for (row in 1:n.specific.components){
  w.k[,row] <- input.parameters[component.type=="specific"][row]$w
  mu.k[,row] <- input.parameters[component.type=="specific"][row]$mu
  sigma2.k[,row] <- input.parameters[component.type=="specific"][row]$sigma2
}

psi_i <- rep(0.5,length(read.count))

iter.count <- 0

#library(profvis)
#library(microbenchmark)
# profvis({ #CODE HERE })
#microbenchmark(CODE HERE)

##### UPDATES
for (round in 1:init.max){
#### E Updates

# for each common component, cacluate prob for each data point, outputs i x c matrix
common <- apply(com.param, 1, function(params){ params['w'] * dnorm(read.count, params['mu'], params['sigma2']^0.5)})

# output i x k matrix
specific <- w.k * dnorm(read.count, mu.k, sigma2.k^0.5)

for (segment in segment.names){
  indexes <- segment.indicies[[segment]]
  psi_i[indexes] <- ( (1-rho[segment]) * rowSums(common[indexes,]) ) / (( (1-rho[segment]) * rowSums(common[indexes,]) ) + rho[segment] * rowSums(specific[indexes,]) )
}

# sometimes dorm gives 0 when no common/specific components nearby, which can result in NaN when /0, so add a tiny probability back in
common[which(rowSums(common)==0),] <- 1e-15
specific[which(rowSums(specific)==0),] <- 1e-15

phi_ic <- common/rowSums(common)

nu_ik <- specific/rowSums(specific)

# check for na values
stopifnot(sum(is.na(phi_ic))+sum(is.na(nu_ik))==0)

#### M Updates

# returns c parameters
topsum <- colSums(phi_ic * psi_i)

com.param$w <- topsum / sum(psi_i)

com.param$mu <- colSums(phi_ic * psi_i * read.count) / topsum

com.param$sigma2 <- colSums(phi_ic * psi_i * outer(read.count,com.param$mu, FUN="-")^2) / topsum

for (segment in segment.names){
  # for each segment j, there are k weightings
  indexes <- segment.indicies[[segment]]
  
  rho[segment] <- 1 - sum(psi_i[indexes])/length(indexes)
  
  temp.sum <- (1-psi_i[indexes]) * nu_ik[indexes,]
  topsum <- colSums(temp.sum) # for each component...
  
  w.k[indexes,] <- rep(topsum / sum(1-psi_i[indexes]), each=length(indexes))

  mu.k[indexes,] <- rep(colSums(temp.sum * read.count[indexes]) /  topsum, each=length(indexes))
  
  sigma2.k[indexes,] <- 0.01 + rep(colSums(temp.sum * (read.count[indexes] - mu.k[indexes,])^2) /  topsum, each=length(indexes))
  # the + 0.01 prevents sigma2 -> 0 which can cause dnorm() to -> Inf e.g dnorn(47,47,0), which causes NaN (inf/inf) which propogates to everything -> NaN!
}

iter.count <- iter.count + 1
print(iter.count)

} # EM updates repetition loop

w.k.u <- unique(w.k)
mu.k.u <- unique(mu.k)
sigma2.k.u <- unique(sigma2.k)

# more raw output
output <- list(w.specific=w.k.u,mu.specific=mu.k.u,sigma2.specific=sigma2.k.u,w.common=com.param$w,mu.common=com.param$mu,sigma2.common=com.param$sigma2,rho=rho,iteration=iter.count)

return(output)

}# end of fit.model function


###### Example code / Tests

# Simulate baisc well seperated gaussian test data, 2 segments, each with 4 components, 2 shared and 2 unique

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

test.data.basic <- bind_rows(comp1.vals,comp2.vals,comp3.vals,comp4.vals,comp5.vals,comp6.vals,comp7.vals,comp8.vals)

# Overall histogram per segment
ggplot(test.data.basic, aes(vals)) +
  geom_density() +
  facet_wrap(~seg,nrow = 2)

# Histogram broken down by true components
ggplot(test.data.basic, aes(vals, colour = comp)) +
  geom_freqpoly(binwidth=0.1) +
  facet_wrap(~seg,nrow = 2)

initial.parameters.basic <- data.table(
  w=c(0.5),
  mu=c(3.5,10,2,12),
  sigma2=c(0.6^2),
  component.type=c("common","common","specific","specific"))

output.test.basic <- fit.model(test.data.basic,rho=0.5,initial.parameters.basic)

plot.components(test.data.basic,output.test.basic)

### Simulate realistic read count data using negative binomial, 4 segments, 2 normal, 1 amp, 1 del. Parameters from real data

test.data.realistic <- bind_rows(
  data.table(comp = "A", vals = rnbinom(5000, size=2.7, mu = 16.5),seg="normal.1"),
  data.table(comp = "A", vals = rnbinom(5000, size=3, mu = 36),seg="amplified.1"),
  data.table(comp = "A", vals = rnbinom(5000, size=2.9, mu = 8.4),seg="deleted.1"),
  data.table(comp = "A", vals = rnbinom(5000, size=2.7, mu = 16.5),seg="normal.2"))

initial.parameters.realistic <- data.table(
  w=c(0.5),
  mu=c(4,7,10,15,25,40,50,300,1,2),
  sigma2=c(8^2),
  component.type=c(rep("common",5),rep("specific",5)))

output.realistic <- fit.model(test.data.realistic,rho=0.5,initial.parameters.realistic)

plot.components(test.data.realistic,output.realistic)
