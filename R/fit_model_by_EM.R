

# initialise segment specific parameter holding matricies
mu.k <- sigma2.k <- w.k <- matrix(nrow=length(read.count),ncol=n.specific.components)
# although these could be a simple j*k matrix, we're already iterating through segments in
# the M step so may as well create the whole matrix required by the dnorm E step below 

for (row in 1:n.specific.components){
  w.k[,row] <- input.parameters[component.type=="specific"][row]$w
  mu.k[,row] <- input.parameters[component.type=="specific"][row]$mean
  sigma2.k[,row] <- input.parameters[component.type=="specific"][row]$sigma2
}

psi_i <- rep(0.5,length(read.count))
likelihood.byiter <- numeric()

iter.count <- 0

#library(profvis)
#library(microbenchmark)
# profvis({ #CODE HERE })
#microbenchmark(CODE HERE)

##### UPDATES
for (round in 1:max.iterations){
#### E Updates

# for each common component, cacluate prob for each data point, outputs i x c matrix
common <- apply(com.param, 1, function(params){ params['w'] * dnorm(read.count, params['mean'], params['sigma2']^0.5)})

# output i x k matrix
specific <- w.k * dnorm(read.count, mu.k, sigma2.k^0.5)

# sometimes dorm gives 0 when no common/specific components nearby, which can result in NaN when /0, so add a tiny probability back in
common[which(rowSums(common)==0),] <- 1e-15
specific[which(rowSums(specific)==0),] <- 1e-15

likelihood.vec <- numeric()

for (segment in segment.names){
  indexes <- segment.indicies[[segment]]
  top.sum <- (1-rho[segment]) * rowSums(common[indexes,,drop=FALSE])
  bottom.sum <- top.sum + (rho[segment] * rowSums(specific[indexes,,drop=FALSE]))
  psi_i[indexes] <- top.sum / ( bottom.sum )
  likelihood.vec <- c(likelihood.vec,bottom.sum)
}

likelihood.byiter <- c(likelihood.byiter,sum(log(likelihood.vec)))


phi_ic <- common/rowSums(common)

nu_ik <- specific/rowSums(specific)

# check for na values
stopifnot(sum(is.na(phi_ic))+sum(is.na(nu_ik))+sum(is.na(phi_ic))==0)

#### M Updates

# returns c parameters
topsum <- colSums(phi_ic * psi_i)

com.param$w <- topsum / sum(psi_i)

com.param$mean <- colSums(phi_ic * psi_i * read.count) / topsum

com.param$sigma2 <- colSums(phi_ic * psi_i * outer(read.count,com.param$mean, FUN="-")^2) / topsum

for (segment in segment.names){
  # for each segment j, there are k weightings
  indexes <- segment.indicies[[segment]]
  
  rho[segment] <- 1 - sum(psi_i[indexes])/length(indexes)
  
  temp.sum <- (1-psi_i[indexes]) * nu_ik[indexes,,drop=FALSE]
  topsum <- colSums(temp.sum) # for each component...
  
  w.k[indexes,] <- rep(topsum / sum(1-psi_i[indexes]), each=length(indexes))

  mu.k[indexes,] <- rep(colSums(temp.sum * read.count[indexes]) /  topsum, each=length(indexes))
  
  sigma2.k[indexes,] <- 0.01 + rep(colSums(temp.sum * (read.count[indexes] - mu.k[indexes,,drop=FALSE])^2) /  topsum, each=length(indexes))
  # the + 0.01 prevents sigma2 -> 0 which can cause dnorm() to -> Inf e.g dnorn(47,47,0), which causes NaN (inf/inf) which propogates to everything -> NaN!
}


iter.count <- iter.count + 1
print(paste("Iteration:",iter.count))

if (iter.count>1){
if (min(abs(diff(likelihood.byiter))) < break.parameter){print("Done, saving parameter estimates")
  break}
}

} # EM updates repetition loop

# parse component specific parameters
w.k.u <- mu.k.u <- sigma2.k.u <- matrix(nrow=length(segment.names),ncol=n.specific.components)

for (segment in 1:length(segment.names)){
  # for each segment j, there are k weightings
  indexes <- segment.indicies[[segment.names[segment]]]
  w.k.u[segment,] <- unique(w.k[indexes,])
  mu.k.u[segment,] <- unique(mu.k[indexes,])
  sigma2.k.u[segment,] <- unique(sigma2.k[indexes,])
}

# more raw output
output <- list(
	specific_parameters = list(
		mix_weights = w.k.u,
		mean = mu.k.u,
		variance = sigma2.k.u),
	common_parameters = data.table(
		mix_weights = com.param$w,
		mean = com.param$mean,
		variance = com.param$sigma2),
	rho = unname(rho),
	iteration = iter.count,
	likelihood = likelihood.byiter,
	segment.names = segment.names)


