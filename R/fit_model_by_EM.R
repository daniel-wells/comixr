EM <- function(segment.indicies, read.count, rho, com.param, n.specific.components, n.common.components, n.segments, input.parameters, max.iterations, break.parameter, quiet, segment.names){

w.k <- matrix(rep(input.parameters[component.type=="specific"]$w, each=n.segments), ncol = n.specific.components)
mu.k <- matrix(rep(input.parameters[component.type=="specific"]$mean, each=n.segments), ncol = n.specific.components)
sigma2.k <- matrix(rep(input.parameters[component.type=="specific"]$sigma2, each=n.segments), ncol = n.specific.components)

rownames(sigma2.k) <- rownames(mu.k) <- rownames(w.k) <- segment.names

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

# sometimes dorm gives 0 when no common/specific components nearby, which can result in NaN when /0, so add a tiny probability back in
common[which(rowSums(common)==0),] <- 1e-15

likelihood.vec <- numeric()
specific <- matrix(nrow=length(read.count),ncol=n.specific.components)
# output i x k matrix
for (segment in segment.names){
  indexes <- segment.indicies[[segment]]
  
  specific[indexes,] <- sweep(exp( sweep(-outer(read.count[indexes],mu.k[segment,],FUN="-")^2, 2, (2*sigma2.k[segment,,drop=F]),"/")), 2, (1/(sigma2.k[segment,]^0.5 * sqrt(2*pi))) * w.k[segment,], "*")
  
  # sometimes dorm gives 0 when no common/specific components nearby, which can result in NaN when /0, so add a tiny probability back in
  specific[which(rowSums(specific[indexes,,drop=F])==0),] <- 1e-15
  
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
  
  w.k[segment,] <- topsum / sum(1-psi_i[indexes])
  
  mu.k[segment,] <- colSums(temp.sum * read.count[indexes]) /  topsum
  
  sigma2.k[segment,] <- 0.01 + colSums(temp.sum * outer(read.count[indexes], mu.k[segment,],FUN="-")^2) / topsum
  # the + 0.01 prevents sigma2 -> 0 which can cause dnorm() to -> Inf e.g dnorn(47,47,0), which causes NaN (inf/inf) which propogates to everything -> NaN!
}


iter.count <- iter.count + 1
print(paste("Iteration:",iter.count))

if (iter.count>1){
  if (min(abs(diff(likelihood.byiter))) < break.parameter){print("Done, saving parameter estimates")
    break}
}

} # EM updates repetition loop

# more raw output
output <- list(
specific_parameters = list(
  mix_weights = w.k,
  mean = mu.k,
  variance = sigma2.k),
common_parameters = data.table(
  mix_weights = com.param$w,
  mean = com.param$mean,
  variance = com.param$sigma2),
rho = unname(rho),
iteration = iter.count,
likelihood = likelihood.byiter,
segment.names = segment.names)

return(output)
}