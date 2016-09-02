
input.parameters <- data.table(
  lambda=c(5),
  m=c(5,2),
  nu=c(2),
  b=c(4),
  c=c(0.125),
  component.type=c(rep("common",1),rep("specific",1)))

#vals.df <- bulkDP[chr==6 & !is.na(segment.no),.(vals=total,seg=segment.no)]
vals.df <- data.table(vals = c(rnorm(2500, mean = 6.1, sd = 0.75),
                               rnorm(2000, mean = 1, sd = 0.3),
                               rnorm(2500, mean = 5.9, sd = 0.75),
                               rnorm(3000, mean = 2, sd = 0.3)), seg = c(rep(1,4500),rep(2,5500)))

hist(vals.df[seg==1]$vals,breaks=300,xlim=c(-1,9))
hist(vals.df[seg==2]$vals,breaks=300,xlim=c(-1,9))
hist(vals.df$vals,breaks=300,xlim=c(-1,9))


#' Fit shared component gaussian mixed model using Variational Bayes
#'
#' \code{fit.model.vb} returns parameters for gassians fitted to the data.
#'
#'
#' @param read.count data.frame with two columns, first being numeric read.counts,
#'  and second being charachter segment names. Column names can be whatever.
#'
#' @param rho.input Numeric value for initial value of rho parameter
#' @param input.parameters
#' @param init.max
#' @param break.parameter Numeric value, when the difference in log likelihood 
#' between two iterations is less than this value the iteration will stop
#' @return A list of parameters
#'
#' @examples
#' # see vignette
#' @export
#' @import data.table

fit.model.vb <- function(data,rho.input,input.parameters,max.iterations=40){
  
  ### INITIALISATION
  
  if (any(input.parameters$component.type=="common") & any(input.parameters$component.type=="specific")){
  } else {
    stop("At least one common and one specific component required")
  }
  
  ## split input read count data frame into 2: 
  # 1) a list of indexes specifying which read count is in which segment
  # 2) a vector of read counts
  stopifnot(ncol(data)==2)
  setnames(data,names(data),c("vals","seg"))

####################################
##### Parse Input Initialise  ######
####################################

read.count <- data$vals

com.param <- input.parameters[component.type=="common"]
com.param$component.type <- NULL

n.specific.components <- nrow(input.parameters[component.type=="specific"])
print(paste(n.specific.components,"segment specific components"))
n.common.components <- nrow(input.parameters[component.type=="common"])
print(paste(n.common.components,"common components"))

# get index ranges of each segment
segment.indicies <- data[,.(index = list(.I)),by=seg]
segment.names <- segment.indicies$seg
segment.indicies <- segment.indicies$index
names(segment.indicies) <- segment.names
n.segments <- length(segment.indicies)
print(paste(length(segment.indicies),"segments"))

# rho parameter for each segment
rho <- rep(rho.input,length(segment.indicies))
names(rho) <- segment.names

t_data_specific <- m_data_specific <- sigma2_specific <- y2_weighted_specific <- y_weighted_specific <- N_specific <- matrix(nrow = n.segments, ncol = n.specific.components)
t_data <- m_data <- sigma2_common <- y2_weighted_common <- y_weighted_common <- N_common <- matrix(nrow = 1, ncol = n.common.components)
kappa <- matrix(rep(1/n.segments, n.segments),nrow = 1, ncol = n.segments)
N_rho <- matrix(rep(1/n.segments, n.segments),nrow = 1, ncol = n.segments)

##### reponsibilities
gamma <- gamma_topsum <- matrix(nrow = length(read.count), ncol = n.specific.components)
phi <- matrix(nrow = length(read.count), ncol = n.common.components)

psi <- numeric(length = length(read.count))

##### hyperparameters
# mixing prior
priors <- vector("list", 5)
names(priors) <- c("lambda","m","nu","b","c")
priors$lambda <- 5

# means prior
priors$m <- mean(read.count)
priors$nu <- (sum(abs(range(read.count)))/3)^2

# precisions prior
priors$b <- 1000
priors$c <- 0.001

##### initialise parameters
# each specific parameter is a j by k matrix
specific_parameters <- vector("list", 6)
names(specific_parameters) <- c("mix_weights","lambda","b","c","nu","m")
specific_parameters$mix_weights <- matrix(rep(1/n.specific.components, n.segments * n.specific.components), ncol = n.specific.components)
specific_parameters$lambda <- matrix(rep(priors$lambda, n.segments * n.specific.components), ncol = n.specific.components)
specific_parameters$b <- matrix(rep(input.parameters[component.type=="specific"]$b, n.segments), ncol = n.specific.components)
specific_parameters$c <- matrix(rep(input.parameters[component.type=="specific"]$c, n.segments), ncol = n.specific.components)
specific_parameters$nu <- matrix(rep(input.parameters[component.type=="specific"]$nu, n.segments), ncol = n.specific.components)
specific_parameters$m <- matrix(rep(input.parameters[component.type=="specific"]$m, each=n.segments), ncol = n.specific.components)

common_parameters <- matrix(nrow = n.common.components, ncol = 6) # vector("list", 6)
colnames(common_parameters) <- c("mix_weights","lambda","b","c","nu","m")
common_parameters <- as.data.table(common_parameters)
common_parameters$mix_weights <- rep((1 / n.common.components), n.common.components )
common_parameters$lambda <- rep(priors$lambda, n.common.components )
common_parameters$b <- com.param$b
common_parameters$c <- com.param$c
common_parameters$nu <- com.param$nu
common_parameters$m <- com.param$m


iter.count <- 0

##### UPDATES
for (round in 1:max.iterations){

#########################
##### E step ############
#########################

# calculate responsibilities for common components

weighting_c = exp(digamma(common_parameters$lambda) - digamma(sum(common_parameters$lambda)))
beta_tilde_c = exp(digamma(common_parameters$c) + log(common_parameters$b))
common_topsum <- weighting_c * beta_tilde_c^0.5 *
  exp(-0.5 * common_parameters$b * common_parameters$c * ( outer(read.count^2, common_parameters$m^2 + common_parameters$nu, FUN="+") - 2 * outer(read.count, common_parameters$m, FUN="*") ) )

phi = common_topsum / rowSums(common_topsum)

# calculate responsibilities for specific components

weighting = exp(digamma(specific_parameters$lambda) - digamma(sum(specific_parameters$lambda)))
beta_tilde = exp(digamma(specific_parameters$c) + log(specific_parameters$b))
beta_bar = specific_parameters$b * specific_parameters$c

##### what should rho be?
rho_weighting = exp(digamma(kappa) - digamma(sum(kappa)))

for (segment in segment.names){
  indexes <- segment.indicies[[segment]]
  
gamma_topsum[indexes, ] <- weighting[segment, ] * beta_tilde[segment, ]^0.5 *
                            exp(-0.5 * beta_bar[segment, ] * 
                                  (outer(read.count[indexes]^2, specific_parameters$m[segment, ]^2 + specific_parameters$nu[segment, ], FUN="+") 
                                    - 2 * outer(read.count[indexes], specific_parameters$m[segment, ], FUN="*")
                                  )
                                )

gamma[indexes, ] <- gamma_topsum[indexes, ]/rowSums(gamma_topsum[indexes, , drop=FALSE])

# calculate probaility for each data point in common or specific component
psi_topsum <- (1-rho_weighting[segment]) * rowSums(common_topsum[indexes, , drop=FALSE])
psi_bottomsum <- psi_topsum + (rho_weighting[segment] * rowSums(gamma_topsum[indexes, , drop=FALSE]))
psi[indexes] <- psi_topsum / psi_bottomsum

hist(phi[indexes],breaks=100)
hist(gamma[indexes],breaks=100)
hist(psi[indexes],breaks=100)
}

hist(phi,breaks=100)
hist(gamma,breaks=100)
hist(psi,breaks=100)

sum(is.na(phi))
sum(is.na(gamma))
sum(is.na(psi))

#########################
#### M Step #############
#########################

# calculate segment specific parameters in a loop
for (segment in segment.names){
  indexes <- segment.indicies[[segment]]
  
  sum.of.1.minus.psi <- sum(1 - psi[indexes, drop=F])
  
  specific_parameters$mix_weights[segment, ] <- colSums(gamma[indexes, , drop=F] * (1 - psi[indexes])) / sum.of.1.minus.psi
  stopifnot(sum(specific_parameters$mix_weights==0)==0) # or mdata is NaN
  rho[segment] <- 1 - sum(psi[indexes]) / length(indexes)
  
  N_specific[segment, ] <- specific_parameters$mix_weights[segment, ] * sum.of.1.minus.psi
  N_rho[segment] <- rho[segment] * length(indexes)
  
  y_weighted_specific[segment, ] <- colSums( (1 - psi[indexes]) * gamma[indexes, , drop=F] * read.count[indexes] ) / sum.of.1.minus.psi
  y2_weighted_specific[segment, ] <- colSums( (1 - psi[indexes]) * gamma[indexes, , drop=F] * read.count[indexes]^2 ) / sum.of.1.minus.psi

  specific_parameters$lambda[segment, ] <- N_specific[segment, ] + priors$lambda
  kappa[segment] <- N_rho[segment] + priors$lambda
  
  sigma2_specific[segment, ] <- y2_weighted_specific[segment, ] + specific_parameters$mix_weights[segment,] * (specific_parameters$m[segment, ]^2 + specific_parameters$nu[segment, ]) - 2 * specific_parameters$m[segment, ] * y_weighted_specific[segment, ]

  specific_parameters$b[segment, ] <- 1/( sigma2_specific[segment, ] * 0.5 * sum.of.1.minus.psi + 1/priors$b)
  
  specific_parameters$c[segment, ] <- N_specific[segment, ] * 0.5 + priors$c
  
  t_data_specific[segment, ] <- N_specific[segment, ] * specific_parameters$b[segment, ] * specific_parameters$c[segment, ]
  
  m_data_specific[segment, ] <- y_weighted_specific[segment, ] / specific_parameters$mix_weights[segment,]
  
  specific_parameters$nu[segment, ] <- 1/ (1/priors$nu + t_data_specific[segment, ])
  
  specific_parameters$m[segment, ] <- (1/priors$nu * priors$m) / (1/specific_parameters$nu[segment, ]) + (t_data_specific[segment, ] * m_data_specific[segment, ]) / (1/specific_parameters$nu[segment, ])
}

# calculate parameters common to all segments

sum.of.psi <- sum(psi)

common_parameters$mix_weights <- colSums(phi * psi) / sum.of.psi
N_common <- common_parameters$mix_weights * sum(psi)
print(paste("sum of N Common:", sum(N_common)))

y_weighted_common <- colSums( psi * phi * read.count ) / sum.of.psi
y2_weighted_common <- colSums( psi * phi * read.count^2 ) / sum.of.psi

common_parameters$lambda <- N_common + priors$lambda

sigma2_common <- y2_weighted_common + common_parameters$mix_weights * (common_parameters$m^2 + common_parameters$nu) - 2 * common_parameters$m * y_weighted_common

common_parameters$b <- 1/ (sigma2_common * 0.5 * sum.of.psi + 1/priors$b)

common_parameters$c <- N_common * 0.5 + priors$c

t_data <- N_common * (common_parameters$b * common_parameters$c)

m_data <- y_weighted_common / common_parameters$mix_weights

common_parameters$nu <- 1/(1/priors$nu + t_data)

common_parameters$m <- (1/priors$nu * priors$m) / (1/common_parameters$nu) + (t_data * m_data) / (1/common_parameters$nu)

iter.count <- iter.count + 1
print(paste("Iteration:",iter.count))

} # EM updates repetition loop

# parse output

output <- list(common_parameters=common_parameters,specific_parameters=specific_parameters,rho=rho,kappa=kappa)

return(output)

}
#########################
##### Done ##############
#########################

out <- fit.model.vb(vals.df,0.5,input.parameters)

y_infered <- rbind(data.table(y=rnorm(ceiling(N_common[1])/n.segments, mean = common_parameters$m[1], sd = 1/(common_parameters$b * common_parameters$c)[1]),source="C1"),
                   data.table(y=rnorm(ceiling(N_specific[1,1]), mean = specific_parameters$m[1,1], sd = 1/(specific_parameters$b * specific_parameters$c)[1,1]),source="S1"))

#data.table(y=rnorm(ceiling(N_specific[2,2]), mean = specific_parameters$m[2,2], sd = 1/(specific_parameters$b * specific_parameters$c)[1,2]),source="S2")
#data.table(y=rnorm(ceiling(N_common[2]), mean = common_parameters$m[2], sd = 1/(common_parameters$b * common_parameters$c)[2]),source="C2"),

toplot <- rbind(data.frame(y=vals.df[seg==1]$vals,source="original"),y_infered)
ggplot(toplot,aes(y,colour=source)) + geom_freqpoly(bins=200)

  
