#' Fit shared component gaussian mixed model using Variational Bayes
#'
#' \code{fit.model.vb} returns parameters for gassians fitted to the data.
#'
#'
#' @param read.count data.frame with two columns, first being numeric read.counts,
#'  and second being charachter segment names. Column names can be whatever.
#'
#' @param rho.input Numeric value for initial value of rho parameter
#' 
#' @param input.parameters Data frame holding initialisation parameters for the prior distribution parameters.
#' With numeric columns "scale" and "shape" the parameters for the gamma
#' prior on precision; "mean" and "nu" for the mean and precision parameters for the normal prior on 
#' the mean; and a charachter column "component.type" with either "common" or "specific"
#' 
#' @param max.iterations maximum number of iterations
#' 
#' @return A list of parameters
#'
#' @examples
#' # see vignette
#' @export
#' @import data.table

fit.model.vb <- function(data, input.parameters, rho.input = 0.5, max.iterations = 40, quiet = FALSE){
  
  ### INITIALISATION
  
  if (any(input.parameters$component.type=="common") & any(input.parameters$component.type=="specific")){
  } else {
    stop("At least one common and one specific component required")
  }
  
  if (all(colnames(input.parameters) != "mean")) stop("Mean values required")
  if (all(colnames(input.parameters) != "nu")) stop("Nu values required")
  if (all(colnames(input.parameters) != "scale")) stop("Scale parameter values required")
  if (all(colnames(input.parameters) != "shape")) stop("Shape parameter values required")
  
  ## split input read count data frame into 2: 
  # 1) a list of indexes specifying which read count is in which segment
  # 2) a vector of read counts
  if(ncol(data) != 2) stop("Two columns of input required")
  setnames(data,names(data),c("vals","seg"))

  if(length(unique(data$seg)) < 2) stop("At least two segments required")

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
kappa <- matrix(rep(500, n.segments),nrow = 1, ncol = n.segments)
N_rho <- matrix(rep(500, n.segments),nrow = 1, ncol = n.segments)

##### reponsibilities
gamma <- gamma_topsum <- matrix(nrow = length(read.count), ncol = n.specific.components)
phi <- matrix(nrow = length(read.count), ncol = n.common.components)

psi <- numeric(length = length(read.count))

##### hyperparameters
# mixing prior
priors <- vector("list", 5)
names(priors) <- c("lambda","mean","nu","scale","shape")
priors$lambda <- 500

# means prior
priors$mean <- mean(input.parameters$mean) #mean(read.count)
priors$nu <- mean(input.parameters$nu) #(sum(abs(range(read.count)))/3)^2

# precisions prior
priors$scale <- mean(input.parameters$scale)
priors$shape <- mean(input.parameters$shape)

##### initialise parameters
# each specific parameter is a j by k matrix
specific_parameters <- vector("list", 6)
names(specific_parameters) <- c("mix_weights","lambda","scale","shape","nu","mean")
specific_parameters$mix_weights <- matrix(rep(1/n.specific.components, n.segments * n.specific.components), ncol = n.specific.components)
specific_parameters$lambda <- matrix(rep(priors$lambda, n.segments * n.specific.components), ncol = n.specific.components)
specific_parameters$scale <- matrix(rep(input.parameters[component.type=="specific"]$scale, n.segments), ncol = n.specific.components)
specific_parameters$shape <- matrix(rep(input.parameters[component.type=="specific"]$shape, n.segments), ncol = n.specific.components)
specific_parameters$nu <- matrix(rep(input.parameters[component.type=="specific"]$nu, n.segments), ncol = n.specific.components)
specific_parameters$mean <- matrix(rep(input.parameters[component.type=="specific"]$mean, each=n.segments), ncol = n.specific.components)

common_parameters <- matrix(nrow = n.common.components, ncol = 6) # vector("list", 6)
colnames(common_parameters) <- c("mix_weights","lambda","scale","shape","nu","mean")
common_parameters <- as.data.table(common_parameters)
common_parameters$mix_weights <- rep((1 / n.common.components), n.common.components )
common_parameters$lambda <- rep(priors$lambda, n.common.components )
common_parameters$scale <- com.param$scale
common_parameters$shape <- com.param$shape
common_parameters$nu <- com.param$nu
common_parameters$mean <- com.param$mean


iter.count <- 0

##### UPDATES
for (round in 1:max.iterations){

#########################
##### E step ############
#########################

# calculate responsibilities for common components

weighting_c = exp(digamma(common_parameters$lambda) - digamma(sum(common_parameters$lambda)))
beta_tilde_c = exp(digamma(common_parameters$shape) + log(common_parameters$scale))
common_topsum <- weighting_c * beta_tilde_c^0.5 *
  exp(-0.5 * common_parameters$scale * common_parameters$shape * ( outer(read.count^2, common_parameters$mean^2 + common_parameters$nu, FUN="+") - 2 * outer(read.count, common_parameters$mean, FUN="*") ) )

phi = common_topsum / rowSums(common_topsum)

# calculate responsibilities for specific components

weighting = exp(digamma(specific_parameters$lambda) - digamma(sum(specific_parameters$lambda)))
beta_tilde = exp(digamma(specific_parameters$shape) + log(specific_parameters$scale))
beta_bar = specific_parameters$scale * specific_parameters$shape

##### what should rho be?
rho_weighting = exp(digamma(kappa) - digamma(sum(kappa)))

for (segment in segment.names){
  indexes <- segment.indicies[[segment]]
  
gamma_topsum[indexes, ] <- weighting[segment, ] * beta_tilde[segment, ]^0.5 *
                            exp(-0.5 * beta_bar[segment, ] * 
                                  (outer(read.count[indexes]^2, specific_parameters$mean[segment, ]^2 + specific_parameters$nu[segment, ], FUN="+") 
                                    - 2 * outer(read.count[indexes], specific_parameters$mean[segment, ], FUN="*")
                                  )
                                )

gamma[indexes, ] <- gamma_topsum[indexes, ]/rowSums(gamma_topsum[indexes, , drop=FALSE])

# calculate probaility for each data point in common or specific component
psi_topsum <- (1-rho_weighting[segment]) * rowSums(common_topsum[indexes, , drop=FALSE])
psi_bottomsum <- psi_topsum + (rho_weighting[segment] * rowSums(gamma_topsum[indexes, , drop=FALSE]))
psi[indexes] <- psi_topsum / psi_bottomsum

#hist(phi[indexes],breaks=100)
#hist(gamma[indexes],breaks=100)
#hist(psi[indexes],breaks=100)
}

#hist(phi,breaks=100)
#hist(gamma,breaks=100)
#hist(psi,breaks=100)

stopifnot(sum(is.na(phi))==0)
stopifnot(sum(is.na(gamma))==0)
stopifnot(sum(is.na(psi))==0)

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
  
  sigma2_specific[segment, ] <- y2_weighted_specific[segment, ] + specific_parameters$mix_weights[segment,] * (specific_parameters$mean[segment, ]^2 + specific_parameters$nu[segment, ]) - 2 * specific_parameters$mean[segment, ] * y_weighted_specific[segment, ]

  specific_parameters$scale[segment, ] <- 1/( sigma2_specific[segment, ] * 0.5 * sum.of.1.minus.psi + 1/priors$scale)
  
  specific_parameters$shape[segment, ] <- N_specific[segment, ] * 0.5 + priors$shape
  
  t_data_specific[segment, ] <- N_specific[segment, ] * specific_parameters$scale[segment, ] * specific_parameters$shape[segment, ]
  
  m_data_specific[segment, ] <- y_weighted_specific[segment, ] / specific_parameters$mix_weights[segment,]
  
  specific_parameters$nu[segment, ] <- 1/ (1/priors$nu + t_data_specific[segment, ])
  
  specific_parameters$mean[segment, ] <- (1/priors$nu * priors$mean) / (1/specific_parameters$nu[segment, ]) + (t_data_specific[segment, ] * m_data_specific[segment, ]) / (1/specific_parameters$nu[segment, ])
}

# calculate parameters common to all segments

sum.of.psi <- sum(psi)

common_parameters$mix_weights <- colSums(phi * psi) / sum.of.psi
N_common <- common_parameters$mix_weights * sum(psi)

y_weighted_common <- colSums( psi * phi * read.count ) / sum.of.psi
y2_weighted_common <- colSums( psi * phi * read.count^2 ) / sum.of.psi

common_parameters$lambda <- N_common + priors$lambda

sigma2_common <- y2_weighted_common + common_parameters$mix_weights * (common_parameters$mean^2 + common_parameters$nu) - 2 * common_parameters$mean * y_weighted_common

common_parameters$scale <- 1/ (sigma2_common * 0.5 * sum.of.psi + 1/priors$scale)

common_parameters$shape <- N_common * 0.5 + priors$shape

t_data <- N_common * (common_parameters$scale * common_parameters$shape)

m_data <- y_weighted_common / common_parameters$mix_weights

common_parameters$nu <- 1/(1/priors$nu + t_data)

common_parameters$mean <- (1/priors$nu * priors$mean) / (1/common_parameters$nu) + (t_data * m_data) / (1/common_parameters$nu)

iter.count <- iter.count + 1
if (quiet == FALSE) print(paste("Iteration:",iter.count))

} # EM updates repetition loop

# parse output

common_parameters$variance = 1/(common_parameters$scale * common_parameters$shape)
specific_parameters$variance = 1/(specific_parameters$scale * specific_parameters$shape)

output <- list(
  common_parameters = common_parameters,
  specific_parameters = specific_parameters,
  rho = rho,
  kappa = kappa,
  iteration = iter.count,
  segment.names = as.character(1:n.segments))

return(output)

}
#########################
##### Done ##############
#########################
