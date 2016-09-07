#' Fit shared component gaussian mixed model.
#'
#' \code{fit.model} returns parameters for gassians fitted to the data.
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
#' @param break.parameter Numeric value, when the difference in log likelihood 
#' between two iterations is less than this value the iteration will stop
#' 
#' @param algorithm Charachter specifying which inference algorithm to use.
#' One of either "EM" (Expectation-Maximisation) or "VB" (Variational-Bayes)
#' 
#' @param quiet, Logical specifying whether the iteration number should be printed for each iteration
#' 
#' @return A list of parameters
#'
#' @examples
#' # see vignette
#' @export
#' @import data.table

fit.model <- function(data, input.parameters, rho = 0.5, max.iterations = 40, break.parameter = 5, algorithm = "EM", quiet = FALSE){

  ### INITIALISATION
  
  if (any(input.parameters$component.type=="common") & any(input.parameters$component.type=="specific")){
  } else {
    stop("At least one common and one specific component required")
  }
  
  if (all(colnames(input.parameters) != "mean")) stop("Mean values required")
  
  ## split input read count data frame into 2: 
  # 1) a list of indexes specifying which read count is in which segment
  # 2) a vector of read counts
  if(ncol(data) != 2) stop("Two columns of input required")
  setnames(data,names(data),c("vals","seg"))
  
  if(length(unique(data$seg)) < 2) stop("At least two segments required")
  
  # get index ranges of each segment
  segment.indicies <- data[,.(index = list(.I)),by=seg]
  segment.names <- segment.indicies$seg
  segment.indicies <- segment.indicies$index
  names(segment.indicies) <- segment.names
  n.segments <- length(segment.indicies)
  
  read.count <- data$vals
  
  ## parse input parameters
  
  # rho parameter for each segment
  if (length(rho)!=1 | class(rho)!="numeric") stop("rho should be a single numeric value")
  if (rho>1 | rho<0) stop("rho should be between 0 and 1")
  rho <- rep(rho,length(segment.indicies))
  stopifnot(length(segment.indicies)==length(rho))
  names(rho) <- segment.names
  
  # set common component parameters, one per common component
  com.param <- input.parameters[component.type=="common"]
  com.param$component.type <- NULL

  n.specific.components <- nrow(input.parameters[component.type=="specific"])
  n.common.components <- nrow(input.parameters[component.type=="common"])
  
  print(paste(n.segments,"segments, each with",n.common.components,"common component(s), and",n.specific.components,"specific component(s)"))
  
if (algorithm == "EM"){
  output <- EM(segment.indicies, read.count, rho, com.param, n.specific.components, n.common.components, n.segments, input.parameters, max.iterations, break.parameter, quiet, segment.names)
  
} else if (algorithm == "VB"){
  output <- VB(segment.indicies, read.count, rho, com.param, n.specific.components, n.common.components, n.segments, input.parameters, max.iterations, break.parameter, quiet, segment.names)
  
} else { stop("Algorithm should be either EM or VB") }

  return(output)
  
}