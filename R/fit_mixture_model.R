#' Fit shared component gaussian mixed model.
#'
#' \code{fit.model} returns parameters for gassians fitted to the data.
#'
#'
#' @param read.count data.frame with two columns, first being numeric read.counts,
#'  and second being charachter segment names. Column names can be whatever.
#'
#' @param rho.input Numeric value for initial value of rho parameter
#' @param input.parameters
#' @param max.iterations
#' @param break.parameter Numeric value, when the difference in log likelihood 
#' between two iterations is less than this value the iteration will stop
#' @return A list of parameters
#'
#' @examples
#' # see vignette
#' @export
#' @import data.table

fit.model <- function(vals.df, input.parameters, rho.input = 0.5, max.iterations = 40, break.parameter = 5, algorithm = "EM"){

  ### INITIALISATION
  
  if (any(input.parameters$component.type=="common") & any(input.parameters$component.type=="specific")){
  } else {
    stop("At least one common and one specific component required")
  }
  
  if (all(colnames(input.parameters) != "mean")) stop("Mean values required")
  
  ## split input read count data frame into 2: 
  # 1) a list of indexes specifying which read count is in which segment
  # 2) a vector of read counts
  if(ncol(vals.df) != 2) stop("Two columns of input required")
  setnames(vals.df,names(vals.df),c("vals","seg"))
  
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
  print(paste(nrow(input.parameters[component.type=="common"]),"common components"))
  
if (algorithm == "EM"){
  source("R/fit_model_by_EM.R", local = TRUE)
  
} else if (algorithm == "VB"){
  source("R/fit_model_by_VB.R", local = TRUE) 
  
} else { stop("Algorithm should be either EM or VB") }

  return(output)
  
}