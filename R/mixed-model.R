library("ggplot2")
library("data.table")

# parse learned/output parameters
parse.output.parameters <- function(output,segment.subset=NULL){

  output.parameters <- data.table(rho=numeric(),w=numeric(),rho_w=numeric(),mu=numeric(),sigma2=numeric(),component.type=character(),segment=character(),iteration=numeric())
  
  for (segment in 1:nrow(output$w.specific)){
    
    output.parameters.temp <- data.frame(
      rho = as.numeric(output$rho[segment]),
      w = c(output$w.common,output$w.specific[segment,]),
      rho_w = c((1-output$rho[segment])*output$w.common,output$rho[segment]*output$w.specific[segment,]),
      mu = c(output$mu.common,output$mu.specific[segment,]),
      sigma2 = c(output$sigma2.common,output$sigma2.specific[segment,]),
      component.type = c(rep("common",length(output$mu.common)),rep("specific",length(output$mu.specific[segment,]))),
      segment = output$segment.names[segment],
      iteration = output$iteration
    )
    
    output.parameters <- rbind(output.parameters,output.parameters.temp)
  }
  
  if(!is.null(segment.subset)){
    output.parameters <- output.parameters[segment %in% segment.subset]
  }

  return(output.parameters)
}

#' Plot fitted shared component gaussian mixed model.
#'
#' \code{plot.components} Returns histograms of fitted gaussians and original data by dataset.
#'
#' @param vals.df original data, identical to fit.model() input
#' @param output.parameters, list of parameters, output of fit.model()
#' @param segment.subset, optional, list of segment names to plot, default uses all segments
#' @param type, optional, if value is "density", a density plot comparing orignal vs model distribution.
#' If value is "QQ", a quantile-quantile plot is output.
#' If left empty the modeled components are plotted (with common and specific in different colours)
#'
#' @examples
#' plot.components(input.data,output.data,segment.subset="segment16")
#' @export
#' @import ggplot2

plot.components <- function(vals.df,output.parameters,segment.subset=NULL,type=NULL){
  output.parameters <- parse.output.parameters(output.parameters,segment.subset)
  
  if(!is.null(segment.subset)){
    vals.df <- vals.df[seg %in% segment.subset]
  }
  
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
  
  if (is.null(type)){
  ggplot(temp[source=="Inferred"], aes(vals, colour=component.type)) +
    geom_freqpoly(binwidth=0.1) +
    ggtitle(paste("iteration:",unique(output.parameters$iteration))) +
    scale_colour_manual(values = c("blue","red","black")) +
    theme(legend.position = "bottom") +
    facet_wrap(~seg,scales="free",nrow=1) # source~seg for original data underneath

  } else if (type=="density"){
  ggplot(temp, aes(vals, colour=source)) +
    geom_density() +
    scale_colour_manual(values = c("red","black")) +
    ggtitle(paste("iteration:",unique(output.parameters$iteration))) +
    theme(legend.position = "bottom") +
    facet_wrap(~seg,scales="free",nrow=1)
    # to plot all components individually
    # geom_density(data=temp[source=="Inferred"],aes(vals,group=comp)) +
  }else if (type=="QQ"){
    
    qq <- data.frame(x=numeric(),y=numeric(),segment=character())
    for (segment in segment.subset){
      d <- as.data.frame(qqplot(temp[seg==segment & source=="Inferred"]$vals,temp[seg==segment & source=="Original Data"]$vals, plot.it=FALSE))
      d$segment <- segment
      qq <- rbind(qq,d)
    }
    qq$segment <- as.factor(as.numeric(qq$segment))
    
    ggplot(qq,aes(x, y)) +
      geom_point(size=0.5) +
      geom_abline(intercept=0,slope=1,colour='red') +
      facet_wrap(~segment,scales="free",nrow=1)
  }
}

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
#' @param init.max
#' @param break.parameter Numeric value, when the difference in log likelihood 
#' between two iterations is less than this value the iteration will stop
#' @return A list of parameters
#'
#' @examples
#' # see vignette
#' @export
#' @import data.table

fit.model <- function(vals.df,rho.input,input.parameters,init.max=40,break.parameter=50){

### INITIALISATION

  if (any(input.parameters$component.type=="common") & any(input.parameters$component.type=="specific")){
  } else {
    stop("At least one common and one specific component required")
  }

## split input read count data frame into 2: 
# 1) a list of indexes specifying which read count is in which segment
# 2) a vector of read counts
stopifnot(ncol(vals.df)==2)
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


# initialise segment specific parameter holding matricies
mu.k <- sigma2.k <- w.k <- matrix(nrow=length(read.count),ncol=n.specific.components)
# although these could be a simple j*k matrix, we're already iterating through segments in
# the M step so may as well create the whole matrix required by the dnorm E step below 

for (row in 1:n.specific.components){
  w.k[,row] <- input.parameters[component.type=="specific"][row]$w
  mu.k[,row] <- input.parameters[component.type=="specific"][row]$mu
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
for (round in 1:init.max){
#### E Updates

# for each common component, cacluate prob for each data point, outputs i x c matrix
common <- apply(com.param, 1, function(params){ params['w'] * dnorm(read.count, params['mu'], params['sigma2']^0.5)})

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

com.param$mu <- colSums(phi_ic * psi_i * read.count) / topsum

com.param$sigma2 <- colSums(phi_ic * psi_i * outer(read.count,com.param$mu, FUN="-")^2) / topsum

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
output <- list(w.specific=w.k.u,mu.specific=mu.k.u,sigma2.specific=sigma2.k.u,w.common=com.param$w,mu.common=com.param$mu,sigma2.common=com.param$sigma2,rho=unname(rho),iteration=iter.count,likelihood=likelihood.byiter,segment.names=segment.names)

return(output)

}# end of fit.model function


# for debugging enter data manually
# vals.df <- chr6.input
# input.parameters <- initial.parameters.realistic
# rho.input <- 0.5

