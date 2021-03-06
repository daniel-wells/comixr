#' Convert fit_comixture output into tidy data format
#'
#' \code{tidy_output} Convert fit_comixture output into tidy data format .e.g for plotting. This function parses learned/output parameters - the output is given as a nested list, but we need a data frame with one line per component (for each segment) to plot with ggplot2
#'
#' @param output the output given by fit_comixture
#' @param segment.subset, optional, list of segment names to plot, default uses all segments
#'
#' @examples
#' 
#' @export
#' @import data.table

tidy_output <- function(output, segment.subset = NULL){

  output.parameters <- data.table(rho = numeric(),
                                  w = numeric(),
                                  rho_w = numeric(),
                                  mu = numeric(),
                                  variance = numeric(),
                                  component.type = character(),
                                  segment = character(),
                                  iteration = numeric())
  
  for (segment in 1:nrow(output$specific_parameters$mix_weights)){
    
    output.parameters.temp <- data.frame(
      rho = as.numeric(output$rho[segment]),
      w = c(output$common_parameters$mix_weights, output$specific_parameters$mix_weights[segment, ]),
      rho_w = c((1 - output$rho[segment]) * output$common_parameters$mix_weights, output$rho[segment] * output$specific_parameters$mix_weights[segment, ]),
      mu = c(output$common_parameters$mean, output$specific_parameters$mean[segment, ]),
      variance = c(output$common_parameters$variance, output$specific_parameters$variance[segment, ]),
      component.type = c(rep("common", length(output$common_parameters$mean)), rep("specific", length(output$specific_parameters$mean[segment, ]))),
      segment = output$segment.names[segment],
      iteration = output$iteration
    )
    
    output.parameters <- rbind(output.parameters, output.parameters.temp)
  }
  
  if(!is.null(segment.subset)){
    output.parameters <- output.parameters[segment %in% segment.subset]
  }

  return(output.parameters)
}

#' Plot fitted shared component gaussian mixed model.
#'
#' \code{plot_comixture} Returns histograms of fitted gaussians and original data by dataset.
#'
#' @param data original data, identical to fit_comixture() input
#' @param output.parameters, list of parameters, output of fit_comixture()
#' @param segment.subset, optional, list of segment names to plot, default uses all segments
#' @param type, optional, if value is "density", a density plot comparing orignal vs model distribution.
#' If value is "QQ", a quantile-quantile plot is output.
#' If left empty the modeled components are plotted (with common and specific in different colours)
#'
#' @examples
#' # plot_comixture(input.data, output.data, segment.subset = "segment16")
#' @export
#' @import ggplot2 data.table

plot_comixture <- function(data, output.parameters, segment.subset = NULL, type = NULL){
  
  output.parameters <- tidy_output(output.parameters, segment.subset)
  
  if(!is.null(segment.subset)){
    data <- data[seg %in% segment.subset]
  }
  
  # create new data based on fitted model parameters
  temp <- data.table(vals = numeric(),
                     comp = character(),
                     seg = character(),
                     component.type = character(),
                     source = character())
  
  for (component in 1:nrow(output.parameters)){
    temp <- rbind(temp,data.table(vals=c(rnorm(70000 * output.parameters[component]$rho_w, 
                                            mean = output.parameters[component]$mu, 
                                            sd = output.parameters[component]$variance^0.5), 0),
                                            comp = component,
                                            seg = output.parameters[component]$segment,
                                            component.type = output.parameters[component]$component.type,
                                            source = "Inferred"))
  }
  
  # fill in extra columns to enable binding below
  data$component.type <- "Original Data"
  data$source <- "Original Data"
  if (!any(names(data) == "comp")){
    data$comp <- NA
  }
  
  # combine original data with new data based from model
  temp <- rbind(temp, data)
  
  # ggplot(temp, aes(vals, group = comp,colour=component.type)) +
  #   geom_freqpoly(binwidth=0.1) +
  #   ggtitle(paste("iteration:",iter.count)) +
  #   scale_colour_manual(values = c("blue","red","black")) +
  #   facet_wrap(~seg,nrow = 2)
  
  if (is.null(type)){ # plot freq poly of fitted model components
                      # one line for common one for specific
  ggplot(temp[source == "Inferred"], aes(vals, colour = component.type)) +
    geom_freqpoly(binwidth = max(temp$vals) / 100) +
    ggtitle(paste("iteration:", unique(output.parameters$iteration))) +
    scale_colour_manual(values = c("blue","red","black")) +
    theme(legend.position = "bottom") +
    facet_wrap(~seg, scales = "free", nrow = 1) # source~seg for original data underneath

  } else if (type == "density"){ # plot total original density vs fitted density
  ggplot(temp, aes(vals, colour = source)) +
    geom_density() +
    scale_colour_manual(values = c("red", "black")) +
    ggtitle(paste("iteration:", unique(output.parameters$iteration))) +
    theme(legend.position = "bottom") +
    facet_wrap(~seg, scales = "free", nrow = 1)
    # to plot all components individually
    # geom_density(data=temp[source=="Inferred"],aes(vals,group=comp)) +
    
  } else if (type == "QQ"){ # plot QQ plot of original vs fitted model quantiles
    
    qq <- data.frame(x = numeric(), y = numeric(), segment = character())
    for (segment in output.parameters$segment){
      d <- as.data.frame(qqplot(temp[seg == segment & source == "Inferred"]$vals, temp[seg == segment & source == "Original Data"]$vals, plot.it = FALSE))
      d$segment <- segment
      qq <- rbind(qq, d)
    }
    qq$segment <- as.factor(as.numeric(qq$segment))
    
    ggplot(qq,aes(x, y)) +
      geom_point(size = 0.5) +
      geom_abline(intercept = 0, slope = 1, colour = 'red') +
      facet_wrap(~segment, scales="free", nrow = 1)
  }
}