
#' @title S3 method
#' @export
plot.drcOrdinal <- function(x, ...){
  object <- x
  dots <- list(...)
  col_pal  <- dots$col_pal
  xlim <- dots$xlim
  
  if(!requireNamespace("scales")){
    stop('package "scales" must be installed to plot drcOrdinal object')
  }
  if(!requireNamespace("dplyr")){
    stop('package "dplyr" must be installed to plot drcOrdinal object')
  }
  if(!requireNamespace("ggplot2")){
    stop('package "ggplot2" must be installed to plot drcOrdinal object')
  }
  
  if(is.null(col_pal)){
    col_pal <- scales::grey_pal(start = 0.9, end = 0)(length(object$levels))
  }
  
  plotData <- tidyr::pivot_longer(object$data, cols = object$levels) %>% #-c(object$dose, object$weights)) %>% 
    dplyr::mutate(dose = eval(parse(text=object$dose)),
           cat = factor(name, levels = object$levels),
           prop = value/eval(parse(text=object$weights)))
  
  plot <- ggplot2::ggplot(plotData) +
    ggplot2::geom_col(aes(x = dose, y = prop, fill = cat), alpha = 0.5) +
    ggplot2::scale_fill_manual(breaks = object$levels, values = col_pal[1:length(object$levels)]) +
    lapply(object$levelsMerged, 
           function(levelsMerged){ 
             ggplot2::geom_function(aes(col = levelsMerged), fun = object$drmList[[levelsMerged]]$curve[[1]], data = data.frame(levelsMerged = levelsMerged))}) +
    ggplot2::scale_color_manual(breaks = object$levelsMerged, values = col_pal[2:length(object$levels)]) +
    ggplot2::scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 1, base = exp(1)), limits = xlim) +
    ggplot2::scale_y_continuous(limits = c(0,1)) +
    ggplot2::labs(x = object$dose, y = "proportion") +
    ggplot2::theme_bw()
  
  plot
}
