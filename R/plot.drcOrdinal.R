plot.drcOrdinal <- function(object, col_pal = NULL, xlim = NULL){
  if(!require("scales")){
    stop('package "scales" must be installed to plot drcOrdinal object')
  }
  
  if(is.null(col_pal)){
    col_pal <- scales::grey_pal(start = 0.9, end = 0)(length(object$levels))
  }
  
  plotData <- tidyr::pivot_longer(object$data, cols = object$levels) %>% #-c(object$dose, object$weights)) %>% 
    mutate(dose = eval(parse(text=object$dose)),
           cat = factor(name, levels = object$levels),
           prop = value/eval(parse(text=object$weights)))
  
  plot <- ggplot(plotData) +
    geom_col(aes(x = dose, y = prop, fill = cat), alpha = 0.5) +
    scale_fill_manual(breaks = object$levels, values = col_pal[1:length(object$levels)]) +
    lapply(object$levelsMerged, 
           function(levelsMerged){ 
             geom_function(aes(col = levelsMerged), fun = object$drmList[[levelsMerged]]$curve[[1]], data = data.frame(levelsMerged = levelsMerged))}) +
    scale_color_manual(breaks = object$levelsMerged, values = col_pal[2:length(object$levels)]) +
    scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 1, base = exp(1)), limits = xlim) +
    scale_y_continuous(limits = c(0,1)) +
    labs(x = object$dose, y = "proportion") +
    theme_bw()
  
  plot
}
