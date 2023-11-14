plot.drcOrdinal <- function(object, col_pal = NULL, xlim = NULL){
  if(is.null(col_pal)){
    col_pal <- scales::grey_pal(start = 0.9, end = 0)(length(object$categories))
  }
  
  plotData <- pivot_longer(object$data, cols = object$categories) %>% #-c(object$dose, object$weights)) %>% 
    mutate(dose = eval(parse(text=object$dose)),
           cat = factor(name, levels = object$categories),
           prop = value/eval(parse(text=object$weights)))
  
  plot <- ggplot(plotData) +
    geom_col(aes(x = dose, y = prop, fill = cat), alpha = 0.5) +
    scale_fill_manual(breaks = object$categories, values = col_pal[1:length(object$categories)]) +
    lapply(object$catMerged, 
           function(catMerged){ 
             geom_function(aes(col = catMerged), fun = object$drmList[[catMerged]]$curve[[1]])}) +
    scale_color_manual(breaks = object$catMerged, values = col_pal[2:length(object$categories)]) +
    scale_x_continuous(trans = scales:::pseudo_log_trans(sigma = 1, base = exp(1)), limits = xlim) +
    scale_y_continuous(limits = c(0,1)) +
    labs(x = object$dose, y = "proportion") +
    theme_bw()
  
  plot
}
