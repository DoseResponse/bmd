drmHetVarSelfStarter <- function(formula, var.formula, data, fct) 
{
  dName <- as.character(formula)[3]
  object0 <- drm(formula, data = data, fct = fct)
  
  data <- cbind(data, fitted = fitted(object0), residuals = residuals(object0))
  data.agg <- dplyr::summarise(dplyr::group_by(data, fitted), 
                               dose0 = mean(.data[[dName]]), 
                               sigma0 = sqrt(mean(residuals^2)))
  colnames(data.agg)[2] <- dName
  # formula <- as.formula(var.formula)
  formula0 <- as.formula(paste0("sigma0 ", paste0(as.character(var.formula), collapse = " ")))
  sigmaMod <- lm(formula0, data = data.agg)
  
  startList <- list(curvePar = coef(object0),
                    tauPar = coef(sigmaMod))
  startList
}

