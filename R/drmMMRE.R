drmMMRE <- function(formula, exp_id, data, fct, type = c("continuous", "binomial")){
  call_expr <- match.call()
  
  # assertions
  if (!inherits(formula, "formula")) {
    stop("Argument 'formula' must be a formula.")
  }
  
  if(missing(exp_id)){
    stop("supply hierarchical structure with argument \"exp_id\"")
  }
  # handling of exp_id variable
  exp_id_sym <- substitute(exp_id)
  exp_id_char <- as.character(exp_id_sym)
  exp_id_unique <- unique(data[[exp_id_char]])
  
  if(missing(data)){
    stop("\"data\" not supplied")
  }
  
  type <- match.arg(type)
  
  if(!requireNamespace("metafor")){
    stop("install package metafor by running\nremotes::install_github(\"wviechtb/metafor\")")
  }
  
  if(!requireNamespace("Matrix")){
    stop("install package \"Matrix\"")
  }
  
  # Fit models for each exp_id
  drc.objList <- lapply(exp_id_unique,
                        function(exp_id){
                          drm(formula, data = subset(data, data[[exp_id_char]] == exp_id), fct = fct, type = type)
                          # eval(substitute(drm(formula = formula0, data = subset(data, data[[exp_id_char]] == x), fct = fct0, type = type0),
                          #            list(formula0 = formula,
                          #                 data0 = subset(data, data[[exp_id_char]] == x),
                          #                 fct0 = fct,
                          #                 type0 = type)))
                        })
  names(drc.objList) <- exp_id_unique
  
  block_diag_matrix <- as.matrix(Matrix::bdiag(lapply(drc.objList, vcov)))
  
  MV_Meta_Data <- do.call(rbind, lapply(drc.objList, 
                                        function(object){
                                          data.frame(Estimate = coef(object),
                                                     Coef = substr(names(coef(object)), 1, 1),
                                                     exp_id = object$origData[[exp_id_char]][1])
                                          }
                                        ))
  
  # Fit MA_MA_model
  MV_MA_model<-tryCatch( expr={metafor::rma.mv(Estimate
                                      ,(block_diag_matrix+t(block_diag_matrix))/2
                                      , mods = ~0+Coef
                                      ,random= ~ 0+Coef|exp_id
                                      , data=MV_Meta_Data
                                      ,struct = "UN"
                                      ,control=list(rel.tol=1e-8))}
                         , error = function(e){
                           cat("MV_MA_Model_failed: ", e$message, "\n")
                           return(NA)})
  
  # Save coefficients from meta model
  par <- coef(MV_MA_model)
  
  # Create drc object and replace selected elements
  object <- drm(formula = formula, data = data, fct = fct, type = type, control = drmc(errorm = FALSE, maxIt = 1, noMessage = TRUE))
  object$fit$par <- par
  coefNames <- names(object$coefficients)
  object$coefficients <- par
  names(object$coefficients) <- coefNames
  object$curve[[1]] <- object$pfFct(parmMat = matrix(par, nrow = 1))
  object$parmMat <- matrix(par, nrow = length(par))
  colnames(object$parmMat) <- 1
  object$pmFct <- function() t(object$parmMat)
  
  object$data$exp_id <- data[[exp_id_char]]
  object$dataList$exp_id <- data[[exp_id_char]]
  
  object$names$exp_idName <- exp_id_char
  object$names$exp_idNames <- exp_id_unique
  
  object$fit <- object$fit[c("par")]
  object$summary <- NULL
  object$start <- NULL
  for(exp_i in exp_id_unique){
    object$predres[data[[exp_id_char]] == exp_i,] <- drc.objList[[exp_i]]$predres
  }
  object$call <- call_expr
  # 
  # object$df.residual <- NA
  # object$sumList$df.residual <- NA
  # object$deriv1 <- NULL
  
  # add MV_MA_model and drc.objList
  object$objList <- drc.objList
  object$MV_MA_model <- MV_MA_model
  
  # assign new class and return object
  class(object) <- c("drcMMRE", "drc") # defined as subclass of drc class
  object
}


vcov.drcMMRE <- function(object, ...){
  vcov(object$MV_MA_model)
}
