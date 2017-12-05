bmd<-function (object, bmr, backg = 0, def = c("excess", "additional", 
    "relative", "hybrid"), interval = c("delta"), ma = FALSE, 
    maList = NULL, display = FALSE) 
{
    def <- match.arg(def)
    respType <- object$type
    if (identical(respType, "binomial")) {
        bmrScaled <- switch(def, excess = bmr * (1 - backg) + 
            backg, additional = bmr + backg)
        typeVal <- "absolute"
    }
    if (identical(respType, "continuous")) {
        bmrScaled0 <- 100 * switch(def, relative = bmr, hybrid = sqrt(summary(object)$resVar) * 
            (qnorm(1 - backg) - qnorm(1 - (backg + bmr)))/(abs(diff(predict(object, 
            data.frame(c(0, Inf)))))))
        bmrScaled <- as.numeric(format(bmrScaled0, digits = 3))
        typeVal <- "relative"
    }
    if (display) {
        cat("Effective relative response level: ", bmrScaled, 
            "\n\n")
    }
    if (ma) {
        if (identical(object$type, "binomial")) {
            maList <- list(LL.2(), LN.2(), W1.2(), W2.2())
        }
        if (identical(object$type, "continuous")) {
            maList <- list(LL.5(), LN.4(), W1.4(), W2.4(), FPL.4(-1, 
                1), FPL.4(-2, 3), FPL.4(-0.5, 0.5))
        }
    }
    if (!is.null(maList)) {
        resMat <- maED(object, maList, bmrScaled, interval = "buckland", 
            level = 0.9, type = typeVal, display = display)[, 
            c(1, 3), drop = FALSE]
    }
    else {
        resMat <- ED(object, bmrScaled, interval = interval, 
            level = 0.9, type = typeVal, display = display)[, 
            c(1, 3), drop = FALSE]
    }
    colnames(resMat) <- c("BMD", "BMDL")
    if (display) {
        cat("\n\n")
    }
    resMat
}