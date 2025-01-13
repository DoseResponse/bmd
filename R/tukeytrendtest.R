.tukeytrendtest <- function(y,x){
  ttw <- .tukeytrendfit(y, x)
  res <- multcomp:::summary.glht(multcomp:::glht(model=ttw$mmm, linfct=ttw$mlf))
  
  res
}

