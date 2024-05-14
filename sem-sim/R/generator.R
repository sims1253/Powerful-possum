generator.get_priors = function(f, p) {
  # ULI identification constraint
  p[p$source == "fact_incpt",1] <- "constant(0)"
  p[p$source == "first_load",1] <- "constant(1)"
  # We will not estimate item intercepts
  p[p$source == "item_incpt",1] <- "constant(0)"
  if(identical(f, "sem01"))
    list(
    wide = (\(){
    p[p$source == "fact_sigma",1] <- "gamma(1,1)"
    p[p$source == "item_load",1] <- "normal(0,10)"
    p[p$source == "item_sigma",1] <- "gamma(0.5,1)"
    p[p$source == "coef_mean",1] <- "normal(0,10)"
    p})()
    )
  else
    stop("No match found")
}