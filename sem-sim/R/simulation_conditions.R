simulation_conditions = \() {
  num_rep = 5
  num_batch = 40 # multiplies num_rep
  num_obs = 75
  num_precon = c(5, 6, 7, 8, 9, 10)
  formulas = list(
    gdsem = # simplified mediation model
     "bf(x1 | mi() ~ 0) +
      bf(y1i1 ~ mi(x1)) +
      bf(y1i2 ~ mi(x1)) +
      bf(y1i3 ~ mi(x1)) +
      bf(x2 | mi() ~ 0 + mi(x1)) +
      bf(y2i1 ~ mi(x2)) +
      bf(y2i2 ~ mi(x2)) +
      bf(y2i3 ~ mi(x2)) +
      bf(x3 | mi() ~ 0 + mi(x1)) +
      bf(y3i1 ~ mi(x3)) +
      bf(y3i2 ~ mi(x3)) +
      bf(y3i3 ~ mi(x3)) +
      bf(x4 | mi() ~ 0 + mi(x1) + mi(x2) + mi(x3),
        sigma ~ mi(x1) + mi(x2) + mi(x3)) +
      bf(y4i1 ~ mi(x4)) +
      bf(y4i2 ~ mi(x4)) +
      bf(y4i3 ~ mi(x4))",
    poldem ='
    # latent variable definitions
      ind60 =~ x1 + fl2ind*x2 + fl3ind*x3
      dem60 =~ y1 + fl2dem1*y2 + fl3dem1*y3 + fl4dem1*y4
      dem65 =~ y5 + fl2dem2*y6 + fl3dem2*y7 + fl4dem2*y8
    # regressions
      dem60 ~ b1dem1*ind60
      dem65 ~ b1dem2*ind60 + b2dem2*dem60
    # residual correlations
      y1 ~~ cory1y5*y5
      y2 ~~ cory2y4*y4 + cory2y6*y6
      y3 ~~ cory3y7*y7
      y4 ~~ cory4y8*y8
      y6 ~~ cory6y8*y8
    # variances
      dem60 ~~ sigmadem1*dem60
      dem65 ~~ sigmadem2*dem65
      ind60 ~~ sigmaind*ind60
      x1 ~~ tauind1*x1
      x2 ~~ tauind2*x2
      x3 ~~ tauind3*x3
      y1 ~~ tau1dem1*y1
      y2 ~~ tau2dem1*y2
      y3 ~~ tau3dem1*y3
      y4 ~~ tau4dem1*y4
      y5 ~~ tau1dem2*y5
      y6 ~~ tau2dem2*y6
      y7 ~~ tau3dem2*y7
      y8 ~~ tau4dem2*y8
'
  )
  ifs_priors = list(
    poldem = NULL,
    gdsem = NULL,
    # replace yourself with prior so targets wont see diff:
    gdsem = \(p) {
      p[p$source == "item_incpt",1] = "normal(0,0.01)"
      p[p$source == "item_load",1] = "constant(1)"
      p[p$source == "coef_mean",1] = "constant(0)"
      p[p$source == "fact_sigma",1] = "constant(1)"
      p[p$source == "fact_dsigma",1] = "constant(0)"
      p[p$source == "coef_sigma",1] = "constant(0)"
      p[p$source == "item_sigma",1] = "constant(1)"
      p 
    }
  )
  # Setting up the dataframe ----
  utils.named_tibble(formulas) %>%
    mutate(
      data_template = list(
        utils.make_template(
          c(paste0("x", 1:4), paste0(rep(paste0("y", 1:4), each = 3) ,"i", 1:3))
        ),
        utils.make_template(
          c(paste0("x", 1:3), paste0("y", 1:8))
        )
      ),
      priors = list(
       veryweak = {p <- brms::get_prior(eval(parse(text = formulas[[1]])) +
          set_rescor(F), data_template[[1]])
        p <- p[!(p$class == "b" & p$coef == "") & p$resp != "", ]
        p[grepl("^x", p$resp) & grepl("Int", p$class), "source"] <- "fact_incpt"
        p[grepl("^x", p$resp) & grepl("sig", p$class), "source"] <- "fact_sigma"
        p[grepl("^x", p$resp) & p$class == "b", "source"] <- "coef_mean"
        p[grepl("sig", p$dpar) & grepl("Int", p$class), "source"] <- "fact_dsigma"
        p[grepl("sig", p$dpar) & p$class == "b", "source"] <- "coef_sigma"
        p[grepl("^y", p$resp) & grepl("mix", p$coef), "source"] <- "item_load"
        p[grepl("^y", p$resp) & grepl("Int", p$class), "source"] <- "item_incpt"
        p[grepl("^y", p$resp) & grepl("sig", p$class), "source"] <- "item_sigma"
        p[grepl(".+i1$", p$resp) & p$source == "item_load", "source"] <- "first_load"
        p[p$source == "item_incpt",1] = "normal(0,32)"
        p[p$source == "fact_incpt",1] = "normal(0,10)"
        p[p$source == "item_load",1] = "normal(0,10)"
        p[p$source == "coef_mean",1] = "normal(0,10)"
        p[p$source == "fact_sigma",1] = "gamma(1,.5)"
        p[p$source == "fact_dsigma",1] = "expgamma(1,.5)"
        p[p$source == "coef_sigma",1] = "normal(0,2.5)"
        p[p$source == "item_sigma",1] = "gamma(1,.5)"
        p[p$source == "first_load",1] = "constant(1)"
        p},
        blavaan = c(
          nu = "normal(0,32)", alpha = "normal(0,10)",
          lambda = "normal(0,10)", beta = "normal(0,10)",
          theta = "gamma(1,.5)[sd]", psi = "gamma(1,.5)[sd]",
          rho = "beta(1,1)"
        ))
      ) %>%
    left_join(
      utils.named_tibble(ifs_priors) %>%
        utils.cross_vec(num_precon) %>%
        mutate(num_precon = num_precon*sapply(ifs_priors, length)) %>%
        unique,
      by = "id") %>%
    mutate(
      ifs_priors = mapply(priors, ifs_priors, FUN = \(p, f) {
        if(!is.null(f)) return(f(p))
    })) %>%
  # Cross over the number of simulations and SBC reps
    utils.cross_vec(num_obs) %>%
    utils.cross_vec(num_rep) %>%
  # Batching: repeats reps
    slice(
      rep(1:n(), each = num_batch)
    ) %>%
    mutate(batch = 1:n(), .by = id)
}