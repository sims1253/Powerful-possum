# Preprocessing -----
utils.named_tibble = function(x) {
  objname <- deparse(substitute(x))
  eval(parse(text = paste0(
    "tibble(id = names(x), ", objname, " = unname(x))"
  )))
}

# cross a vector's values with all of the rows in d
utils.cross_vec = function(d, v) {
  objname <- deparse(substitute(v))
  eval(parse(text = paste0(
    "d %>%
    bind_cols(tibble(", objname, " = list(v))) %>%
    unnest(", objname, ")"
  )))
}

# get response variables in brms formula
utils.get_resp = function(x, pattern) {
  if(inherits(x, "bform"))
    r = x$responses
  if(inherits(x, "brmsprior"))
    r = x$resp
  unique(r[grepl(pattern, r)])
}

# make template dataframe from names
utils.make_template = \(v) {
  setNames(as.data.frame(matrix(0, ncol = length(v))), v) 
}

# generate data with additional handling:
# - set vars matching na_pattern to NA
# - cache the result and return filepath
utils.make_data = \(condition) {
  stopifnot(identical(nrow(condition), 1L))
  filepath = with(condition, {
    if(id %in% "gdsem") {
      formula = eval(parse(text = formulas[[1]])) +
        set_rescor(FALSE)
      if(num_precon == 0) {
        generator = SBC::SBC_generator_brms(
          formula, prior = priors[[1]],
          data = data_template[[1]][rep(1, num_obs), ],
          stanvars = bdlvm_stanvars, generate_lp = TRUE,
          thin = 50, warmup = 10000, refresh = 0
        )
        d = SBC::generate_datasets(generator, num_rep)
        d = utils.format_sbc_data(d, formula,
          na_pattern = "^x[1-4]$", drop = "^$")
      } else {
        generator = brm(
          formula, prior = ifs_priors[[1]],
          data = data_template[[1]][rep(1, num_precon), ],
          warmup = 10000, iter = 10001, chains = 1, refresh = 0,
          backend = "cmdstanr", sample_prior = "only",
          stanvars = bdlvm_stanvars
        )
        ifs_data = SBC:::brms_full_ppred(generator)[[1]]
        # for each dataset fit and sample
        ifs_generator = brm(
          formula, prior = priors[[1]], data = ifs_data,
          warmup = 10000, iter = 10000 + 50*num_rep, thin = 50,
          chains = 1, refresh = 0, backend = "cmdstanr",
          stanvars = bdlvm_stanvars
        )
        generated = SBC:::brms_full_ppred(ifs_generator,
          newdata = data_template[[1]][rep(1, num_obs), ])
        d = structure(list(
          variables = as_draws_matrix(ifs_generator),
          generated = generated
        ), class = "SBC_datasets")
        d = utils.format_sbc_data(d, formula,
          na_pattern = "^x[1-4]$", drop = "^$")
        d$generated = lapply(d$generated, bind_rows, ifs_data)
        d
      }
    } else if(id %in% "poldem") {
      generator = bsem(
        formulas[[1]], dp = priors[[1]],
        data = apply(data_template[[1]][rep(1, num_obs), ],
          c(1,2), \(x) rnorm(1,0)),
        target = "stan", prisamp = TRUE,
        n.chains = 1, burnin = 10000, sample = 50*num_rep,
        bcontrol = list(thin = 50)
      )
      generated = sampleData(generator)
      d = list(
        variables = utils.format_blav_draws(generator),
        generated = lapply(
          generated,
          \(x) as.data.frame(x[[1]]) |>
            setNames(colnames(lavaan::PoliticalDemocracy))
        )
      )
    } else {
      stop("No match found")
    }
      filepath = here(
        tar_path_store(), "user",
        paste0(
          id, "_n", num_obs, "_", names(priors),
          ".ifs", num_precon, "_batch", batch, ".data")
      )
      qsave(d, filepath)
      filepath
  })
  condition %>% mutate(data = filepath)
}

utils.format_sbc_data = \(d, formula, ..., na_pattern, drop) {
  data_ok = sapply(d$generated, \(x) all(!is.na(x)))
  # drop simulations with any NA/NaN values
  d$variables = d$variables[data_ok,]
  d$generated = d$generated[data_ok]
  # drop user-specified columns
  d$variables = d$variables[, !grepl(drop, colnames(d$variables))]
  # put missing data where needed
  make_na = utils.get_resp(formula, na_pattern)
  if(!is.null(make_na)) {
    d$na_values = lapply(d$generated,\(x)x[, make_na])
    d$generated = lapply(d$generated,\(x){
      x[, make_na] = NA_real_
      x
    })
  }
  return(d)
}

utils.format_blav_draws = \(blavobj) {
  samples = blavInspect(blavobj, "mcmc")
  as_draws_matrix(bind_rows( # puts everything in 1 chain, hopefully ok
    lapply(samples, \(x) as.data.frame(x[, !duplicated(colnames(x))]))
  ))
}

# Fitting -----
utils.sbc_compute = function(x, ...) {
  if("gdsem" == x$id) {
    return(utils.compute_gdsem(x, ...))
  } else if("poldem" == x$id) {
    return(utils.compute_poldem(x, ...))
  } else {
    stop("No compute step defined for ", x)
  }
}

utils.compute_gdsem = function(condition, ...) {
  with(condition, {
    formula = eval(parse(text = formulas[[1]])) +
      set_rescor(FALSE)
    backend = SBC_backend_brms(
      formula,
      template_data = data_template[[1]],
      prior = priors[[1]],
      stanvars = bdlvm_stanvars,
      chains = 4, thin = 1, warmup = 1000, iter = 2000
    )
    data = qread(data[[1]])
    SBC::compute_SBC(data, backend, ...)
  })
}

utils.compute_poldem = function(condition, ...) {
  with(condition, {
    data = qread(data[[1]])
    fits = list()
    stats = list()
    for(i in seq_along(data$generated)) {
      singlefit = bsem(
        formulas[[1]], dp = priors[[1]],
        data = data$generated[[i]],
        target = "stan", prisamp = TRUE,
        n.chains = 4, burnin = 1000, sample = 1000,
        bcontrol = list(thin = 1)
      )
      singlestat = SBC::SBC_statistics_from_single_fit(
        utils.format_blav_draws(singlefit),
        data$variables[i,],
        NULL, 10, 2, NULL, NULL
      )
      singlestat = bind_cols(sim_id = i, singlestat)
      stats = c(stats, list(singlestat))
      fits = c(fits, list(singlefit))
    }
    return(list(
      stats = bind_rows(stats),
      fits = fits
    ))
  })
}

# Postprocessing ----
## Extract summaries ----
utils.extract_stats = \(each_model, full_results) {
  bind_cols(
    each_model[, c("id", "batch", "num_obs", "num_precon")],
    full_results$stats %>%
    filter(
      !grepl("bsp_y.i1", variable),
      !variable %in% c("lprior", "lp__"),
      !grepl("^b_.+_Intercept$", variable)
    )
  ) %>% utils.alt_names
}
utils.extract_backend = \(each_model, full_results) {
  x = left_join(
    full_results$backend_diagnostics,
    full_results$stats %>%
      filter(
        !grepl("bsp_y.i1", variable),
        !variable %in% c("lprior", "lp__")
      ) %>%
      summarise(
        mean_rhat = mean(rhat, na.rm = TRUE),
        .by = sim_id
      ),
    by = "sim_id"
  )
  bind_cols(each_model[, c("id", "batch", "num_obs", "num_precon")], x)
}
### ^Helpers ----
utils.alt_names = function(df) {
  mutate(df,
    param_type = case_when(
      grepl("^Intercept_y.+i[0-9]+$", variable) ~
        "Item intercept",
      grepl("^bsp_y.+_mix[0-9]+$", variable) ~
        "Factor loading",
      grepl("^sigma_y[0-9]+i[0-9]+$", variable) ~
        "Item std. dev.",
      grepl("^(sigma_x[0-9]+)|(Intercept_sigma_x[0-9]+)$", variable) ~
        "Factor std. dev.",
      grepl("^bsp_.*x[0-9]+_.*mix.+$", variable) ~
        "Factor slopes",
      # grepl("^bsp_x[0-9]+_.*mix.+$", variable) ~
      #   "Slope on mean",
      # grepl("^bsp_sigma_x[0-9]+_.*mix.+$", variable) ~
      #   "Slope on std. dev."
    ) %>%
    factor(levels = c(
      "Factor slopes",
      #"Slope on mean",
      #"Factor std. dev.",
      #"Slope on std. dev.",
      "Item intercept",
      "Factor loading",
      "Item std. dev."
    ), ordered = TRUE)
  )
}

# Plotting ----
utils.scholzplot = function(df) {
  df %>%
    summarise(
      log_gamma = log(utils.gammadisc(rank, M = 399)),
      .by = c(num_precon, variable)
    ) %>%
    mutate(
      isifs = ifelse(num_precon == 0,
        "Original", "Implicit"),
      num_precon = factor(num_precon,
        levels = unique(num_precon), ordered = TRUE),
      log_gamma = ifelse(is.infinite(log_gamma),
        -1000, log_gamma)
    ) %>%
    ggplot(aes(x = num_precon, y = log_gamma)) +
    stat_pointinterval(aes(color = isifs),
      position = position_dodge(width = 0.6),
      .width = c(0.66, 0.90),
      fatten_point = 1) +
    scale_size_continuous(range = c(1, 4))+
    theme_bw() +
    geom_hline(yintercept = log(SBC:::adjust_gamma(N = 399, L = 1))) +
    scale_y_continuous(trans = "pseudo_log",
      breaks = c(0, -10, -100, -1000),
      labels = c(0, -10, -100, "≤ -1000")) +
    theme_bw(base_size = 12) +
    theme(
      strip.text.x = element_text(size = 9),
      strip.text.y = element_text(size = 9),
      axis.ticks.y = element_blank()
    ) +
    xlab("Preconditioning Sample Size") +
    ylab("log(γ) score") +
    labs(color = "Prior") +
    scale_color_manual(values = c("#4477AA", "#EE6677"))
}

utils.rhatplot = \(df, opt) {
  df %>%
    complete(num_precon, sim_id) %>%
    mutate(
      num_precon = factor(num_precon,
        levels = unique(num_precon), ordered = TRUE),
      cat_fit = factor(case_when(
        mean_rhat < 1.01 ~ "[1, 1.01)",
        mean_rhat < 1.05 ~ "[1.01, 1.05)",
        mean_rhat < 1.5 ~ "[1.05, 1.50)",
        !is.na(mean_rhat) ~ "1.50+",
        TRUE ~ "Failed fit"
      ), levels = c("Failed fit", "1.50+", "[1.05, 1.50)",
        "[1.01, 1.05)", "[1, 1.01)"))
    ) %>%
    ggplot(aes(x = num_precon, fill = cat_fit)) + 
    geom_bar(position = "fill") +
    scale_fill_viridis_d(
      name = "Rhat",
      direction = -1,
      option = opt
    ) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme_bw() + labs(y = "Proportion", x = "Preconditioning sample size")
}

utils.sbc_ecdf = function(df) {
  nifs = unique(df$num_precon)
  grouping = df[, c("variable", "param_type")] %>%
    unique %>%
    summarise(
      param_type = setNames(list(variable), unique(param_type)),
      .by = param_type
    ) %>% pull(param_type)
  plot_list = list()
  for(n in nifs) {
    title_string = ifelse(n == 0, "No preconditioning", "Preconditioning")
    subtitle_string = ifelse(n == 0, "", paste0("(n = ",n,")"))
    p = df %>% filter(num_precon == n) %>%
    SBC::plot_ecdf_diff(combine_variables = grouping) +
    facet_wrap(facets = "group", ncol = 1, nrow = 2,
      scales = "free") +
    ggtitle(title_string, subtitle_string) +
    theme_bw() + theme(
      #legend.position = "none",
      axis.text.x =
        element_text(angle = 35, size = 8, hjust = 1)
    )
    plot_list = c(plot_list, list(p))
  }
  wrap_plots(plot_list, ncol = length(nifs), nrow = 1) +
    plot_layout(guides = "collect")
}
## Helpers ^ -----
utils.gammadisc = function (ranks, M) {
  if (any(is.na(ranks))) return(NA)
  iset = 1:(M + 1)
  S = length(ranks)
  R_i <- sapply(iset, \(i) sum(ranks <  i))
  z_i <- sapply(iset, \(i) i/(M + 1))
  x1 <- pbinom(q = R_i, size = S, prob = z_i)
  x2 <- 1 - pbinom(q = R_i - 1, size = S, prob = z_i)
  2*min(x1, x2)
}