# Setup ----
library(bdlvm)
library(blavaan)
library(brms)
library(dplyr)
library(future)
library(ggplot2)
library(ggdist)
library(here)
library(latex2exp)
library(purrr)
library(SBC)
library(tarchetypes)
library(targets)
library(tidyr)
library(patchwork)
library(qs)
plan(multisession)
options(
  SBC.min_chunk_size = 5,
  brms.backend = "cmdstanr"
)
tar_option_set(
  format = "qs", memory = "transient",
  garbage_collection = TRUE
)
tar_source("R")
# Simulation code ----
tar_plan(
  sim_df = simulation_conditions() %>%
    filter(id == "gdsem", num_precon != 5),
  all_generators = sim_df %>%
    select(
      id, batch,
      formulas, priors,
      ifs_priors, num_precon,
      data_template, num_obs, num_rep
    ),
  tar_target(
    each_generator,
    all_generators,
    pattern = map(all_generators)
  ),
  tar_target(
    datasets,
    utils.make_data(each_generator),
    pattern = map(each_generator)
  ),
  tar_target(
    all_models,
    sim_df %>%
      left_join(datasets) %>%
      select(
        id, batch,
        formulas, priors,
        data, data_template,
        num_obs, ifs_priors, num_precon
      )
  ),
  tar_target(
    each_model,
    all_models,
    pattern = map(all_models)
  ),
  tar_target(
    full_results,
    utils.sbc_compute(
      each_model,
      keep_fits = TRUE
    ),
    pattern = map(each_model),
    iteration = "list"
  ),
  tar_target(
    each_statstable,
    utils.extract_stats(each_model, full_results),
    pattern = map(each_model, full_results)
  ),
  tar_target(
    each_backtable,
    utils.extract_backend(each_model, full_results),
    pattern = map(each_model, full_results)
  ),
  full_statstable = each_statstable %>% mutate(
      sim_id = as.numeric(as.factor(paste0(batch, sim_id))),
      .by = c(id, num_obs, num_precon)
    ) %>% select(-batch) %>% arrange(num_precon, sim_id),
  full_backtable = each_backtable %>% mutate(
      sim_id = as.numeric(as.factor(paste0(batch, sim_id))),
      .by = c(id, num_obs, num_precon)
    ) %>% select(-batch) %>% arrange(num_precon, sim_id),
  # what now? plots!
  plots = {
    p1 = full_statstable %>%
      filter(
        num_precon %in% c(0, 6, 10),
        param_type %in%
          c("Factor slopes", "Item std. dev.")) %>%
      utils.sbc_ecdf()
    p2 = utils.rhatplot(full_backtable, "plasma")
    p3 = utils.scholzplot(full_statstable)
    (p1 / (p3 | p2)) + plot_layout(heights = c(0.6,0.4))
  },
  tar_target(
    exported_plot, {
      path = here("plots", "sem_plot.pdf")
      ggsave(path, plot = plots,
        width = 8, height = 6, device = cairo_pdf)
    },
    format = "file"
  )
)
