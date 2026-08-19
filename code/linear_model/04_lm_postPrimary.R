#updated by Deus & Dan
#06/17/2026

# ============================================================================

#set seed for reproducibility
set.seed(20250530)

serotypes <- c("4", "6B", "9V", "14", "18C", "19F", "23F", "1", "3",  "5",  "6A", "7F", "19A")

#locate the project data
dat_postPrimary <-
  rio::import(here::here('data', 'df_colonization_7m.7_24.csv')) %>%
  dplyr::mutate(kid_id = Study_key,
                visit = Age_in_months,
                Serotype = factor(Serotype),
                st = if_else(Serotype=='4', 1, if_else(Serotype=='6B', 2, if_else(Serotype=='9V', 3,
                     if_else(Serotype=='14', 4, if_else(Serotype=='18C', 5, if_else(Serotype=='19F', 6,
                     if_else(Serotype=='23F', 7, if_else(Serotype=='1', 8, if_else(Serotype=='3', 9,
                     if_else(Serotype=='5', 10, if_else(Serotype=='6A', 11, if_else(Serotype=='7F', 12, if_else(Serotype=='19A', 13, NA_integer_))))))))))))),
                bedouin = if_else(Ethnicity=='Bedouin', 1L, 0L),
                log_GMC = log_GMC,
                Colonization = Colonization) %>%
  dplyr::filter(visit == 7 | visit == 12) %>%
  dplyr::select(kid_id, visit, Serotype, st, bedouin, log_GMC, Colonization)

#summary stats about the data
print(dat_postPrimary %>%
        group_by(Serotype) %>%
        summarise(N = n(), n_col = sum(Colonization), .groups = "drop"))

#compile and fit the Stan model with NUTS
stan_data <- list(
  N       = nrow(dat_postPrimary),
  J       = length(serotypes),
  y       = dat_postPrimary$Colonization,
  st      = dat_postPrimary$st,
  log_GMC = dat_postPrimary$log_GMC,
  bedouin = dat_postPrimary$bedouin
)

#compile stan code
#modigg <- cmdstanr::cmdstan_model(here::here('code', 'linear_model', '00_linear_model.stan'))

# #fit the stan model (this take considerable 30 minutes to runs)
# fit <- modigg$sample(
#   data            = stan_data,
#   chains          = 4,
#   parallel_chains = 4,
#   iter_warmup     = 500,
#   iter_sampling   = 1000,
#   seed            = 20250530,
#   refresh         = 100,
#   adapt_delta     = 0.95,
#   max_treedepth   = 12
# )
# 
# #save the fitted model file
# fit$save_object(here::here('results', "linear_postPrimary_fit.rds"))

#read the fitted model
linear_postPrimary_fit <- rio::import(here::here('results', 'linear_postPrimary_fit.rds'))

#diagnostics
diag_pars <- c("mu_b0", "mu_b1", "mu_b3", "sigma_b0", "sigma_b1", "sigma_b3")
sero_pars <- c(paste0("b0[",  1:13, "]"), paste0("b1[",  1:13, "]"), paste0("b3[",  1:13, "]"))

#print summary estimates
draws <- linear_postPrimary_fit$draws(variables = c(diag_pars, sero_pars), format = "draws_array")
summ  <- posterior::summarise_draws(draws, mean, median, sd, ~ posterior::quantile2(.x, probs = c(0.025, 0.975)), rhat, ess_bulk, ess_tail)
print(summ, n = Inf)

#HMC-specific diagnostics
linear_postPrimary_fit$diagnostic_summary()

#trace plots
bayesplot::color_scheme_set("brightblue")
trace_plot <- mcmc_trace(draws, pars = diag_pars, facet_args = list(ncol = 1))
trace_plot
ggsave(here::here('output', 'postPrimary', "lm_trace_hypers.png"), trace_plot, width = 10, height = 10, dpi = 150)

dens_plot <- mcmc_dens_overlay(draws, pars = diag_pars)
dens_plot
ggsave(here::here('output', 'postPrimary', "lm_dens_hypers.png"), dens_plot, width = 10, height = 6, dpi = 150)

#density plots
for (par in c("b0", "b1", "b3")) {
  pars_par <- paste0(par, "[", 1:13, "]")

  ggsave(here::here('output', 'postPrimary', paste0("lm_trace_", par, ".png")),
         mcmc_trace(draws, pars = pars_par,
                    facet_args = list(ncol = 4)),
         width = 12, height = 8, dpi = 150)

  ggsave(here::here('output', 'postPrimary', paste0("lm_dens_", par, ".png")),
         mcmc_dens_overlay(draws, pars = pars_par,
                           facet_args = list(ncol = 4)),
         width = 12, height = 8, dpi = 150)
}

#pairs plot for a couple of serotypes to inspect funnel or correlations
pairs_plot <- mcmc_pairs(draws, pars = c("b0[1]", "b1[1]", "b3[1]"))
print(pairs_plot)
ggsave(here::here('output', 'postPrimary', "lm_pairs_st1.png"), pairs_plot, width = 8, height = 8, dpi = 150)

#posterior predictive risk curves
post_mat <- as.matrix(linear_postPrimary_fit$draws(variables = c("b0", "b1", "b3"), format = "draws_matrix"))
grid <- expand_grid(
  st = seq_len(13),
  log_GMC = seq(min(dat_postPrimary$log_GMC), max(dat_postPrimary$log_GMC), length.out = 80),
  bedouin = c(0L, 1L)) %>%
  dplyr::mutate(serotype = factor(serotypes[st], levels = serotypes),
                ethnicity = ifelse(bedouin == 1, "Bedouin", "Jewish"))

b0_cols <- paste0("b0[", grid$st, "]")
b1_cols <- paste0("b1[", grid$st, "]")
b3_cols <- paste0("b3[", grid$st, "]")
stopifnot(all(c(b0_cols, b1_cols, b3_cols) %in% colnames(post_mat)))

risk_for_draw <- function(k) {
  b0 <- post_mat[k, b0_cols]
  b1 <- post_mat[k, b1_cols]
  b3 <- post_mat[k, b3_cols]
  pmin(exp(b0 + grid$bedouin * b3 + b1 * grid$log_GMC), 1 - 1e-6)
}

idx <- sample.int(nrow(post_mat), min(1000, nrow(post_mat)))
risk_mat <- vapply(idx, risk_for_draw, numeric(nrow(grid)))
stopifnot(!anyNA(risk_mat))

grid <-
  grid %>%
  dplyr::mutate(median = apply(risk_mat, 1, median),
                lower  = apply(risk_mat, 1, quantile, 0.025),
                upper  = apply(risk_mat, 1, quantile, 0.975))

p_risk <-
  ggplot(grid, aes(x = log_GMC, y = median, colour = ethnicity, fill = ethnicity)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ serotype, scales = "free_y", ncol = 5) +
  scale_y_log10() +
  labs(x = "log(IgG concentration)", y = "Predicted risk of colonization",
       colour = "Ethnicity", fill = "Ethnicity",
       title = "Posterior risk of colonization by serotype (linear model, post-primary)",
       subtitle = "Median (line) and 95% credible interval (ribbon)") +
  theme_minimal(base_size = 14) +
  theme(strip.background = element_rect(fill = "grey95", colour = NA))

print(p_risk)
ggsave(here::here('output', 'postPrimary', "lm_risk_curve.png"), p_risk, width = 12, height = 8, dpi = 150)

#serotype-specific slope summary
slope_summary <- summ %>%
  dplyr::filter(grepl("^b1\\[", variable)) %>%
  dplyr::mutate(st = as.integer(gsub("b1\\[|\\]", "", variable)),
                serotype = factor(serotypes[st], levels = serotypes))

p_slope <-
  ggplot(slope_summary, aes(x = serotype, y = median, ymin = q2.5, ymax = q97.5)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
  labs(x = "Serotype", y = "log-log slope b1",
       title = "Posterior serotype-specific slope on log IgG (post-primary)") +
  theme_minimal(base_size = 11)

print(p_slope)
ggsave(here::here('output', 'postPrimary', "lm_slope.png"), p_slope, width = 8, height = 5, dpi = 150)
