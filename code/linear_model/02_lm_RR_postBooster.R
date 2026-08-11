#updated by Deus & Dan
#06/17/2026

# ============================================================================

# LOG-LOG LINEAR model version of 02_RR_postBooster.R.
#
# For each head-to-head trial (PCV13 vs PCV7, PCV13 vs PCV10, PCV14 vs PCV13,
# PCV15 vs PCV13, PCV20 vs PCV13):
# for each posterior draw,
# sample a log-GMC for each PCV from N(log(GMC), SE_logGMC^2) and predict the
# log-risk of colonization from the LINEAR (no change point) model;
# compute log-RR = log-risk(higher PCV) - log-risk(lower PCV) per draw,
# summarise as median + 95% credible interval + SD;
# pool study-level estimates within a comparison with inverse-variance
# weights (1 / SD^2).

# For the linear model this reduces analytically to
#   log(RR)_s = b1[s] * (log_GMC_high - log_GMC_low) = b1[s] * log(GMR)
# but we keep the same posterior-sample -> per-draw log-RR machinery as the
# change-point code so that trial-level GMC uncertainty is propagated in
# exactly the same way.

# ============================================================================

#set seed for reproducibility
set.seed(20250530)

#serotype map (must match the index used when fitting the Stan model)
pcv7_st  <- c("4", "6B", "9V", "14", "18C", "19F", "23F")
pcv13_extra <- c("1", "3", "5", "6A", "7F", "19A")
serotypes <- c(pcv7_st, pcv13_extra)
df_map <- tibble(serotype = serotypes, st = seq_along(serotypes), in_pcv7 = serotype %in% pcv7_st)

#posterior parameter matrix (rows = draws, cols = b0[s], b1[s], b3[s])
post_mat <- as.matrix(linear_postBooster_fit$draws(variables = c("b0", "b1", "b3"), format = "draws_matrix"))

n_draws <- nrow(post_mat)
J <- length(serotypes)

#read a head-to-head xlsx and convert GMC/CI to log scale
#se on the log scale comes from the width of the (log) CI / (2 * 1.96)
log_transform_igg <- function(path) {
  df <- readxl::read_excel(path) %>%
    dplyr::mutate(across(c(GMC, upper_limit, lower_limit),
                         ~ ifelse(.x == 0, 0.01, .x))) %>%
    dplyr::mutate(log_GMC      = log(GMC),
                  log_upper    = log(upper_limit),
                  log_lower    = log(lower_limit),
                  se_log_GMC   = abs(log_upper - log_lower) / (1.96 * 2)) %>%
    dplyr::inner_join(df_map, by = "serotype")
  df
}

#predicted log-risk for one serotype and one (GMC, SE_log_GMC),
#vectorised across posterior draws.  For each draw the log-GMC is itself
#sampled from N(log_GMC_mean, se_log_GMC^2), propagating trial-level GMC
#uncertainty into the posterior of log-risk.
predict_log_risk <- function(st, log_GMC_mean, se_log_GMC, ethnic = c("jewish", "bedouin")) {
  ethnic <- match.arg(ethnic)
  b0 <- post_mat[, paste0("b0[", st, "]")]
  b1 <- post_mat[, paste0("b1[", st, "]")]
  b3 <- post_mat[, paste0("b3[", st, "]")]

  lg <- rnorm(n_draws, mean = log_GMC_mean, sd = se_log_GMC)

  if (ethnic == "bedouin") b0 + b3 + b1 * lg
  else                     b0 +      b1 * lg
}

#for one head-to-head dataframe, compute per-study log-RR
#(higher-valency vs lower-valency) by serotype
calc_RR <- function(df, pcv_low, pcv_high, ethnic = "jewish") {
  studies   <- unique(df$study_id)

  out <- list()
  for (st in seq_along(serotypes)) {
    sero_name <- serotypes[st]
    sub <- df %>% filter(serotype == sero_name)
    if (nrow(sub) == 0) next
    for (s in unique(sub$study_id)) {
      d  <- sub %>% filter(study_id == s)
      d1 <- d %>% filter(vaccine == pcv_low)
      d2 <- d %>% filter(vaccine == pcv_high)
      if (nrow(d1) != 1 || nrow(d2) != 1) next

      lr1 <- predict_log_risk(st, d1$log_GMC, d1$se_log_GMC, ethnic)
      lr2 <- predict_log_risk(st, d2$log_GMC, d2$se_log_GMC, ethnic)
      log_RR <- lr2 - lr1

      out[[length(out) + 1]] <- tibble(
        st        = st,
        serotype  = sero_name,
        study_id  = s,
        logRR     = median(log_RR),
        logRR_lci = quantile(log_RR, 0.025),
        logRR_uci = quantile(log_RR, 0.975),
        sd_logRR  = sd(log_RR)
      )
    }
  }
  bind_rows(out)
}

#inverse-variance pooled mean across studies, per serotype.
get_wtRR <- function(df_RR, comparison) {
  df_RR %>%
    mutate(inv_var = 1 / (sd_logRR^2)) %>%
    group_by(serotype) %>%
    summarise(
      total_inv_var = sum(inv_var),
      pooled_logRR  = sum(logRR * inv_var) / total_inv_var,
      se_logRR      = sqrt(1 / total_inv_var),
      lci_logRR     = pooled_logRR - 1.96 * se_logRR,
      uci_logRR     = pooled_logRR + 1.96 * se_logRR,
      .groups       = "drop"
    ) %>%
    mutate(Comparison = comparison,
           serotype   = factor(serotype, levels = serotypes)) %>%
    arrange(serotype)
}

#run the five head-to-head comparisons
comparisons <- list(
  "PCV13 vs PCV7"  = list(file = here::here("data", "df_13v7_postboost_n5.xlsx"),  pcv_low = "PCV7",  pcv_high = "PCV13"),
  "PCV13 vs PCV10" = list(file = here::here("data", "df_13v10_postboost_n1.xlsx"), pcv_low = "PCV10", pcv_high = "PCV13"),
  "PCV14 vs PCV13" = list(file = here::here("data", "df_14v13_postboost_n2.xlsx"), pcv_low = "PCV13", pcv_high = "PCV14"),
  "PCV15 vs PCV13" = list(file = here::here("data", "df_15v13_postboost_n11.xlsx"),pcv_low = "PCV13", pcv_high = "PCV15"),
  "PCV20 vs PCV13" = list(file = here::here("data", "df_20v13_postboost_n4.xlsx"), pcv_low = "PCV13", pcv_high = "PCV20")
)

per_study <- list()
pooled    <- list()
for (cmp in names(comparisons)) {
  cfg <- comparisons[[cmp]]
  stopifnot(file.exists(cfg$file))
  df  <- log_transform_igg(cfg$file)
  per_study[[cmp]] <- calc_RR(df, cfg$pcv_low, cfg$pcv_high) %>% mutate(Comparison = cmp)
  pooled[[cmp]]    <- get_wtRR(per_study[[cmp]], cmp)
  cat(sprintf("  %s : %d studies, %d serotype-rows\n", cmp, length(unique(per_study[[cmp]]$study_id)), nrow(per_study[[cmp]])))
}

df_per_study <- bind_rows(per_study)
df_pooled <- bind_rows(pooled)

write_csv(df_per_study, here::here("output", 'postBooster', "lm_rVE_per_study.csv"))
write_csv(df_pooled,    here::here("output", 'postBooster', "lm_rVE_pooled.csv"))

#plot pooled RR for shared 7 PCV7 serotypes and for PCV13 extra 6 serotypes
df_plot <-
  df_pooled %>%
  dplyr::mutate(RR = exp(pooled_logRR),
                lci_RR = exp(lci_logRR),
                uci_RR = exp(uci_logRR),
                Comparison = factor(Comparison, levels = c("PCV13 vs PCV7", "PCV13 vs PCV10", "PCV14 vs PCV13", "PCV15 vs PCV13", "PCV20 vs PCV13")))

base_theme <-
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey95", colour = NA), panel.grid.minor = element_blank(), legend.position  = "top")

p_shared <-
  df_plot %>%
  filter(serotype %in% pcv7_st) %>%
  ggplot(aes(x = serotype, y = RR, ymin = lci_RR, ymax = uci_RR)) +
  geom_pointrange(colour = "steelblue4") +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  facet_wrap(~ Comparison, ncol = 5, scales = 'free_x') +
  labs(x = "Serotype", y = "Relative risk of colonization", title = "Shared PCV7 serotypes (linear model)") +
  base_theme
print(p_shared)

p_extra <-
  df_plot %>%
  filter(serotype %in% pcv13_extra, Comparison != "PCV13 vs PCV7") %>%
  ggplot(aes(x = serotype, y = RR, ymin = lci_RR, ymax = uci_RR)) +
  geom_pointrange(colour = "brown4") +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  facet_wrap(~ Comparison, ncol = 4, scales = 'free_x') +
  labs(x = "Serotype", y = "Relative risk of colonization", title = "Additional PCV13 serotypes (1, 3, 5, 6A, 7F, 19A) - linear model") +
  base_theme
print(p_extra)

p_combined <-
  p_shared / p_extra +
  plot_annotation(
    title    = "Relative risk of colonization, higher- vs lower-valency PCV (linear model)",
    subtitle = "Posterior median and 95% CI, pooled across head-to-head trials with inverse-variance weights",
    caption  = "RR > 1 indicates higher risk with the higher-valency PCV")
print(p_combined)

ggsave(here::here("output", 'postBooster', "lm_relative_VE.png"), p_combined, width = 12, height = 9, dpi = 150)

#per-study forest plots
for (cmp in names(comparisons)) {
  d <-
    df_per_study %>%
    filter(Comparison == cmp) %>%
    mutate(serotype = factor(serotype, levels = serotypes),
           RR  = exp(logRR),
           lci = exp(logRR_lci),
           uci = exp(logRR_uci))

  pooled_d <-
    df_pooled %>% filter(Comparison == cmp) %>%
    mutate(study_id = "Pooled",
           RR  = exp(pooled_logRR),
           lci = exp(lci_logRR),
           uci = exp(uci_logRR))

  pp <-
    ggplot(d, aes(x = study_id, y = RR, ymin = lci, ymax = uci)) +
    geom_pointrange(colour = "grey30") +
    geom_pointrange(data = pooled_d, aes(x = study_id, y = RR, ymin = lci, ymax = uci), shape = 18, size = 0.9, colour = "firebrick") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    facet_wrap(~ serotype, scales = "free_y") +
    scale_y_log10() +
    coord_flip() +
    labs(x = NULL, y = "Relative risk of colonization", title = paste0("Per-study RR (linear): ", cmp)) +
    theme_bw(base_size = 9) +
    theme(strip.background = element_rect(fill = "grey95", colour = NA))

  ggsave(here::here("output", 'postBooster', paste0("lm_per_study_RR_", gsub(" ", "_", cmp), ".png")), pp, width = 11, height = 8, dpi = 150)
}

#pooled estimates
df_pooledNS <- print(df_plot %>% select(Comparison, serotype, RR, lci_RR, uci_RR), n = Inf)
write_csv(df_pooledNS, here::here("output", 'postBooster', "lm_rVE_pooled_exp.csv"))
