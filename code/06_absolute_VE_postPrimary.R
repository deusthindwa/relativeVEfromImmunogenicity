#updated by Deus & Dan 
#06/17/2026

# ============================================================================
#
# convert the pooled relative-risk estimates produced by relative_VE.R into
# absolute vaccine efficacy (VE) against colonization for PCV10, PCV13, PCV14, PCV15, and
# PCV20, anchored at an assumed VE(PCV7) = 60% and absolute VE(PCV13) from the data.
#
# chain rule:
#
#   For serotypes in PCV7 (4, 6B, 9V, 14, 18C, 19F, 23F): VE0(PCV7) = 0.60
#     VE_PCV13[s] = 1 - (1 - VE0) * RR_{13v7}[s]
#     VE_PCV15[s] = 1 - (1 - VE0) * RR_{13v7}[s] * RR_{15v13}[s]
#     VE_PCV20[s] = 1 - (1 - VE0) * RR_{13v7}[s] * RR_{20v13}[s]
#
#   For PCV13-only serotypes (1, 3, 5, 6A, 7F, 19A): PCV7 confers no
#   antibody-mediated protection because the serotype is not in PCV7,
#   so VE0(PCV7) = 0 and the change-point model's RR_{13v7}[s] (computed
#   from the head-to-head data, where PCV7-arm IgG is at baseline for
#   these serotypes) gives the absolute VE for PCV13 directly:
#     VE_PCV13[s] = 1 - RR_{13v7}[s]
#     VE_PCV15[s] = 1 - RR_{13v7}[s] * RR_{15v13}[s]
#     VE_PCV20[s] = 1 - RR_{13v7}[s] * RR_{20v13}[s]
#
# Uncertainty is propagated by Monte-Carlo: for each serotype × comparison we
# draw log-RR ~ N(pooled_logRR, se_logRR) and chain through; reported
# intervals are the empirical 2.5 / 97.5 quantiles.

# ============================================================================

#set seed for reproducibility
set.seed(20260605L)

#set VE0
cli <- commandArgs(trailingOnly = TRUE)
ve0_arg <- sub("^--ve0=", "", grep("^--ve0=", cli, value = TRUE))
VE0 <- if (length(ve0_arg) == 1) as.numeric(ve0_arg) else 0.60
stopifnot(VE0 > 0, VE0 < 1)

#read pooled log-RR
pcv7_st <- c("4", "6B", "9V", "14", "18C", "19F", "23F")
pcv13_extra <- c("1", "3", "5", "6A", "7F", "19A")
serotypes <- c(pcv7_st, pcv13_extra)

rr <- read_csv(here::here('output', 'postPrimary','rVE_pooled.csv'), show_col_types = FALSE) %>%
  mutate(serotype = factor(serotype, levels = serotypes),
         Comparison = factor(Comparison,
                             levels = c("PCV13 vs PCV7",
                                        "PCV13 vs PCV10",
                                        "PCV14 vs PCV13",
                                        "PCV15 vs PCV13",
                                        "PCV20 vs PCV13")))

#wide: one row per serotype, one column per comparison
rr_wide <- rr %>%
  select(serotype, Comparison, pooled_logRR, se_logRR) %>%
  pivot_wider(names_from  = Comparison,
              values_from = c(pooled_logRR, se_logRR))

#monte-Carlo chaining
M <- 20000L     #monte-Carlo draws per serotype

draw_logRR <- function(mu, se, M) {
  if (is.na(mu) || is.na(se)) rep(NA_real_, M) else rnorm(M, mu, se)
}

abs_ve <- lapply(seq_len(nrow(rr_wide)), function(i) {
  row <- rr_wide[i, ]
  s   <- as.character(row$serotype)
  in_pcv7 <- s %in% pcv7_st
  
  #sample log-RRs (independent across comparisons -- treating the
  #head-to-head trials as independent meta-analytic inputs)
  rr_13v7  <- exp(draw_logRR(row$`pooled_logRR_PCV13 vs PCV7`, row$`se_logRR_PCV13 vs PCV7`,  M))
  rr_13v10 <- exp(draw_logRR(row$`pooled_logRR_PCV13 vs PCV10`, row$`se_logRR_PCV13 vs PCV10`,  M))
  rr_14v13 <- exp(draw_logRR(row$`pooled_logRR_PCV14 vs PCV13`, row$`se_logRR_PCV14 vs PCV13`, M))
  rr_15v13 <- exp(draw_logRR(row$`pooled_logRR_PCV15 vs PCV13`, row$`se_logRR_PCV15 vs PCV13`, M))
  rr_20v13 <- exp(draw_logRR(row$`pooled_logRR_PCV20 vs PCV13`, row$`se_logRR_PCV20 vs PCV13`, M))
  
  if (in_pcv7) {
    ve13 <- 1 - (1 - VE0) * rr_13v7
    ve10 <- 1 - (1 - VE0) * rr_13v7 * rr_13v10
    ve14 <- 1 - (1 - VE0) * rr_13v7 * rr_14v13
    ve15 <- 1 - (1 - VE0) * rr_13v7 * rr_15v13
    ve20 <- 1 - (1 - VE0) * rr_13v7 * rr_20v13
  } else {
    # PCV13-only serotype: PCV7 provides no antibody-mediated protection,
    # so VE(PCV7) = 0 and the model-estimated RR_{13v7} already reflects
    # the absolute risk reduction conferred by PCV13.
    ve13 <- 1-rr_13v7
    ve10 <- 1-rr_13v7 * rr_13v10
    ve14 <- 1-rr_13v7 * rr_14v13
    ve15 <- 1-rr_13v7 * rr_15v13
    ve20 <- 1-rr_13v7 * rr_20v13
  }
  
  summ <- function(v, lbl) {
    if (all(is.na(v))) return(NULL)
    tibble(serotype = s, vaccine = lbl,
           median = median(v, na.rm = TRUE),
           lci    = quantile(v, 0.025, na.rm = TRUE, names = FALSE),
           uci    = quantile(v, 0.975, na.rm = TRUE, names = FALSE))
  }
  bind_rows(summ(ve10, "PCV10"), summ(ve13, "PCV13"), summ(ve14, "PCV14"), summ(ve15, "PCV15"), summ(ve20, "PCV20"))
}) %>% bind_rows()


abs_ve <- abs_ve %>%
  mutate(serotype = factor(serotype, levels = serotypes),
         vaccine  = factor(vaccine,  levels = c("PCV10", "PCV13", "PCV14", "PCV15", "PCV20")),
         in_pcv7  = serotype %in% pcv7_st) %>%
  arrange(serotype, vaccine)

write_csv(abs_ve, here::here('output', 'postPrimary', sprintf("absolute_VE_ve0_%02d.csv", round(VE0 * 100))))

#absolute VE table
print(abs_ve %>% transmute(serotype, 
                           vaccine, 
                           VE_pct = round(100 * median, 1), 
                           CI_95 = sprintf("(%.1f, %.1f)", 100 * lci, 100 * uci)), n=Inf)


#plot the absolute VE
p_absolute <- 
  abs_ve %>%
  dplyr::mutate(pcvind = factor(if_else(in_pcv7 == TRUE, 'PCV7 serotypes', 'PCV13-only serotypes'), levels = c('PCV7 serotypes', 'PCV13-only serotypes'))) %>%
  ggplot(aes(x = serotype, y = 100 * median, ymin = 100 * lci, ymax = 100 * uci, colour = vaccine)) +
  geom_hline(yintercept = 100 * VE0, linetype = "dashed", colour = "grey50") +
  facet_grid(. ~ pcvind, scales = "free_x", space = "free_x") +
  geom_pointrange(position = position_dodge(width = 0.55), size = 0.4) +
  scale_colour_manual(values = c("PCV10" = "green4", "PCV13" = "steelblue4", "PCV14" = "gray40", "PCV15" = "darkorange2", "PCV20" = "firebrick")) +
  labs(x = "Serotype", y = "Absolute vaccine efficacy \nagainst colonization (%)", colour = "Vaccine", title   = "Absolute VE against colonization, chained from RR estimates", subtitle = sprintf("Assumed VE(PCV7) = %.0f%% for PCV7 serotypes; VE(PCV13) = 0%% for the six PCV13-only serotypes", 100 * VE0, 100 * VE0)) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position  = "top")

print(p_absolute)
ggsave(here::here('output', 'postPrimary', sprintf("absolute_VE_ve0_%02d.png", round(VE0 * 100))), p_absolute, width = 11, height = 6, dpi = 150)
