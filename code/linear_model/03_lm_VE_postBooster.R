#updated by Deus & Dan
#06/17/2026

# ============================================================================
# Convert the LINEAR-model pooled relative-risk estimates into absolute VE
# for PCV10, PCV13, PCV14, PCV15, PCV20 anchored at assumed VE(PCV7) = 60%
# for PCV7 serotypes and 0% for the six PCV13-only serotypes.
# ============================================================================

#set seed for reproducibility
set.seed(20250530)

#set VE0
cli <- commandArgs(trailingOnly = TRUE)
ve0_arg <- sub("^--ve0=", "", grep("^--ve0=", cli, value = TRUE))
VE0 <- if (length(ve0_arg) == 1) as.numeric(ve0_arg) else 0.60
stopifnot(VE0 > 0, VE0 < 1)

#read pooled log-RR
pcv7_st     <- c("4", "6B", "9V", "14", "18C", "19F", "23F")
pcv13_extra <- c("1", "3", "5", "6A", "7F", "19A")
serotypes   <- c(pcv7_st, pcv13_extra)

rr <- read_csv(here::here('output', 'postBooster', 'lm_RR_pooled.csv'), show_col_types = FALSE) %>%
  mutate(serotype = factor(serotype, levels = serotypes),
         Comparison = factor(Comparison,
                             levels = c("PCV13 vs PCV7",
                                        "PCV13 vs PCV10",
                                        "PCV14 vs PCV13",
                                        "PCV15 vs PCV13",
                                        "PCV20 vs PCV13")))

#wide: one row per serotype, one column per comparison
rr_wide <- rr %>%
  dplyr::select(serotype, Comparison, pooled_logRR, se_logRR) %>%
  pivot_wider(names_from  = Comparison,
              values_from = c(pooled_logRR, se_logRR))

#monte-Carlo chaining
M <- 20000L

draw_logRR <- function(mu, se, M) {
  if (is.na(mu) || is.na(se)) rep(NA_real_, M) else rnorm(M, mu, se)
}

abs_ve <- lapply(seq_len(nrow(rr_wide)), function(i) {
  row <- rr_wide[i, ]
  s   <- as.character(row$serotype)
  in_pcv7 <- s %in% pcv7_st

  #sample log-RRs (independent across comparisons)
  rr_13v7  <- exp(draw_logRR(row$`pooled_logRR_PCV13 vs PCV7`,  row$`se_logRR_PCV13 vs PCV7`,  M))
  rr_13v10 <- exp(draw_logRR(row$`pooled_logRR_PCV13 vs PCV10`, row$`se_logRR_PCV13 vs PCV10`, M))
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
    rr_13v7 = 1 #delete this to compute absolute VE of pcv13-only serotypes
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
  bind_rows(summ(ve13, "PCV13"), summ(ve10, "PCV10"), summ(ve14, "PCV14"), summ(ve15, "PCV15"), summ(ve20, "PCV20"))
}) %>% bind_rows()


abs_ve <- abs_ve %>%
  dplyr::mutate(serotype = factor(serotype, levels = serotypes),
         vaccine  = factor(vaccine,  levels = c("PCV10", "PCV13", "PCV14", "PCV15", "PCV20")),
         in_pcv7  = serotype %in% pcv7_st) %>%
  dplyr::arrange(serotype, vaccine)

write_csv(abs_ve, here::here('output', 'postBooster', sprintf("lm_abs_rel_VE.csv", round(VE0 * 100))))

#absolute VE table
print(abs_ve %>% transmute(serotype,
                           vaccine,
                           VE_pct = round(100 * median, 1),
                           CI_95 = sprintf("(%.1f, %.1f)", 100 * lci, 100 * uci)), n=Inf)

#plot the absolute VE
p_absolute <-
  abs_ve %>%
  dplyr::filter(in_pcv7 == TRUE) %>%
  dplyr::mutate(pcvind = factor(if_else(in_pcv7 == TRUE, 'PCV7 serotypes', 'PCV13-only serotypes'), levels = c('PCV7 serotypes', 'PCV13-only serotypes'))) %>%
  ggplot(aes(x = serotype, y = median, ymin = lci, ymax = uci, colour = vaccine)) +
  geom_hline(yintercept = VE0, linetype = "dashed", colour = "grey50") +
  facet_grid(. ~ pcvind, scales = "free_x", space = "free_x") +
  geom_pointrange(position = position_dodge(width = 0.55), size = 0.4) +
  scale_y_continuous(limits = c(0.35, 0.75), labels = scales::percent, n.breaks = 6) +
  scale_colour_manual(values = c("PCV10" = "green4", "PCV13" = "steelblue4", "PCV14" = "gray40", "PCV15" = "darkorange2", "PCV20" = "firebrick")) +
  labs(x = "Serotype", y = "Absolute vaccine efficacy \nagainst colonization (%)", colour = "Vaccine",
       title   = "VE against colonization (linear model), chained from RR estimates (linear model, post-booster)",
       subtitle = sprintf("Assumed VE(PCV7) = %.0f%% for PCV7 serotypes and 0%% for the six PCV13-only serotypes", 100 * VE0)) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position  = "top")
print(p_absolute)

p_relative <-
  abs_ve %>%
  dplyr::filter(in_pcv7 == FALSE) %>%
  dplyr::mutate(pcvind = factor(if_else(in_pcv7 == TRUE, 'PCV7 serotypes', 'PCV13-only serotypes'), levels = c('PCV7 serotypes', 'PCV13-only serotypes'))) %>%
  ggplot(aes(x = serotype, y = median, ymin = lci, ymax = uci, colour = vaccine)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  facet_grid(. ~ pcvind, scales = "free_x", space = "free_x") +
  geom_pointrange(position = position_dodge(width = 0.55), size = 0.4) +
  scale_y_continuous(limits = c(-0.3, 0.3), labels = scales::percent, n.breaks = 6) +
  scale_colour_manual(values = c("PCV10" = "green4", "PCV13" = "steelblue4", "PCV14" = "gray40", "PCV15" = "darkorange2", "PCV20" = "firebrick")) +
  labs(x = "Serotype", y = "Vaccine efficacy against colonization \nrelative to PCV13 VE (%)", colour = "Vaccine", title = "") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position  = "none")

print(p_relative)

print(p_absolute/p_relative)
ggsave(here::here('output', 'postBooster', sprintf("lm_abs_rel_VE.png", round(VE0 * 100))), p_absolute, width = 12, height = 8, dpi = 150)
