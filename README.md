# Change-Point and Linear Models of PCV effectiveness against targeted serotypes 

### (Update based on original work by Wong et al., JID 2025)

Stan + NUTS implementation of the hierarchical Bayesian change-point model (& equiv Linear model)

## model

For child–visit–serotype observation *i* with serotype index *s(i)*:

```
y_i ~ Bernoulli(p_i)
log(q_i) = b0[s] + bedouin_i · b3[s] + b1[s] · max(log_GMC_i − cp[s], 0)
log(q_i) = b0[s] + bedouin_i · b3[s] + b1[s] · log_GMC_i
p_i      = min(q_i, 1 − ε)
```

Serotype-level vectors `b0`, `b1`, `b3`, `cp` are pooled across the 13
PCV13 serotypes via Normal hierarchical priors. The same informative prior on `mu_cp` 
is used e.g., `mu_cp ~ N(mean(log_GMC), 1 / (2·var(log_GMC)))` for the change-point mean
- but may toggle with `flat_cp_prior = 1` to switch to a diffuse variant
- precisions are half-Cauchy SDs, which are well-behaved under HMC.
- every hierarchical parameter uses a **non-centered parameterization** facilitating smooth NUTS geometry around small `sigma_*`.


## running

```r
# requires primary packages including: cmdstanr (+ cmdstan), posterior, bayesplot, readxl, patchwork, tidyverse
# stage 1: fit the change-point or linear model with NUTS 
# stage 2: generate colonisation risk curves wrt change-point or linear model
# stage 3: compute the relative risk (RR) following vaccination with higher and lower-valency PCV (based on using fitted posterior & head-to-head GMC summaries)
# stage 4: compute the absolute or relative VE against carriage for PCV13/PCV15/PCV20/PCV14/PCV10 versus lower-valency comparator
```

#stages 1-2: sampler settings
- 4 chains × 1000 warmup + 1000 sampling iterations
- `adapt_delta = 0.95`, `max_treedepth = 12` (hinge is non-smooth for change point model, higher target acceptance keeps divergences down)
- Diagnostics reported: R-hat, bulk/tail ESS, divergent transitions, treedepth saturation.

#stage 3
#for each head-to-head trial (PCV13 vs PCV7, PCV15 vs PCV13, PCV20 vs PCV13, PCV14 vs PCV13, PCV25 vs PCV13):
#for each posterior draw, 
#sample a log-GMC for each PCV from N(log(GMC), SE_logGMC^2) and predict the log-risk of colonization from change-point or linear model
#compute log-RR = log-risk(higher PCV) - log-risk(lower PCV) per draw,
#summarise as median + 95% credible interval + SD;
#pool study-level estimates within a comparison with inverse-variance
#weights (1 / SD^2).

#stage 4
# convert the pooled relative-risk estimates produced in stage 3
# absolute vaccine efficacy (VE) against colonization for PCV10, PCV13, PCV14, PCV15, PCV20
# anchored at an assumed VE(PCV7) = 60% and absolute VE(PCV13) from the data.
#
# chain rule:
#
# for serotypes in PCV7 (4, 6B, 9V, 14, 18C, 19F, 23F): VE0(PCV7) = 0.60
#   VE_PCV13[s] = 1 - (1 - VE0) * RR_{13v7}[s]
#   VE_PCV15[s] = 1 - (1 - VE0) * RR_{13v7}[s] * RR_{15v13}[s]
#   VE_PCV20[s] = 1 - (1 - VE0) * RR_{13v7}[s] * RR_{20v13}[s]
#
# for PCV13-only serotypes (1, 3, 5, 6A, 7F, 19A): PCV7 confers no
#   antibody-mediated protection because the serotype is not in PCV7, so VE0(PCV7) = 0 and the change-point model's RR_{13v7}[s] (computed
#   from the head-to-head data, where PCV7-arm IgG is at baseline for these serotypes) gives the absolute VE for PCV13 directly:
#     VE_PCV13[s] = 1 - RR_{13v7}[s]
#     VE_PCV15[s] = 1 - RR_{13v7}[s] * RR_{15v13}[s]
#     VE_PCV20[s] = 1 - RR_{13v7}[s] * RR_{20v13}[s]
#
# Uncertainty is propagated by Monte-Carlo: for each serotype × comparison we draw log-RR ~ N(pooled_logRR, se_logRR) 
# and chain through, reported intervals are the empirical 2.5 / 97.5 quantiles
