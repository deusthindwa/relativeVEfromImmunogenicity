# A Change-Point Logistic Model and Linear Model of PCV effetiveness against targeted serotypes 

### (Worked based on original work by Wong et al., JID 2025)

Stan + NUTS implementation of the hierarchical Bayesian Logistic change-point model (& equiv Linear model)

## model

For child–visit–serotype observation *i* with serotype index *s(i)*:

```
y_i ~ Bernoulli(p_i)
log(q_i) = b0[s] + bedouin_i · b3[s] + b1[s] · max(log_GMC_i − cp[s], 0)
p_i      = min(q_i, 1 − ε)
```

Serotype-level vectors `b0`, `b1`, `b3`, `cp` are pooled across the 13
PCV13 serotypes via Normal hierarchical priors. The same informative prior on `mu_cp` 
is used e.g., `mu_cp ~ N(mean(log_GMC), 1 / (2·var(log_GMC)))` for the change-point mean
- but may toggle with `flat_cp_prior = 1` to switch to a diffuse variant
- precisions are half-Cauchy SDs, which are well-behaved under HMC.
- every hierarchical parameter uses a **non-centered parameterization**
- so NUTS does not fight the funnel geometry that a centered version would create at small `sigma_*`.

## running

```r
# Requires: cmdstanr (+ cmdstan), posterior, bayesplot, readxl, patchwork, tidyverse

# fit the change-point model with NUTS - 01_run_changepoint.R
# compute relative VE against carriage for PCV13/PCV15/PCV20/PCV14/PCV25 versus lower-valency comparator, using fitted posterior & head-to-head IgG summaries
```

## sampler settings

- 4 chains × 1000 warmup + 1000 sampling iterations
- `adapt_delta = 0.95`, `max_treedepth = 12` (hinge is non-smooth at the change point; the higher target acceptance keeps divergences down)
- Diagnostics reported: R-hat, bulk/tail ESS, divergent transitions, treedepth saturation.
