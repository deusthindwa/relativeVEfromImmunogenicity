// updated by Deus & Dan
// 06/17/2026

/* ============================================================================

 Hierarchical Bayesian LOG-LOG LINEAR model for the risk of pneumococcal
 colonization as a function of vaccine-induced serum IgG (log scale).

 For each observation i (i = 1..N):

   y_i ~ Bernoulli(pi_i)

   log(q_i) = b0[s_i] + bedouin_i * b3[s_i] + b1[s_i] * log_GMC_i

   pi_i = min(q_i, 1 - eps)  // log-link constraint

 Serotype-level parameters (s = 1..J) are drawn from common hierarchical priors:
   b0[s] ~ N(mu_b0, sigma_b0)   // serotype intercept
   b1[s] ~ N(mu_b1, sigma_b1)   // serotype slope on log_GMC
   b3[s] ~ N(mu_b3, sigma_b3)   // ethnicity-specific intercept shift

  b0 - "Intercept"
  b1 - "Slope"
  b3 - "Ethnicity-specific intercept"

 This is the log-log linear counterpart of the change-point model: there is no change point, so log-risk depends on log_GMC over the entire
 observed range (no hinge). All other structural choices including the hierarchical pooling across serotypes, non-centered parameterization for NUTS, and
 weakly-informative half-Cauchy priors on the population SDs -- mirror the change-point specification.

 ============================================================================
*/

data {
  int<lower=1>                 N;         // total observations
  int<lower=1>                 J;         // number of serotypes
  array[N] int<lower=0,upper=1> y;        // colonization indicator
  array[N] int<lower=1,upper=J> st;       // serotype index
  vector[N]                    log_GMC;   // log IgG concentration
  array[N] int<lower=0,upper=1> bedouin;  // ethnicity flag
}

parameters {
  // Population means
  real mu_b0;
  real mu_b1;
  real mu_b3;

  // Population SDs (half-Cauchy via <lower=0>)
  real<lower=0> sigma_b0;
  real<lower=0> sigma_b1;
  real<lower=0> sigma_b3;

  // Non-centered serotype deviations
  vector[J] z_b0;
  vector[J] z_b1;
  vector[J] z_b3;
}

transformed parameters {
  vector[J] b0 = mu_b0 + sigma_b0 * z_b0;
  vector[J] b1 = mu_b1 + sigma_b1 * z_b1;
  vector[J] b3 = mu_b3 + sigma_b3 * z_b3;
}

model {
  // Hyperpriors
  mu_b0 ~ normal(0, 100);
  mu_b1 ~ normal(0, 100);
  mu_b3 ~ normal(0, 100);

  sigma_b0 ~ cauchy(0, 2.5);
  sigma_b1 ~ cauchy(0, 2.5);
  sigma_b3 ~ cauchy(0, 2.5);

  // Standard-normal deviations -> non-centered hierarchical priors
  z_b0 ~ std_normal();
  z_b1 ~ std_normal();
  z_b3 ~ std_normal();

  // Likelihood
  for (i in 1:N) {
    real eta = b0[st[i]] + bedouin[i] * b3[st[i]] + b1[st[i]] * log_GMC[i];
    real q   = exp(eta);
    real p   = fmin(q, 1 - 1e-6);
    target += bernoulli_lpmf(y[i] | p);
  }
}

generated quantities {
  // Log-likelihood for LOO / WAIC
  vector[N] log_lik;
  for (i in 1:N) {
    real eta = b0[st[i]] + bedouin[i] * b3[st[i]] + b1[st[i]] * log_GMC[i];
    real q   = exp(eta);
    real p   = fmin(q, 1 - 1e-6);
    log_lik[i] = bernoulli_lpmf(y[i] | p);
  }
}
