// updated by Deus & Dan 
// 06/17/2026

/* ============================================================================

 Hierarchical Bayesian change-point model for the risk of pneumococcal
 colonization as a function of vaccine-induced serum IgG (log scale)

 For each observation i (i = 1..N):

   y_i ~ Bernoulli(pi_i)
   
   log(q_i) = b0[s_i] + bedouin_i * b3[s_i] + b1[s_i] * max(log_GMC_i - cp[s_i], 0)
   
   pi_i = min(q_i, 1 - eps) // log-link constraint

 Serotype-level parameters (s = 1..J) are drawn from common hierarchical priors:
   b0[s] ~ N(mu_b0, sigma_b0)
   b1[s] ~ N(mu_b1, sigma_b1)
   b3[s] ~ N(mu_b3, sigma_b3)
   cp[s] ~ N(mu_cp, sigma_cp)
   
  b0 - "Intercept"
  b1 -"Slope"
  b3 - "Ethnicity-specific intercept"
  cp - "Change point"

 To avoid the funnel pathology that wrecks NUTS sampling of
 centered hierarchical models, every group-level vector uses a
 non-centered (Matt-trick) parameterization. Cauchy(0, 2.5) priors
 on the SDs are weakly-informative replacements for the
 InverseGamma(0.01, 0.01) precisions used in the original JAGS code.
 
 ============================================================================
*/

data {
  int<lower=1>                 N;            // total observations
  int<lower=1>                 J;            // number of serotypes
  array[N] int<lower=0,upper=1> y;           // colonization indicator
  array[N] int<lower=1,upper=J> st;          // serotype index
  vector[N]                    log_GMC;      // log IgG concentration
  array[N] int<lower=0,upper=1> bedouin;     // ethnicity flag
  real                         log_gmc_mean; // prior mean for mu_cp
  real<lower=0>                log_gmc_sd;   // prior SD   for mu_cp
  int<lower=0,upper=1>         flat_cp_prior;// 1 = use diffuse N(0,100)
}

parameters {
  // Population means
  real mu_b0;
  real mu_b1;
  real mu_b3;
  real mu_cp;

  // Population SDs (half-Cauchy via <lower=0>)
  real<lower=0> sigma_b0;
  real<lower=0> sigma_b1;
  real<lower=0> sigma_b3;
  real<lower=0> sigma_cp;

  // Non-centered serotype deviations
  vector[J] z_b0;
  vector[J] z_b1;
  vector[J] z_b3;
  vector[J] z_cp;
}

transformed parameters {
  vector[J] b0 = mu_b0 + sigma_b0 * z_b0;
  vector[J] b1 = mu_b1 + sigma_b1 * z_b1;
  vector[J] b3 = mu_b3 + sigma_b3 * z_b3;
  vector[J] cp = mu_cp + sigma_cp * z_cp;
}

model {
  // Hyperpriors
  mu_b0 ~ normal(0, 100);
  mu_b1 ~ normal(0, 100);
  mu_b3 ~ normal(0, 100);

  if (flat_cp_prior == 1)
    mu_cp ~ normal(0, 100);                  // diffuse prior
  else
    mu_cp ~ normal(log_gmc_mean, log_gmc_sd); // informative prior

  sigma_b0 ~ cauchy(0, 2.5);
  sigma_b1 ~ cauchy(0, 2.5);
  sigma_b3 ~ cauchy(0, 2.5);
  sigma_cp ~ cauchy(0, 2.5);

  // Standard-normal deviations -> non-centered hierarchical priors
  z_b0 ~ std_normal();
  z_b1 ~ std_normal();
  z_b3 ~ std_normal();
  z_cp ~ std_normal();

  // Likelihood
  for (i in 1:N) {
    real hinge = fmax(log_GMC[i] - cp[st[i]], 0);
    real eta   = b0[st[i]] + bedouin[i] * b3[st[i]] + b1[st[i]] * hinge;
    real q     = exp(eta);
    real p     = fmin(q, 1 - 1e-6);
    target += bernoulli_lpmf(y[i] | p);
  }
}

generated quantities {
  // Log-likelihood for LOO / WAIC
  vector[N] log_lik;
  for (i in 1:N) {
    real hinge = fmax(log_GMC[i] - cp[st[i]], 0);
    real eta   = b0[st[i]] + bedouin[i] * b3[st[i]] + b1[st[i]] * hinge;
    real q     = exp(eta);
    real p     = fmin(q, 1 - 1e-6);
    log_lik[i] = bernoulli_lpmf(y[i] | p);
  }
}
