// generated with brms 2.20.4
functions {
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  array[N] int<lower=-1,upper=2> cens;  // indicates censoring
  vector[N] rcens;  // right censor points for interval censoring
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int<lower=1> Kc;  // number of population-level effects after centering
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  array[N] int<lower=1> J_1;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  // data for group-level effects of ID 2
  int<lower=1> N_2;  // number of grouping levels
  int<lower=1> M_2;  // number of coefficients per level
  array[N] int<lower=1> J_2;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_2_1;
  // data for group-level effects of ID 3
  int<lower=1> N_3;  // number of grouping levels
  int<lower=1> M_3;  // number of coefficients per level
  array[N] int<lower=1> J_3;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_3_1;
  // data for group-level effects of ID 4
  int<lower=1> N_4;  // number of grouping levels
  int<lower=1> M_4;  // number of coefficients per level
  array[N] int<lower=1> J_4;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_4_1;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] b;  // regression coefficients
  real<lower=0> Intercept;  // temporary intercept for centered predictors
  real<lower=0> sigma;  // dispersion parameter
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  array[M_1] vector[N_1] z_1;  // standardized group-level effects
  vector<lower=0>[M_2] sd_2;  // group-level standard deviations
  array[M_2] vector[N_2] z_2;  // standardized group-level effects
  vector<lower=0>[M_3] sd_3;  // group-level standard deviations
  array[M_3] vector[N_3] z_3;  // standardized group-level effects
  vector<lower=0>[M_4] sd_4;  // group-level standard deviations
  array[M_4] vector[N_4] z_4;  // standardized group-level effects
}
transformed parameters {
  vector[N_1] r_1_1;  // actual group-level effects
  vector[N_2] r_2_1;  // actual group-level effects
  vector[N_3] r_3_1;  // actual group-level effects
  vector[N_4] r_4_1;  // actual group-level effects
  real lprior = 0;  // prior contributions to the log posterior
  r_1_1 = (sd_1[1] * (z_1[1]));
  r_2_1 = (sd_2[1] * (z_2[1]));
  r_3_1 = (sd_3[1] * (z_3[1]));
  r_4_1 = (sd_4[1] * (z_4[1]));
  lprior += normal_lpdf(b | 0,25);
  lprior += gamma_lpdf(Intercept | 3.65, 0.01);
  lprior += normal_lpdf(sigma | 0,15)
    - 1 * normal_lccdf(0 | 0,15);
  lprior += normal_lpdf(sd_1 | 0,9)
    - 1 * normal_lccdf(0 | 0,9);
  lprior += normal_lpdf(sd_2 | 0,9)
    - 1 * normal_lccdf(0 | 0,9);
  lprior += normal_lpdf(sd_3 | 0,9)
    - 1 * normal_lccdf(0 | 0,9);
  lprior += normal_lpdf(sd_4 | 0,9)
    - 1 * normal_lccdf(0 | 0,9);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    mu += Intercept + Xc * b;
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_2_1[J_2[n]] * Z_2_1[n] + r_3_1[J_3[n]] * Z_3_1[n] + r_4_1[J_4[n]] * Z_4_1[n];
    }
    for (n in 1:N) {
    // special treatment of censored data
      if (cens[n] == 0) {
        target += normal_lpdf(Y[n] | mu[n], sigma);
      } else if (cens[n] == 1) {
        target += normal_lccdf(Y[n] | mu[n], sigma);
      } else if (cens[n] == -1) {
        target += normal_lcdf(Y[n] | mu[n], sigma);
      } else if (cens[n] == 2) {
        target += log_diff_exp(
          normal_lcdf(rcens[n] | mu[n], sigma),
          normal_lcdf(Y[n] | mu[n], sigma)
        );
      }
    }
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(z_1[1]);
  target += std_normal_lpdf(z_2[1]);
  target += std_normal_lpdf(z_3[1]);
  target += std_normal_lpdf(z_4[1]);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
  // additionally sample draws from priors
  real prior_b = normal_rng(0,25);
  real prior_Intercept = gamma_rng(3.65,0.01);
  real prior_sigma = normal_rng(0,15);
  real prior_sd_1 = normal_rng(0,9);
  real prior_sd_2 = normal_rng(0,9);
  real prior_sd_3 = normal_rng(0,9);
  real prior_sd_4 = normal_rng(0,9);
  // use rejection sampling for truncated priors
  while (prior_Intercept < 0) {
    prior_Intercept = gamma_rng(3.65,0.01);
  }
  while (prior_sigma < 0) {
    prior_sigma = normal_rng(0,15);
  }
  while (prior_sd_1 < 0) {
    prior_sd_1 = normal_rng(0,9);
  }
  while (prior_sd_2 < 0) {
    prior_sd_2 = normal_rng(0,9);
  }
  while (prior_sd_3 < 0) {
    prior_sd_3 = normal_rng(0,9);
  }
  while (prior_sd_4 < 0) {
    prior_sd_4 = normal_rng(0,9);
  }
}
