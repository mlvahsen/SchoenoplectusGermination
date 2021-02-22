# Vahsen et al. Germination paper

# STAN Models 1-4: 1 = Binomial without zero-inflation, 2 = Zero-inflated
# binomial, 3 = Beta-binomial, 4 = Zero-inflated beta-binomial

# Also includes STAN models for Models 1-4 that also contain Assay_Num as a
# random effect. Results not shown for these, but Assay_Num did not explain
# considerable variance for all models so the term was dropped from the models.

# Model 1 (binomial model WITHOUT zero inflation) ##
# Fixed effects (binomial): age, temperature, media, pretreatment
# Random effects (binomial): location

# Load libraries
library(here)

model1_stan <- "
data{
int N; //number of rows in the dataset
int Nbeta; //number of coefficients for the polynomial calibration
int L; //number of locations
int<lower=0> y[N]; //number of successfully germinated seeds
int s[N]; //number of seeds planted
real x[N]; //depth
int t1[N]; //maps row to temperature 1
int t2[N]; //maps row to temperature 2
int t3[N]; //maps row to temperature 3
int m1[N]; //maps row to media 1
int m2[N]; //maps row to media 2
int l[N]; //maps row to location
int h[N]; //maps row to pre-treatment
int k1[N]; //maps row to photoperiod 1
int k2[N]; //maps row to photoperiod 2
vector[Nbeta] prior_b_mean; // hyperparameter mean for calibration
matrix[Nbeta,Nbeta] prior_b_covar; // hyperparameter covariance for calibration
real prior_sigma_alpha; // hyperparameter alpha for calibration error
real prior_sigma_beta; // hyperparameter beta for calibration error
}

parameters{
real b0;
real b1;
real b3[L];
real b4; // coefficient for m1
real b5; // coefficient for m2
real b6; // coefficient for t1
real b7; // coefficient for t2
real b8; // coefficient for t3
real b9; // coefficient for h
real b10; // coefficient for k1
real b11; // coefficient for k2
real <lower=0> sigma3;
real sigma_pred;
real <lower=0> sigma_notch;
vector [Nbeta] b;
}

transformed parameters{
vector [N] xa;
vector [N] mu;
real bar_xa;
real xa_sd;
vector [N] xa_std;

for (n in 1:N){
xa[n] = b[1] + b[2]*x[n] + b[3]*x[n]^2 + sigma_pred;
}

bar_xa = mean(xa);
xa_sd = sd(xa);
xa_std = (xa - bar_xa)/xa_sd;

for (n in 1:N){
mu[n] = inv_logit(b0+b1*xa_std[n]+b3[l[n]]+b4*m1[n]+b5*m2[n]+b6*t1[n]+b7*t2[n]+b8*t3[n]+b9*h[n]+b10*k1[n]+b11*k2[n]);
}
}

model{
b0~normal(0, 10);
b1~normal(0, 10);
b3~normal(0, sigma3);
b4~normal(0, 10);
b5~normal(0, 10);
b6~normal(0, 10);
b7~normal(0, 10);
b8~normal(0, 10);
b9~normal(0, 10);
b10~normal(0, 10);
b11~normal(0, 10);
b ~ multi_normal(prior_b_mean, prior_b_covar);
sigma_pred ~ normal(0, sigma_notch);
sigma_notch ~ gamma(prior_sigma_alpha, prior_sigma_beta);

y ~ binomial(s, mu);
}

generated quantities{
vector [N] log_lik;
vector [N] y_new;

for (n in 1:N){
           y_new[n] = binomial_rng(s[n],mu[n]);
         }

for (n in 1:N){
log_lik[n] = binomial_lpmf(y[n] |s[n], mu[n]);
}
}
"



# Model 2 (full model with zero inflation) ##
# Fixed effects (binomial): age, temperature, media, pretreatment
# Fixed effects (bernoulli): age
# Random effects (binomial): location

model2_stan <- "
data{
int N; //number of rows in the dataset
int Nbeta; //number of coefficients for the polynomial calibration
int L; //number of locations
int<lower=0> y[N]; //number of successfully germinated seeds
int s[N]; //number of seeds planted
real x[N]; //depth
int t1[N]; //maps row to temperature 1
int t2[N]; //maps row to temperature 2
int t3[N]; //maps row to temperature 3
int m1[N]; //maps row to media 1
int m2[N]; //maps row to media 2
int l[N]; //maps row to location
int h[N]; //maps row to pre-treatment
int k1[N]; //maps row to photoperiod 1
int k2[N]; //maps row to photoperiod 2
vector[Nbeta] prior_b_mean; // hyperparameter mean for calibration
matrix[Nbeta,Nbeta] prior_b_covar; // hyperparameter covariance for calibration
real prior_sigma_alpha; // hyperparameter alpha for calibration error
real prior_sigma_beta; // hyperparameter beta for calibration error
}

parameters{
real b0;
real b1;
real b3[L];
real b4; // coefficient for m1
real b5; // coefficient for m2
real b6; // coefficient for t1
real b7; // coefficient for t2
real b8; // coefficient for t3
real b9; // coefficient for h
real b10; // coefficient for k1
real b11; // coefficient for k2
real g0;
real g1;
real <lower=0> sigma3;
real sigma_pred;
real <lower=0> sigma_notch;
vector [Nbeta] b;
}

transformed parameters{
vector [N] xa;
vector [N] mu;
vector [N] theta;
real bar_xa;
real xa_sd;
vector [N] xa_std;

for (n in 1:N){
xa[n] = b[1] + b[2]*x[n] + b[3]*x[n]^2 + sigma_pred;
}

bar_xa = mean(xa);
xa_sd = sd(xa);
xa_std = (xa - bar_xa)/xa_sd;

for (n in 1:N){
mu[n]=exp(b0+b1*xa_std[n]+b3[l[n]]+b4*m1[n]+b5*m2[n]+b6*t1[n]+b7*t2[n]+b8*t3[n]+b9*h[n]+b10*k1[n]+b11*k2[n])/(1+exp(b0+b1*xa_std[n]+b3[l[n]]+b4*m1[n]+b5*m2[n]+b6*t1[n]+b7*t2[n]+b8*t3[n]+b9*h[n]+b10*k1[n]+b11*k2[n]));
theta[n]=exp(g0+ g1*xa_std[n])/(1+exp(g0+ g1*xa_std[n]));
}
}

model{
b0~normal(0, 10);
b1~normal(0, 10);
b3~normal(0, sigma3);
b4~normal(0, 10);
b5~normal(0, 10);
b6~normal(0, 10);
b7~normal(0, 10);
b8~normal(0, 10);
b9~normal(0, 10);
b10~normal(0, 10);
b11~normal(0, 10);
g0~normal(0, 10);
g1~normal(0, 10);
b ~ multi_normal(prior_b_mean, prior_b_covar);
sigma_pred ~ normal(0, sigma_notch);
sigma_notch ~ gamma(prior_sigma_alpha, prior_sigma_beta);

for (n in 1:N) {
if (y[n] == 0)
target += log_sum_exp(bernoulli_lpmf(1 | theta[n]),
bernoulli_lpmf(0 | theta[n]) + binomial_lpmf(y[n] |s[n], mu[n]));
else
target += bernoulli_lpmf(0 | theta[n]) + binomial_lpmf(y[n] | s[n], mu[n]);
}
}

generated quantities{
vector [N] log_lik;
vector [N] y_new;

for (n in 1:N){
if (bernoulli_rng(theta[n]) == 1) {
             y_new[n] = 0;
         } else {
           y_new[n] = binomial_rng(s[n],mu[n]);
             }
         }

for (n in 1:N){ if (y[n]==0)
log_lik[n] = log_sum_exp(bernoulli_lpmf(1 | theta[n]),
bernoulli_lpmf(0 | theta[n]) + binomial_lpmf(y[n] |s[n], mu[n]));
else
log_lik[n] = bernoulli_lpmf(0 | theta[n]) + binomial_lpmf(y[n] | s[n], mu[n]);	}	}
"


# Beta-binomial
model3_stan <- "

functions {

  real beta_binomial2_lpmf(int y, real mu, real phi, int T) {
    return beta_binomial_lpmf(y | T, mu * phi, (1 - mu) * phi);
  }
  int beta_binomial2_rng(real mu, real phi, int T) {
    return beta_binomial_rng(T, mu * phi, (1 - mu) * phi);
  }

}

data{
  int N; //number of rows in the dataset
  int Nbeta; //number of coefficients for the polynomial calibration
  int L; //number of locations
  int<lower=0> y[N]; //number of successfully germinated seeds
  int s[N]; //number of seeds planted
  real x[N]; //depth
  int t1[N]; //maps row to temperature 1
  int t2[N]; //maps row to temperature 2
  int t3[N]; //maps row to temperature 3
  int m1[N]; //maps row to media 1
  int m2[N]; //maps row to media 2
  int l[N]; //maps row to location
  int h[N]; //maps row to pre-treatment
  int k1[N]; //maps row to photoperiod 1
  int k2[N]; //maps row to photoperiod 2
  vector[Nbeta] prior_b_mean; // hyperparameter mean for calibration
  matrix[Nbeta,Nbeta] prior_b_covar; // hyperparameter covariance for calibration
  real prior_sigma_alpha; // hyperparameter alpha for calibration error
  real prior_sigma_beta; // hyperparameter beta for calibration error
}

parameters{
  real <lower=0> phi;
  real b0;
  real b1;
  real b3[L];
  real b4; // coefficient for m1
  real b5; // coefficient for m2
  real b6; // coefficient for t1
  real b7; // coefficient for t2
  real b8; // coefficient for t3
  real b9; // coefficient for h
  real b10; // coefficient for k1
  real b11; // coefficient for k2
  real <lower=0> sigma3;
  real sigma_pred;
  real <lower=0> sigma_notch;
  vector [Nbeta] b;
}

transformed parameters{
  vector [N] xa;
  vector [N] mu;
  real bar_xa;
  real xa_sd;
  vector [N] xa_std;
  
  for (n in 1:N){
    xa[n] = b[1] + b[2]*x[n] + b[3]*x[n]^2 + sigma_pred;
  }
  
  bar_xa = mean(xa);
  xa_sd = sd(xa);
  xa_std = (xa - bar_xa)/xa_sd;
  
  for (n in 1:N){
    mu[n] = inv_logit(b0+b1*xa_std[n]+b3[l[n]]+b4*m1[n]+b5*m2[n]+b6*t1[n]+b7*t2[n]+b8*t3[n]+b9*h[n]+b10*k1[n]+b11*k2[n]);
  }
  
}

model{
  b0~normal(0, 10);
  b1~normal(0, 10);
  b3~normal(0, sigma3);
  b4~normal(0, 10);
  b5~normal(0, 10);
  b6~normal(0, 10);
  b7~normal(0, 10);
  b8~normal(0, 10);
  b9~normal(0, 10);
  b10~normal(0, 10);
  b11~normal(0, 10);
  phi~gamma(0.01, 0.01);
  b ~ multi_normal(prior_b_mean, prior_b_covar);
  sigma_pred ~ normal(0, sigma_notch);
  sigma_notch ~ gamma(prior_sigma_alpha, prior_sigma_beta);
  
  for (n in 1:N) {
      target += beta_binomial2_lpmf(y[n] | mu[n], phi, s[n]);
  }
  
}

generated quantities{
  vector [N] y_new;
  vector [N] log_lik;
  
  for (n in 1:N){
    y_new[n] = beta_binomial2_rng(mu[n], phi, s[n]); 
  } 
  
  for (n in 1:N){
    log_lik[n] = beta_binomial2_lpmf(y[n] | mu[n], phi, s[n]);
 }
  

  }"

model4_stan <- "
functions {

  real beta_binomial2_lpmf(int y, real mu, real phi, int T) {
    return beta_binomial_lpmf(y | T, mu * phi, (1 - mu) * phi);
  }
  int beta_binomial2_rng(real mu, real phi, int T) {
    return beta_binomial_rng(T, mu * phi, (1 - mu) * phi);
  }

}

data{
int N; //number of rows in the dataset
int Nbeta; //number of coefficients for the polynomial calibration
int L; //number of locations
int<lower=0> y[N]; //number of successfully germinated seeds
int s[N]; //number of seeds planted
real x[N]; //depth
int t1[N]; //maps row to temperature 1
int t2[N]; //maps row to temperature 2
int t3[N]; //maps row to temperature 3
int m1[N]; //maps row to media 1
int m2[N]; //maps row to media 2
int l[N]; //maps row to location
int h[N]; //maps row to pre-treatment
int k1[N]; //maps row to photoperiod 1
int k2[N]; //maps row to photoperiod 2
vector[Nbeta] prior_b_mean; // hyperparameter mean for calibration
matrix[Nbeta,Nbeta] prior_b_covar; // hyperparameter covariance for calibration
real prior_sigma_alpha; // hyperparameter alpha for calibration error
real prior_sigma_beta; // hyperparameter beta for calibration error
}

parameters{
real b0;
real b1;
real b3[L];
real b4; // coefficient for m1
real b5; // coefficient for m2
real b6; // coefficient for t1
real b7; // coefficient for t2
real b8; // coefficient for t3
real b9; // coefficient for h
real b10; // coefficient for k1
real b11; // coefficient for k2
real g0;
real g1;
real <lower=0> phi;
real <lower=0> sigma3;
real sigma_pred;
real <lower=0> sigma_notch;
vector [Nbeta] b;
}

transformed parameters{
vector [N] xa;
vector [N] mu;
vector [N] theta;
real bar_xa;
real xa_sd;
vector [N] xa_std;

for (n in 1:N){
xa[n] = b[1] + b[2]*x[n] + b[3]*x[n]^2 + sigma_pred;
}

bar_xa = mean(xa);
xa_sd = sd(xa);
xa_std = (xa - bar_xa)/xa_sd;

for (n in 1:N){
mu[n]=exp(b0+b1*xa_std[n]+b3[l[n]]+b4*m1[n]+b5*m2[n]+b6*t1[n]+b7*t2[n]+b8*t3[n]+b9*h[n]+b10*k1[n]+b11*k2[n])/(1+exp(b0+b1*xa_std[n]+b3[l[n]]+b4*m1[n]+b5*m2[n]+b6*t1[n]+b7*t2[n]+b8*t3[n]+b9*h[n]+b10*k1[n]+b11*k2[n]));
theta[n]=exp(g0+ g1*xa_std[n])/(1+exp(g0+ g1*xa_std[n]));
}
}

model{
b0~normal(0, 10);
b1~normal(0, 10);
b3~normal(0, sigma3);
b4~normal(0, 10);
b5~normal(0, 10);
b6~normal(0, 10);
b7~normal(0, 10);
b8~normal(0, 10);
b9~normal(0, 10);
b10~normal(0, 10);
b11~normal(0, 10);
g0~normal(0, 10);
g1~normal(0, 10);
b ~ multi_normal(prior_b_mean, prior_b_covar);
sigma_pred ~ normal(0, sigma_notch);
sigma_notch ~ gamma(prior_sigma_alpha, prior_sigma_beta);

for (n in 1:N) {
if (y[n] == 0)
target += log_sum_exp(bernoulli_lpmf(1 | theta[n]),
bernoulli_lpmf(0 | theta[n]) + beta_binomial2_lpmf(y[n] | mu[n], phi, s[n]));
else
target += bernoulli_lpmf(0 | theta[n]) + beta_binomial2_lpmf(y[n] | mu[n], phi, s[n]);
}
}

generated quantities{
vector [N] log_lik;
vector [N] y_new;

for (n in 1:N){
if (bernoulli_rng(theta[n]) == 1) {
             y_new[n] = 0;
         } else {
           y_new[n] = beta_binomial2_rng(mu[n], phi, s[n]);
             }
         }

for (n in 1:N){ if (y[n]==0)
log_lik[n] = log_sum_exp(bernoulli_lpmf(1 | theta[n]),
bernoulli_lpmf(0 | theta[n]) + beta_binomial2_lpmf(y[n] | mu[n], phi, s[n]));
else
log_lik[n] = bernoulli_lpmf(0 | theta[n]) + beta_binomial2_lpmf(y[n] | mu[n], phi, s[n]);	}	}
"

# Write stan objects that can be pulled into main script for model fitting
writeLines(model1_stan, here("outputs", "model1_stan.stan"))
writeLines(model2_stan, here("outputs", "model2_stan.stan"))
writeLines(model3_stan, here("outputs", "model3_stan.stan"))
writeLines(model4_stan, here("outputs", "model4_stan.stan"))


## Assay models #### 

# These models contain assay number (experiment) as a random intercept. Results
# not published.

# Model 1 (binomial model WITHOUT zero inflation) ##
# Fixed effects (binomial): age, temperature, media, pretreatment
# Random effects (binomial): location, experiment

model1_assay_stan <- "
data{
int N; //number of rows in the dataset
int Nbeta; //number of coefficients for the polynomial calibration
int L; //number of locations
int M; // number of experiments
int<lower=0> y[N]; //number of successfully germinated seeds
int s[N]; //number of seeds planted
real x[N]; //depth
int t1[N]; //maps row to temperature 1
int t2[N]; //maps row to temperature 2
int t3[N]; //maps row to temperature 3
int m1[N]; //maps row to media 1
int m2[N]; //maps row to media 2
int l[N]; //maps row to location
int m[N]; // maps row to experiment
int h[N]; //maps row to pre-treatment
int k1[N]; //maps row to photoperiod 1
int k2[N]; //maps row to photoperiod 2
vector[Nbeta] prior_b_mean; // hyperparameter mean for calibration
matrix[Nbeta,Nbeta] prior_b_covar; // hyperparameter covariance for calibration
real prior_sigma_alpha; // hyperparameter alpha for calibration error
real prior_sigma_beta; // hyperparameter beta for calibration error
}

parameters{
real b0;
real b1;
real b2[M];
real b3[L];
real b4; // coefficient for m1
real b5; // coefficient for m2
real b6; // coefficient for t1
real b7; // coefficient for t2
real b8; // coefficient for t3
real b9; // coefficient for h
real b10; // coefficient for k1
real b11; // coefficient for k2
real <lower=0> sigma2;
real <lower=0> sigma3;
real sigma_pred;
real <lower=0> sigma_notch;
vector [Nbeta] b;
}

transformed parameters{
vector [N] xa;
vector [N] mu;
real bar_xa;
real xa_sd;
vector [N] xa_std;

for (n in 1:N){
xa[n] = b[1] + b[2]*x[n] + b[3]*x[n]^2 + sigma_pred;
}

bar_xa = mean(xa);
xa_sd = sd(xa);
xa_std = (xa - bar_xa)/xa_sd;

for (n in 1:N){
mu[n] = inv_logit(b0+b1*xa_std[n]+b2[m[n]]+b3[l[n]]+b4*m1[n]+b5*m2[n]+b6*t1[n]+b7*t2[n]+b8*t3[n]+b9*h[n]+b10*k1[n]+b11*k2[n]);
}
}

model{
b0~normal(0, 10);
b1~normal(0, 10);
b2~normal(0, sigma2);
b3~normal(0, sigma3);
b4~normal(0, 10);
b5~normal(0, 10);
b6~normal(0, 10);
b7~normal(0, 10);
b8~normal(0, 10);
b9~normal(0, 10);
b10~normal(0, 10);
b11~normal(0, 10);
b ~ multi_normal(prior_b_mean, prior_b_covar);
sigma_pred ~ normal(0, sigma_notch);
sigma_notch ~ gamma(prior_sigma_alpha, prior_sigma_beta);

y ~ binomial(s, mu);
}

generated quantities{
vector [N] log_lik;
vector [N] y_new;

for (n in 1:N){
           y_new[n] = binomial_rng(s[n],mu[n]);
         }

for (n in 1:N){
log_lik[n] = binomial_lpmf(y[n] |s[n], mu[n]);
}
}
"


# Model 2 (full model with zero inflation) ##
# Fixed effects (binomial): age, temperature, media, pretreatment
# Fixed effects (bernoulli): age
# Random effects (binomial): location, experiment

model2_assay_stan <- "
data{
int N; //number of rows in the dataset
int Nbeta; //number of coefficients for the polynomial calibration
int L; //number of locations
int M; // number of experiments
int<lower=0> y[N]; //number of successfully germinated seeds
int s[N]; //number of seeds planted
real x[N]; //depth
int t1[N]; //maps row to temperature 1
int t2[N]; //maps row to temperature 2
int t3[N]; //maps row to temperature 3
int m1[N]; //maps row to media 1
int m2[N]; //maps row to media 2
int l[N]; //maps row to location
int m[N]; // maps row to experiment
int h[N]; //maps row to pre-treatment
int k1[N]; //maps row to photoperiod 1
int k2[N]; //maps row to photoperiod 2
vector[Nbeta] prior_b_mean; // hyperparameter mean for calibration
matrix[Nbeta,Nbeta] prior_b_covar; // hyperparameter covariance for calibration
real prior_sigma_alpha; // hyperparameter alpha for calibration error
real prior_sigma_beta; // hyperparameter beta for calibration error
}

parameters{
real b0;
real b1;
real b2[M];
real b3[L];
real b4; // coefficient for m1
real b5; // coefficient for m2
real b6; // coefficient for t1
real b7; // coefficient for t2
real b8; // coefficient for t3
real b9; // coefficient for h
real b10; // coefficient for k1
real b11; // coefficient for k2
real g0;
real g1;
real <lower=0> sigma2;
real <lower=0> sigma3;
real sigma_pred;
real <lower=0> sigma_notch;
vector [Nbeta] b;
}

transformed parameters{
vector [N] xa;
vector [N] mu;
vector [N] theta;
real bar_xa;
real xa_sd;
vector [N] xa_std;

for (n in 1:N){
xa[n] = b[1] + b[2]*x[n] + b[3]*x[n]^2 + sigma_pred;
}

bar_xa = mean(xa);
xa_sd = sd(xa);
xa_std = (xa - bar_xa)/xa_sd;

for (n in 1:N){
mu[n]=exp(b0+b1*xa_std[n]+b2[m[n]]+b3[l[n]]+b4*m1[n]+b5*m2[n]+b6*t1[n]+b7*t2[n]+b8*t3[n]+b9*h[n]+b10*k1[n]+b11*k2[n])/(1+exp(b0+b1*xa_std[n]+b2[m[n]]+b3[l[n]]+b4*m1[n]+b5*m2[n]+b6*t1[n]+b7*t2[n]+b8*t3[n]+b9*h[n]+b10*k1[n]+b11*k2[n]));
theta[n]=exp(g0+ g1*xa_std[n])/(1+exp(g0+ g1*xa_std[n]));
}
}

model{
b0~normal(0, 10);
b1~normal(0, 10);
b2~normal(0, sigma2);
b3~normal(0, sigma3);
b4~normal(0, 10);
b5~normal(0, 10);
b6~normal(0, 10);
b7~normal(0, 10);
b8~normal(0, 10);
b9~normal(0, 10);
b10~normal(0, 10);
b11~normal(0, 10);
g0~normal(0, 10);
g1~normal(0, 10);
b ~ multi_normal(prior_b_mean, prior_b_covar);
sigma_pred ~ normal(0, sigma_notch);
sigma_notch ~ gamma(prior_sigma_alpha, prior_sigma_beta);

for (n in 1:N) {
if (y[n] == 0)
target += log_sum_exp(bernoulli_lpmf(1 | theta[n]),
bernoulli_lpmf(0 | theta[n]) + binomial_lpmf(y[n] |s[n], mu[n]));
else
target += bernoulli_lpmf(0 | theta[n]) + binomial_lpmf(y[n] | s[n], mu[n]);
}
}

generated quantities{
vector [N] log_lik;
vector [N] y_new;

for (n in 1:N){
if (bernoulli_rng(theta[n]) == 1) {
             y_new[n] = 0;
         } else {
           y_new[n] = binomial_rng(s[n],mu[n]);
             }
         }

for (n in 1:N){ if (y[n]==0)
log_lik[n] = log_sum_exp(bernoulli_lpmf(1 | theta[n]),
bernoulli_lpmf(0 | theta[n]) + binomial_lpmf(y[n] |s[n], mu[n]));
else
log_lik[n] = bernoulli_lpmf(0 | theta[n]) + binomial_lpmf(y[n] | s[n], mu[n]);	}	}
"


# Beta-binomial
model3_assay_stan <- "

functions {

  real beta_binomial2_lpmf(int y, real mu, real phi, int T) {
    return beta_binomial_lpmf(y | T, mu * phi, (1 - mu) * phi);
  }
  int beta_binomial2_rng(real mu, real phi, int T) {
    return beta_binomial_rng(T, mu * phi, (1 - mu) * phi);
  }

}

data{
  int N; //number of rows in the dataset
  int Nbeta; //number of coefficients for the polynomial calibration
  int L; //number of locations
  int M; // number of experiments
  int<lower=0> y[N]; //number of successfully germinated seeds
  int s[N]; //number of seeds planted
  real x[N]; //depth
  int t1[N]; //maps row to temperature 1
  int t2[N]; //maps row to temperature 2
  int t3[N]; //maps row to temperature 3
  int m1[N]; //maps row to media 1
  int m2[N]; //maps row to media 2
  int l[N]; //maps row to location
  int m[N]; //maps row to experiment
  int h[N]; //maps row to pre-treatment
  int k1[N]; //maps row to photoperiod 1
  int k2[N]; //maps row to photoperiod 2
  vector[Nbeta] prior_b_mean; // hyperparameter mean for calibration
  matrix[Nbeta,Nbeta] prior_b_covar; // hyperparameter covariance for calibration
  real prior_sigma_alpha; // hyperparameter alpha for calibration error
  real prior_sigma_beta; // hyperparameter beta for calibration error
}

parameters{
  real <lower=0> phi;
  real b0;
  real b1;
  real b2[M];
  real b3[L];
  real b4; // coefficient for m1
  real b5; // coefficient for m2
  real b6; // coefficient for t1
  real b7; // coefficient for t2
  real b8; // coefficient for t3
  real b9; // coefficient for h
  real b10; // coefficient for k1
  real b11; // coefficient for k2
  real <lower=0> sigma2;
  real <lower=0> sigma3;
  real sigma_pred;
  real <lower=0> sigma_notch;
  vector [Nbeta] b;
}

transformed parameters{
  vector [N] xa;
  vector [N] mu;
  real bar_xa;
  real xa_sd;
  vector [N] xa_std;
  
  for (n in 1:N){
    xa[n] = b[1] + b[2]*x[n] + b[3]*x[n]^2 + sigma_pred;
  }
  
  bar_xa = mean(xa);
  xa_sd = sd(xa);
  xa_std = (xa - bar_xa)/xa_sd;
  
  for (n in 1:N){
    mu[n] = inv_logit(b0+b1*xa_std[n]+b2[m[n]]+b3[l[n]]+b4*m1[n]+b5*m2[n]+b6*t1[n]+b7*t2[n]+b8*t3[n]+b9*h[n]+b10*k1[n]+b11*k2[n]);
  }
  
}

model{
  b0~normal(0, 10);
  b1~normal(0, 10);
  b2~normal(0, sigma2);
  b3~normal(0, sigma3);
  b4~normal(0, 10);
  b5~normal(0, 10);
  b6~normal(0, 10);
  b7~normal(0, 10);
  b8~normal(0, 10);
  b9~normal(0, 10);
  b10~normal(0, 10);
  b11~normal(0, 10);
  phi~gamma(0.01, 0.01);
  b ~ multi_normal(prior_b_mean, prior_b_covar);
  sigma_pred ~ normal(0, sigma_notch);
  sigma_notch ~ gamma(prior_sigma_alpha, prior_sigma_beta);
  
  for (n in 1:N) {
      target += beta_binomial2_lpmf(y[n] | mu[n], phi, s[n]);
  }
  
}

generated quantities{
  vector [N] y_new;
  vector [N] log_lik;
  
  for (n in 1:N){
    y_new[n] = beta_binomial2_rng(mu[n], phi, s[n]); 
  } 
  
  for (n in 1:N){
    log_lik[n] = beta_binomial2_lpmf(y[n] | mu[n], phi, s[n]);
 }
  

  }"

model4_assay_stan <- "
functions {

  real beta_binomial2_lpmf(int y, real mu, real phi, int T) {
    return beta_binomial_lpmf(y | T, mu * phi, (1 - mu) * phi);
  }
  int beta_binomial2_rng(real mu, real phi, int T) {
    return beta_binomial_rng(T, mu * phi, (1 - mu) * phi);
  }

}

data{
int N; //number of rows in the dataset
int Nbeta; //number of coefficients for the polynomial calibration
int L; //number of locations
int M; //number of experiments
int<lower=0> y[N]; //number of successfully germinated seeds
int s[N]; //number of seeds planted
real x[N]; //depth
int t1[N]; //maps row to temperature 1
int t2[N]; //maps row to temperature 2
int t3[N]; //maps row to temperature 3
int m1[N]; //maps row to media 1
int m2[N]; //maps row to media 2
int l[N]; //maps row to location
int m[N]; //maps row to experiment
int h[N]; //maps row to pre-treatment
int k1[N]; //maps row to photoperiod 1
int k2[N]; //maps row to photoperiod 2
vector[Nbeta] prior_b_mean; // hyperparameter mean for calibration
matrix[Nbeta,Nbeta] prior_b_covar; // hyperparameter covariance for calibration
real prior_sigma_alpha; // hyperparameter alpha for calibration error
real prior_sigma_beta; // hyperparameter beta for calibration error
}

parameters{
real b0;
real b1;
real b2[M];
real b3[L];
real b4; // coefficient for m1
real b5; // coefficient for m2
real b6; // coefficient for t1
real b7; // coefficient for t2
real b8; // coefficient for t3
real b9; // coefficient for h
real b10; // coefficient for k1
real b11; // coefficient for k2
real g0;
real g1;
real <lower=0> phi;
real <lower=0> sigma2;
real <lower=0> sigma3;
real sigma_pred;
real <lower=0> sigma_notch;
vector [Nbeta] b;
}

transformed parameters{
vector [N] xa;
vector [N] mu;
vector [N] theta;
real bar_xa;
real xa_sd;
vector [N] xa_std;

for (n in 1:N){
xa[n] = b[1] + b[2]*x[n] + b[3]*x[n]^2 + sigma_pred;
}

bar_xa = mean(xa);
xa_sd = sd(xa);
xa_std = (xa - bar_xa)/xa_sd;

for (n in 1:N){
mu[n]=exp(b0+b1*xa_std[n]+b2[m[n]]+b3[l[n]]+b4*m1[n]+b5*m2[n]+b6*t1[n]+b7*t2[n]+b8*t3[n]+b9*h[n]+b10*k1[n]+b11*k2[n])/(1+exp(b0+b1*xa_std[n]+b2[m[n]]+b3[l[n]]+b4*m1[n]+b5*m2[n]+b6*t1[n]+b7*t2[n]+b8*t3[n]+b9*h[n]+b10*k1[n]+b11*k2[n]));
theta[n]=exp(g0+ g1*xa_std[n])/(1+exp(g0+ g1*xa_std[n]));
}
}

model{
b0~normal(0, 10);
b1~normal(0, 10);
b2~normal(0, sigma2);
b3~normal(0, sigma3);
b4~normal(0, 10);
b5~normal(0, 10);
b6~normal(0, 10);
b7~normal(0, 10);
b8~normal(0, 10);
b9~normal(0, 10);
b10~normal(0, 10);
b11~normal(0, 10);
g0~normal(0, 10);
g1~normal(0, 10);
b ~ multi_normal(prior_b_mean, prior_b_covar);
sigma_pred ~ normal(0, sigma_notch);
sigma_notch ~ gamma(prior_sigma_alpha, prior_sigma_beta);

for (n in 1:N) {
if (y[n] == 0)
target += log_sum_exp(bernoulli_lpmf(1 | theta[n]),
bernoulli_lpmf(0 | theta[n]) + beta_binomial2_lpmf(y[n] | mu[n], phi, s[n]));
else
target += bernoulli_lpmf(0 | theta[n]) + beta_binomial2_lpmf(y[n] | mu[n], phi, s[n]);
}
}

generated quantities{
vector [N] log_lik;
vector [N] y_new;

for (n in 1:N){
if (bernoulli_rng(theta[n]) == 1) {
             y_new[n] = 0;
         } else {
           y_new[n] = beta_binomial2_rng(mu[n], phi, s[n]);
             }
         }

for (n in 1:N){ if (y[n]==0)
log_lik[n] = log_sum_exp(bernoulli_lpmf(1 | theta[n]),
bernoulli_lpmf(0 | theta[n]) + beta_binomial2_lpmf(y[n] | mu[n], phi, s[n]));
else
log_lik[n] = bernoulli_lpmf(0 | theta[n]) + beta_binomial2_lpmf(y[n] | mu[n], phi, s[n]);	}	}
"

