## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----run example model--------------------------------------------------------
library(rstan)
library(bennu)
library(bayesplot)
library(ggplot2)

rstan_options(auto_write = TRUE)
options(mc.cores = 2)

## basic example code
d <- generate_model_data()
# note iter should be at least 2000 to generate a reasonable posterior sample
fit <- est_naloxone(d,iter=100)
mcmc_pairs(fit, pars = c("sigma","mu0","zeta"),
           off_diag_args = list(size = 1, alpha = 0.5))

## ----summarize draws----------------------------------------------------------
library(posterior)

summarise_draws(fit, default_summary_measures())

## ----plot draws---------------------------------------------------------------
plot_kit_use(model = fit,data=d)

## ----run prior----------------------------------------------------------------
# note iter should be at least 2000 to generate a reasonable posterior sample
prior_fit <- est_naloxone(d,
                          run_estimation = FALSE,
                          iter=100)
mcmc_pairs(prior_fit, pars = c("sigma","mu0","zeta"),
           off_diag_args = list(size = 1, alpha = 0.5))

## ----plot compare prior posterior---------------------------------------------
plot_kit_use(prior = prior_fit,posterior = fit, data=d)

## ----rw order 2 model---------------------------------------------------------
# note iter should be at least 2000 to generate a reasonable posterior sample
# note use `rw_type` to specify order of random walk
rw2_fit <- est_naloxone(d,
                    rw_type = 2,
                    iter=100)
mcmc_pairs(rw2_fit, pars = c("sigma","mu0","zeta"),
           off_diag_args = list(size = 1, alpha = 0.5))

## ----rw order 2 summary draws-------------------------------------------------
summarise_draws(rw2_fit, default_convergence_measures())

## ----plot draws comparison----------------------------------------------------
plot_kit_use(rw_1 = fit,rw_2 = rw2_fit,data=d)

## ----data generation lower frequency------------------------------------------
## basic example code
d_missing <- generate_model_data(reporting_freq = 3)
ggplot(aes(x=times,y=Reported_Used,color=as.factor(regions)),data=d_missing) +
  geom_point()


## ----missing data inference---------------------------------------------------
missing_fit <- est_naloxone(d_missing,
                    rw_type = 2,
                    iter=100)
mcmc_pairs(missing_fit, pars = c("sigma","mu0","zeta"),
           off_diag_args = list(size = 1, alpha = 0.5))


## ----plot draws with missing data---------------------------------------------
plot_kit_use(model = missing_fit, data=d_missing)

## ----example change units of data, eval = FALSE-------------------------------
#  # monthly reporting delay distribution
#  psi_vec <- c(0.7, 0.2, 0.1)
#  # convert to weeks using interpolation
#  weekly_psi_vec <- rep(psi_vec,4) / sum(psi_vec)
#  
#  # properties of order to distributed delay distribution in months
#  max_delays <- 3
#  delay_alpha <- 2
#  delay_beta <- 1
#  
#  # convert to weeks
#  weekly_max_delays <- max_delays*4
#  weekly_delay_alpha <- delay_alpha
#  weekly_delay_beta <- 0.25 * delay_beta
#  
#  result <- est_naloxone(d,
#   psi_vec = weekly_psi_vec,
#   max_delays = weekly_max_delays,
#   delay_alpha = weekly_delay_alpha,
#   delay_beta = weekly_delay_beta)
#  

