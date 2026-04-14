library(coda)
library(dplyr)
library(ggplot2)
library(forecast)

wait_times <- faithful$waiting # need to scale data
K <- 2
n <- length(wait_times)
n_iter <- max_iter <- 10000
n_burn <- 2000
sigma_mu2 <- 100

mu_init <- c(50,80)
s2_init <- c(5,5)
phi_init <- matrix(data = 1/K, 
                   nrow = n,
                   ncol = K)

# Var Bayes
var_results <- var_gauss_mix(
  max_iter = max_iter,
  tol_criteria = 1e-5,
  s2_init = s2_init,
  m_init = mu_init,
  phi_init = phi_init,
  data = wait_times,
  sigma_mu2 = sigma_mu2
)

# Gibbs
n_chains <- 3
gibbs_results <- gibbs_gaussian_mix(
    n_iter = 10000,
    n_chains = n_chains,
    n_burn = 2000,
    mu_init = mu_init,
    y = wait_times,
    sigma_mu2 = 100
  )

saveRDS(var_results,
        file = "./output/faithful_var_results.rds")
saveRDS(gibbs_results,
        file = "./output/faithful_gibbs_results.rds")

