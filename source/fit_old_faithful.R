wait_times <- faithful$waiting
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
  n = n,
  max_iter = max_iter,
  tol_criteria = 1e-5,
  s2_init = s2_init,
  m_init = mu_init,
  phi_init = phi_init,
  y = wait_times,
  sigma_mu2 = sigma_mu2
)
