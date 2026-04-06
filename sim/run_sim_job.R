# Will go back and forth between sampling c_i and mu_k. Sample all c_is at once or 
# one at a time? 
# nt total iterations

K <- 4
n <- 100
n_iter <- 10000
n_burn <- 2000
#C_init <- sample(1:K, n, replace = TRUE)
sigma_mu2 <- 100


# sim data
set.seed(12)
mu_true <- c(0, 5, 10, 20)
C_label_true <- sample(1:K, n, replace = TRUE)
y <- as.numeric(lapply(1:n, function(i) { 
  rnorm(n = 1, 
        mean = mu_true[C_label_true[i]],
        sd = sqrt(1)
  )
}))

#needs to be a n x K matrix with 1s in kth component with initial guess
# of which mixture y came from or actually just vector of length n with components in 
# 1 - k

mu_init <- c(0, 0, 0, 0)

one_sim <- gibbs_gaussian_mix(
  n = n, 
  n_iter = n_iter, 
  n_burn = n_burn, 
  mu_init = mu_init, 
  y = y, 
  sigma_mu2 = sigma_mu2)




