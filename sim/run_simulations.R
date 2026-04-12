suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))

source(here::here("source","gibbs_sampler_gaus_kmix.R"))
source(here::here("source","variational_gaussian_mix.R"))
source(here::here("source","sim_gauss_mix.R"))

nsim <- 500
lower_sim <- seq(1, 500, by = 50)
upper_sim <- seq(50,500, by = 50)
sim_chunks <- cbind(lower_sim,upper_sim)

n_scen <- c(100, 1000, 10000)

job_table <- expand.grid(n = n_scen, chunk = 1:nrow(sim_chunks))
job_table$lower_sim <- sim_chunks[job_table$chunk, 1]
job_table$upper_sim <- sim_chunks[job_table$chunk, 2]
job_table$chunk <- NULL

job <- as.numeric(commandArgs(trailingOnly = TRUE))[1]
seed <- floor(runif(nsim, 1, 10000))
results <- as.list(rep(NA, nsim))
mu_true <- c(0, 5, 10, 20)
max_iter <- 10000
K <- length(mu_true)

lower_sim <- job_table$lower_sim[job]
upper_sim <- job_table$upper_sim[job]
results <- vector(mode = "list", length = 0)
for(sim in lower_sim:upper_sim){
  print(sim)
  set.seed(seed[sim])
  
  n <- job_table$n[job]
  C_true <- sample(1:K, n, replace = TRUE)
  y <- sim_gauss_mix(n = n, K = K, mu_true = mu_true, C_true = C_true)
  
  
  mu_init <- seq(min(y), max(y), length.out = K)
  print(mu_init)
  gibbs_results <- gibbs_gaussian_mix(n_iter = 10000,
                     n_burn = 2000,
                     mu_init = mu_init,
                     y = y,
                     sigma_mu2 = 100,
                     n_chains = 3)
  
  print("finished gibbs")
  
  s2_init <- c(5, 5, 5, 5)
  phi_init <- matrix(data = 1/K, 
                     nrow = n,
                     ncol = K)
  vi_results <- var_gauss_mix(
    max_iter = max_iter,
    tol_criteria = 1e-5,
    s2_init = s2_init,
    m_init = mu_init,
    phi_init = phi_init,
    data = y,
    sigma_mu2 = 100
  )
  
  print("finished VI")

  res_label <- paste0("n_", n, "_sim_", sim)
  results[[res_label]] <- list("Gibbs" = gibbs_results,
                            "VI" = vi_results,
                            "mu_init" = mu_init)
  
}

Date <- gsub("-", "", Sys.Date())
results_path <- file.path(here::here("results"))
if(!file.exists(results_path)){
  dir.create(results_path, Date)
}

filename <- paste0(results_path, "/",
                  "n_", n,"_sim_", lower_sim, "-", upper_sim, ".RDA")

save(results,
     file = filename)

