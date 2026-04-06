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


gibbs_gaussian_mix <- function(n, n_iter, n_burn, C_init, mu_init, y, sigma_mu2){

  comp_time <- system.time({
  K <- length(mu_init)
  mu_hist <- matrix(data = NA, nrow = n_iter, ncol = K)
  C_hist <- array(data = NA, dim = c(n_iter, n, K))
  mu <- mu_init

  for(iter in 1:n_iter){
    
    C <- matrix(0, nrow = n, ncol = K)

    for(i in 1:n){
    
      log_w <- rep(NA, K)
      
      for(k in 1:K){
        log_w[k] <- -(1/2)*(y[i] - mu[k])^2  + (-1/(2*sigma_mu2)) * mu[k]^2
      }
      
      log_w <- log_w - max(log_w)
      
      w <- exp(log_w) / sum(exp(log_w))


      c_index <- sample(1:K, 1, prob = w)
      
      C[i,c_index] <- 1  
      
      
    }
    
    for(k in 1:K){
      n_k <- sum(C[,k])
      
      if(n_k > 0){
      
        y_bar_k <- t(y) %*% C[,k] / n_k
        
        mu_post_mean <- n_k * y_bar_k / (n_k + 1/sigma_mu2)
        
        mu_post_var <- 1/(n_k + 1/sigma_mu2)
        
        mu[k] <- rnorm(1, mean = mu_post_mean, sd = sqrt(mu_post_var))
      }
    }
    

    mu_hist[iter,] <- mu
    C_hist[iter, , ] <- C

    
  }  
  
  post_means <- colMeans(mu_hist[(n_burn+1):n_iter, ], na.rm = TRUE)
  mu_CIs <- apply(mu_hist[(n_burn+1):n_iter, ], MARGIN = 2, 
                  function(x) quantile(x, c(0.025, 0.975))
  )
  names(post_means) <- 1:K
  colnames(mu_CIs) <- 1:K
  sorted_post_means <- sort(post_means)
  sorted_mu_CIs <- mu_CIs[ , as.numeric(names(sorted_post_means))]
  post_labels <- rep(NA, length = n) 
  sorted_post_labels <- rep(NA, length = n) 
    for(i in 1:n){
      subj_slice <- C_hist[(n_burn+1):n_iter, i, ] #(niter-nburn) x k
      subj_summary <- colSums(subj_slice) # 1 x k
      if(length(which(subj_summary == max(subj_summary))) != 1){
       print(i)
       print(subj_summary)
      }
      post_labels[i] <- which(subj_summary == max(subj_summary))[1]
      sorted_post_labels[i] <- which(names(sorted_post_means) == post_labels[i])[1]
    }
  })


  results_df <- data.frame(t(sorted_post_means))
  names(results_df) <- paste0("mu_", 1:K)
  results_tib <- as_tibble(results_df) %>%
    mutate(comp_time = comp_time[3],
           labels = list(sorted_post_labels),
           CIs = list(sorted_mu_CIs))
  
  return(list(
    results_tib = results_tib)
  )

}
