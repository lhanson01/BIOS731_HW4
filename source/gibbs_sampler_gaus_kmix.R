
gibbs_gaussian_mix <- function(n_iter, n_burn, mu_init, y, sigma_mu2, n_chains){

  
    chain_results <- vector(mode = "list", length = n_chains)
    mu_hist_chains <- vector(mode = "list", length = n_chains)
    
    for(chain in 1:n_chains){
    comp_time <- system.time({
    n <- length(y)
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
        post_labels[i] <- which(subj_summary == max(subj_summary))[1]
        sorted_post_labels[i] <- which(names(sorted_post_means) == post_labels[i])[1]
      }
    })
    results_df <- data.frame(t(sorted_post_means))
    names(results_df) <- paste0("mu_", 1:K)
    chain_results[[chain]] <- as_tibble(results_df) %>%
      mutate(comp_time = comp_time[3],
             labels = list(sorted_post_labels),
             CIs = list(sorted_mu_CIs))
    mu_hist_chains[[chain]] <- mu_hist
  }


  results_tib <- bind_rows(chain_results)
  
  return(list(
    results_tib = results_tib,
    mu_hist = mu_hist_chains)
  )

}
