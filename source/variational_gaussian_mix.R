var_gauss_mix <- function(
    max_iter,
    tol_criteria,
    s2_init,
    m_init,
    phi_init,
    data,
    sigma_mu2){
  
  comp_time <- system.time(
    {
      
  K <- length(m_init)
  n <- length(data)
  m_hist <- matrix(NA, nrow = max_iter, ncol = K)
  s2_hist <- matrix(NA, nrow = max_iter, ncol = K)
  phi_hist <- array(data = NA, dim = c(max_iter, n, K))
  ELBO_hist <- rep(NA, max_iter)
  
  y <- data

  
  m_hist[1,] <- m_init   
  s2_hist[1,] <- s2_init
  phi_hist[1, , ] <- phi_init
  ELBO_hist[1] <- 0
  
  m <- m_init #- y_center #) / y_scale 
  s2 <- s2_init #/ y_scale^2
  phi <- phi_init
  

  tol = 100
  iter = 2
  
  while(iter <= max_iter & abs(tol) > tol_criteria){ #"abs" not correct

    ELBO_iter <- 0
    
    for(i in 1:n){
      log_phi <- y[i]*m - (m^2 + s2) / 2
      log_phi <- log_phi - max(log_phi)
      phi[i,] <- exp(log_phi)
      phi[i,] <- phi[i,] / sum(phi[i,])
    }
    

    for(k in 1:K){
      s2[k] <- 1 / (sum(phi[,k]) + 1 / sigma_mu2)
      m[k] <- sum(y * phi[,k]) / (sum(phi[,k]) + 1 / sigma_mu2)
      
      
      ELBO_t1_stor <- rep(NA, n)
      ELBO_t4_stor <- rep(NA, n)
      for(i in 1:n){
        ELBO_t1_stor[i] <- sum(phi[i,k] * (y^2 - 2*y*m[k] + m[k]^2 + s2[k]))
        ELBO_t4_stor[i] <- sum(log(phi[i,k]))
      }
      
      
      ELBO_t1 <- sum(ELBO_t1_stor)
      ELBO_t4 <- sum(ELBO_t4_stor)
      
      ELBO_iter <- ELBO_iter + 
        ELBO_t1 - 
        (1 / (2*sigma_mu2))*(m[k]^2 + s2[k]) + 
        (1/2) * (log(2*pi*s2[k]) + 1) -
        ELBO_t4
      


    }
    

    
    ELBO_hist[iter] <- ELBO_iter
    phi_hist[iter, , ] <- phi
    s2_hist[iter, ] <- s2
    m_hist[iter, ] <- m
    
    ELBO_iter <- 0
    tol <- ELBO_hist[iter] - ELBO_hist[iter-1] 
    iter <- iter + 1
  }
  print(m)
  
  labels <- apply(phi_hist[iter-1, , ], MARGIN = 1, function(x){
    which(x == max(x))
    } )
  })
  
  #m <- m * y_scale + y_center
  #s2 <- s2 * y_scale^2
  
  return(list(
    means = m,
    variances = s2,
    labels = labels,
    n_iter = iter,
    comp_time = comp_time[3]
  ))
  
}

#ELBO not monotonically increasing?

