# Will go back and forth between sampling c_i and mu_k. Sample all c_is at once or 
# one at a time? 
# nt total iterations

K <- 4
n <- 100
n_iter <- 10000
n_burn <- 2000
#C_init <- sample(1:K, n, replace = TRUE)
sigma_mu2 <- 100
tol_criteria <- 1e-05
max_iter <- 10000


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

m_init <- runif(4, range(y)[1], range(y)[2]) # cannot have equal values, cannot have wild values
s2_init <- c(10, 5, 3, 4)
phi_init <- matrix(data = 1/K, 
                   nrow = n,
                   ncol = K)

var_gauss_mix <- function(
    n,
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
  print(K)
  m_hist <- matrix(NA, nrow = max_iter, ncol = K)
  s2_hist <- matrix(NA, nrow = max_iter, ncol = K)
  phi_hist <- array(data = NA, dim = c(max_iter, n, K))
  ELBO_hist <- rep(NA, max_iter)
  
  #y <- scale(data)
  #y_center <- attributes(y)$`scaled:center`
  #y_scale <- attributes(y)$`scaled:scale`
  
  m_hist[1,] <- m_init#-y_center) #/ y_scale   
  s2_hist[1,] <- s2_init #/ y_scale^2 
  phi_hist[1, , ] <- phi_init
  ELBO_hist[1] <- 0
  
  m <- m_init #-y_center) #/ y_scale 
  s2 <- s2_init #/ y_scale^2
  phi <- phi_init
  

  tol = 100
  iter = 2
  
  while(iter <= max_iter & abs(tol) > tol_criteria){ #"abs" not correct
    print(iter)
    
    ELBO_iter <- 0
    
    for(i in 1:n){
      phi[i,] <- exp(y[i]*m - (m^2 + s2) / 2)
      phi[i,] <- phi[i,] / sum(phi[i,])
    }
    
    for(k in 1:K){
      print(paste("k:",k))
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
      
      print(ELBO_iter)


    }
    

    
    
    ELBO_hist[iter] <- ELBO_iter
    phi_hist[iter, , ] <- phi
    s2_hist[iter, ] <- s2
    m_hist[iter, ] <- m
    
    ELBO_iter <- 0
    tol <- ELBO_hist[iter] - ELBO_hist[iter-1] 
    print(tol)
    iter <- iter + 1
  }
  
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
    n_iter = iter
  ))
  
}

#ELBO not monotonically increasing?

