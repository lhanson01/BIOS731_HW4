sim_gauss_mix <- function(n, K, mu_true, C_true){
  y <- as.numeric(lapply(seq_len(n), function(i) { 
    rnorm(n = 1, 
          mean = mu_true[C_true[i]],
          sd = sqrt(1)
    )
  })
  )
    return(y)
}
