base_gibbs <- function(T, S, a_p, b_p, a_lambda, b_lambda, a_psi, b_psi,
                       print_interval = 50) {
  # import necessary libraries
  library(dplyr)
  
  # load data
  N1 <- dim(T)[1]
  N2 <- dim(T)[2]
  N3 <- dim(T)[3]
  
  # store the values during the Gibbs sampling
  P <- matrix(NA, nrow = S, ncol = N3)
  LAMBDA <- matrix(NA, nrow = S, ncol = N2)
  PSI <- matrix(NA, nrow = S, ncol = N2)
  Y <- array(rep(NA, S*N1*N3), dim = c(S, N1, N3))
  P_MAT <- array(rep(NA, S*N1*N3), dim = c(S, N1, N3))
  
  # initialize model parameters
  p <- matrix(rbeta(n = N3, shape1 = a_p, shape2 = b_p))
  y <- matrix(rbinom(N1*N3, 1, a_p/(a_p+b_p)), N1, N3)
  p_mat <- matrix(NA, nrow = N1, ncol = N3)
  lambda <- matrix(rbeta(n = N2, shape1 = a_lambda, shape2 = b_lambda))
  psi <- matrix(rbeta(n = N2, shape1 = a_psi, shape2 = b_psi))
  
  t1 <- Sys.time()
  # Gibbs sampling
  for (t in 1:(2*S)) {
    # -------------------------------------------------------------------------- #
    ## update p_{k}, k = 1,...,N_3
    p <- rbeta(N3, shape1 = colSums(y) + a_p,
               shape2 = N1 - colSums(y) + b_p)
    
    
    # -------------------------------------------------------------------------- #
    # update y_{ik}, i = 1,...,N_1; k = 1,...,N_3
    for (i in 1:N1) {
      for (k in 1:N3) {
        prob_one <- p[k]*(((lambda**T[i, ,k])*((1-lambda)**(1-T[i, ,k]))) %>% prod(na.rm = TRUE))
        prob_zero <- (1 - p[k])*(((psi**T[i, ,k])*((1-psi)**(1-T[i, ,k]))) %>% prod(na.rm = TRUE))
        p_ik <- prob_one / (prob_one + prob_zero)
        y[i, k] <- rbinom(1, 1, p_ik)
        # keep track of the values of p_ik is even more valuable than y_ik
        p_mat[i, k] <- p_ik
      }
    }
    
    
    # -------------------------------------------------------------------------- #
    # update lambda_j & psi_j, j = 1,...,N_2
    for (j in 1:N2) {
      # lambda
      a_lambda_t <- a_lambda + (y*T[, j, ]) %>% sum(na.rm = TRUE)
      b_lambda_t <- b_lambda + (y*(1-T[, j, ])) %>% sum(na.rm = TRUE)
      lambda[j] <- rbeta(1, a_lambda_t, b_lambda_t)
      
      # psi
      a_psi_t <- a_psi + ((1-y)*T[, j, ]) %>% sum(na.rm = TRUE)
      b_psi_t <- b_psi + ((1-y)*(1-T[, j, ])) %>% sum(na.rm = TRUE)
      psi[j] <- rbeta(1, a_psi_t, b_psi_t)
    }
    
    
    # -------------------------------------------------------------------------- #
    # store the posterior samples after burn-in
    if (t > S) {
      P[t - S, ] <- p

      Y[t - S, , ] <- y
      P_MAT[t - S, , ] <- p_mat
      
      LAMBDA[t - S, ] <- lambda
      PSI[t - S, ] <- psi
    }
    
    if (t %% print_interval == 0) {
      t2 <- Sys.time()
      cat(paste0("Iteration ", t, ": ", difftime(t2, t1, units = "secs"), " seconds\n"))
      t1 <- t2
    }
  }
  
  return(list(P = P, LAMBDA = LAMBDA, PSI = PSI,
              Y = Y, P_MAT = P_MAT))
}
