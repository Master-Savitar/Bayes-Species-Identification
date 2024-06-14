base_hierachy_gibbs <- function(T, rater_species_mat, S, a_p, b_p, mu_lambda, phi_lambda, phi_lambda_star,
                                mu_psi, phi_psi, phi_psi_star, burnin = 2000, print_interval = 50) {
  # import necessary libraries
  library(dplyr)
  
  # define related function
  sigmoid <- function(x) {
    return(1 / (1 + exp(-x)))
  }
  
  # load data
  N1 <- dim(T)[1]
  N2 <- dim(T)[2]
  N3 <- dim(T)[3]
  
  # store the values during the Gibbs sampling
  P <- matrix(NA, nrow = S, ncol = N3)
  LAMBDA <- matrix(NA, nrow = S, ncol = N2)
  LMABDA_MAT <- array(rep(NaN, S*N2*N3), dim = c(S, N2, N3))
  PSI <- matrix(NA, nrow = S, ncol = N2)
  PSI_MAT <- array(rep(NaN, S*N2*N3), dim = c(S, N2, N3))
  Y <- array(rep(NA, S*N1*N3), dim = c(S, N1, N3))
  P_MAT <- array(rep(NA, S*N1*N3), dim = c(S, N1, N3))
  
  # initialize model parameters
  p <- matrix(rbeta(n = N3, shape1 = a_p, shape2 = b_p))
  y <- matrix(rbinom(N1*N3, 1, a_p / (a_p+b_p)), N1, N3)
  p_mat <- matrix(NA, nrow = N1, ncol = N3)
  lambda <- matrix(rnorm(n = N2, mean = mu_lambda, sd = phi_lambda))
  lambda_mat <- matrix(rnorm(
    N2*N3, mean = matrix(rep(lambda, each = N3), nrow = N2, ncol = N3, byrow = TRUE), 
    sd = phi_lambda_star), nrow = N2, ncol = N3)
  lambda_mat[rater_species_mat == 0] <- NA
  psi <- matrix(rnorm(n = N2, mean = mu_psi, sd = phi_psi))
  psi_mat <- matrix(rnorm(
    N2*N3, mean = matrix(rep(psi, each = N3), nrow = N2, ncol = N3, byrow = TRUE), 
    sd = phi_psi_star), nrow = N2, ncol = N3)
  psi_mat[rater_species_mat == 0] <- NA
  
  t1 <- Sys.time()
  # Gibbs sampling
  for (t in 1:(S+burnin)) {
    # -------------------------------------------------------------------------- #
    ## update p_{k}, k = 1,...,N_3
    p <- rbeta(N3, shape1 = colSums(y) + a_p,
               shape2 = N1 - colSums(y) + b_p)
    
    
    # -------------------------------------------------------------------------- #
    # update y_{ik}, i = 1,...,N_1; k = 1,...,N_3
    for (i in 1:N1) {
      for (k in 1:N3) {
        prob_one <- p[k]*(((sigmoid(lambda_mat[, k])**T[i, , k]) * 
                             ((1 - sigmoid(lambda_mat[, k]))**(1 - T[i, , k]))) %>% prod(na.rm = TRUE))
        prob_zero <- (1 - p[k])*(((sigmoid(psi_mat[, k])**T[i, , k]) * 
                                    ((1 - sigmoid(psi_mat[, k]))**(1 - T[i, , k]))) %>% prod(na.rm = TRUE))
        p_ik <- prob_one / (prob_one + prob_zero)
        y[i, k] <- rbinom(1, 1, p_ik)
        # keep track of the values of p_ik is even more valuable than y_ik
        p_mat[i, k] <- p_ik
      }
    }
    
    
    # -------------------------------------------------------------------------- #
    ## sample \lambda_1, \lambda_2, \dots, \lambda_{N_2} together
    tau_lambda_t <- 1 / (1 / phi_lambda^2 + N3 / phi_lambda_star^2)
    mu_lambda_t <- (mu_lambda/phi_lambda^2 + (rowSums(lambda_mat, na.rm = TRUE))/phi_lambda_star^2) * tau_lambda_t
    lambda <- matrix(rnorm(n = N2, mean = mu_lambda_t, sd = sqrt(tau_lambda_t)))
    
    ## sample \lambda_{j, k}, j = 1, \dots, N_2, k = 1, \dots, N_3
    delta_lambda_p <- (aperm(array(y, dim = c(dim(y), dim(T)[2])), c(1, 3, 2)) * T) %>% 
      apply(MARGIN = c(2, 3), FUN = sum, na.rm = TRUE)
    delta_lambda_n <- (aperm(array(y, dim = c(dim(y), dim(T)[2])), c(1, 3, 2)) * (1 - T)) %>% 
      apply(MARGIN = c(2, 3), FUN = sum, na.rm = TRUE)
    
    # mu_lambda_mat_t <- matrix(NA, nrow = N2, ncol = N3)
    # tau_lambda_mat_t <- matrix(NA, nrow = N2, ncol = N3)
    for (j in 1:N2) {
      for (k in 1:N3) {
        if (rater_species_mat[j, k] == 1) {
          pg_p <- rpg(num = delta_lambda_p[j, k], 1, lambda_mat[j, k]) %>% sum()
          pg_n <- rpg(num = delta_lambda_n[j, k], 1, lambda_mat[j, k]) %>% sum()
          tau_lambda_jk <- 1 / (1/phi_lambda_star^2 + pg_p + pg_n)
          mean_lambda_jk <- ((lambda[j]/phi_lambda_star^2)+(delta_lambda_p[j, k]-delta_lambda_n[j, k])/2)*tau_lambda_jk
          lambda_mat[j, k] <- rnorm(n = 1, mean = mean_lambda_jk, sd = sqrt(tau_lambda_jk))
        }
      }
    }
    phi_lambda_star <- sqrt((((lambda_mat - lambda %>% as.vector())^2) %>% sum(na.rm = TRUE)) / sum(rater_species_mat))
    # lambda_mat <- matrix(rnorm(n = N2*N3, mean = mu_lambda_mat_t, sd = sqrt(tau_lambda_mat_t)),
    #                      nrow = N2, ncol = N3)
    
    
    # -------------------------------------------------------------------------- #
    ## sample \psi_1, \psi_2, \dots, \psi_{N_2} together
    tau_psi_t <- 1 / (1 / phi_psi^2 + N3 / phi_psi_star^2)
    mu_psi_t <- (mu_psi/phi_psi^2+(rowSums(psi_mat, na.rm = TRUE))/phi_psi_star^2) * tau_psi_t
    psi <- matrix(rnorm(n = N2, mean = mu_psi_t, sd = sqrt(tau_psi_t)))
    
    ## sample \psi_{j, k}, j = 1, \dots, N_2, k = 1, \dots, N_3
    delta_psi_p <- (aperm(array(1 - y, dim = c(dim(y), dim(T)[2])), c(1, 3, 2)) * T) %>% 
      apply(MARGIN = c(2, 3), FUN = sum, na.rm = TRUE)
    delta_psi_n <- (aperm(array(1 - y, dim = c(dim(y), dim(T)[2])), c(1, 3, 2)) * (1 - T)) %>% 
      apply(MARGIN = c(2, 3), FUN = sum, na.rm = TRUE)
    
    # mu_psi_mat_t <- matrix(NA, nrow = N2, ncol = N3)
    # tau_psi_mat_t <- matrix(NA, nrow = N2, ncol = N3)
    for (j in 1:N2) {
      for (k in 1:N3) {
        if (rater_species_mat[j, k] == 1) {
          pg_p <- rpg(num = delta_psi_p[j, k], 1, psi_mat[j, k]) %>% sum()
          pg_n <- rpg(num = delta_psi_n[j, k], 1, psi_mat[j, k]) %>% sum()
          tau_psi_jk <- 1 / (1/phi_psi_star^2 + pg_p + pg_n)
          mean_psi_jk <- ((psi[j]/phi_psi_star^2)+(delta_psi_p[j, k]-delta_psi_n[j, k])/2)*tau_psi_jk
          psi_mat[j, k] <- rnorm(n = 1, mean = mean_psi_jk, sd = sqrt(tau_psi_jk)) 
        }
      }
    }
    phi_psi_star <- sqrt((((psi_mat - psi %>% as.vector())^2) %>% sum(na.rm = TRUE)) / sum(rater_species_mat))
    # psi_mat <- matrix(rnorm(n = N2*N3, mean = mu_psi_mat_t, sd = sqrt(tau_psi_mat_t)),
    #                   nrow = N2, ncol = N3)
    
    
    # -------------------------------------------------------------------------- #
    # store the posterior samples after burn-in
    if (t > burnin) {
      P[t - burnin, ] <- p

      Y[t - burnin, , ] <- y
      P_MAT[t - burnin, , ] <- p_mat

      LAMBDA[t - burnin, ] <- lambda
      LMABDA_MAT[t - burnin, , ] <- lambda_mat

      PSI[t - burnin, ] <- psi
      PSI_MAT[t - burnin, , ] <- psi_mat
    }
    
    
    if (t %% print_interval == 0) {
      t2 <- Sys.time()
      cat(paste0("Iteration ", t, ": ", difftime(t2, t1, units = "secs"), " seconds\n"))
      cat("phi_lambda_star: ", phi_lambda_star, "\n")
      cat("phi_psi_star: ", phi_psi_star, "\n")
      t1 <- t2
    }
  }
  
  return(list(P = P, LAMBDA = LAMBDA, LMABDA_MAT = LMABDA_MAT,
              PSI = PSI, PSI_MAT = PSI_MAT, Y = Y, P_MAT = P_MAT))
}
