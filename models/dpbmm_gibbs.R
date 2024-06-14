dpbmm_gibbs <- function(T, S, a_p, b_p, a_lambda, b_lambda, a_psi, b_psi,
                        u1, u2, dp_print_interval = 500, print_interval = 50) {
  # import necessary libraries
  library(dplyr)
  
  # load data
  N1 <- dim(T)[1]
  N2 <- dim(T)[2]
  N3 <- dim(T)[3]
  
  # store the values during the Gibbs sampling
  GAMMA_DP <- vector("numeric", length = S)
  COMPONENT_NUM <- vector("numeric", length = S)
  Y <- array(rep(NaN, S*N1*N3), dim = c(S, N1, N3))
  PI_MAT <- array(rep(NaN, S*N1*N3), dim = c(S, N1, N3))
  LAMBDA <- matrix(NA, nrow = S, ncol = N2)
  PSI <- matrix(NA, nrow = S, ncol = N2)
  
  # initialize model parameters
  gamma_dp <- rgamma(1, shape = u1, rate = u2)
  y <- matrix(rbinom(N1*N3, 1, a_p/(a_p+b_p)), N1, N3)
  pi_mat <- matrix(NA, nrow = N1, ncol = N3)
  lambda <- matrix(rbeta(n = N2, shape1 = a_lambda, shape2 = b_lambda))
  psi <- matrix(rbeta(n = N2, shape1 = a_psi, shape2 = b_psi))
  z <- vector("numeric", N1)
  z <- sample(1:2, size = N1, replace = TRUE)
  accepted_gamma_dp <- 0; log_s_gamma <- 0
  
  t1 <- Sys.time()
  # Gibbs sampling
  for (t in 1:(2*S)) {
    # -------------------------------------------------------------------------- #
    ## sample each recording's latent membership z_i, i = 1, ..., N_1
    for (i in 1:N1) {
      z[i] <- 0
      component_list <- z[z > 0] %>% unique()
      log_prob_z <- c()
      
      for (r in component_list) {
        ### calculate the first term in z_i's full conditional
        n_r_exclude <- (z == r) %>% sum()
        log_first_term <- n_r_exclude %>% log()
        # cat(r, "term1: ", log_first_term, "\n")
        
        ### calculate the second term in z_i's full conditional
        n_r_k_exclude_1_list <- ((z == r) & (y == 1)) %>% colSums()
        n_r_k_exclude_0_list <- ((z == r) & (y == 0)) %>% colSums()
        log_second_term <- (
          mapply(beta, 
                 a_p + n_r_k_exclude_1_list + y[i, ], 
                 b_p + n_r_k_exclude_0_list + (1 - y[i, ])) %>% log() -
            mapply(beta, 
                   a_p + n_r_k_exclude_1_list, 
                   b_p + n_r_k_exclude_0_list) %>% log()
        ) %>% sum()
        # cat(r, "term2: ", log_second_term, "\n")
        
        log_prob_z <- c(log_prob_z, log_first_term + log_second_term)
      }
      
      ### generate potential new component for audio recording i
      log_first_term <- gamma_dp %>% log()
      log_second_term <- (
        mapply(beta, a_p + y[i, ], b_p + (1 - y[i, ])) %>% log() - 
          mapply(beta, a_p, b_p) %>% log()
      ) %>% sum()
      
      # cat("new's term1: ", log_first_term, "\n")
      # cat("new's term2:", log_second_term, "\n")
      log_prob_z <- c(log_prob_z, log_first_term + log_second_term)
      
      ### normalize
      prob_z <- (log_prob_z - max(log_prob_z)) %>% exp()
      prob_z <- prob_z / sum(prob_z)
      
      ### draw recording's latent membership
      index <- sample(length(prob_z), size = 1, prob = prob_z)
      if (index > length(component_list)) {
        # create a new component
        z[i] <- max(component_list) + 1
      } else {
        z[i] <- component_list[index]
      }
      
      if (i %% dp_print_interval == 0) {
        cat(paste0("Iteration ", t, " - ", "Sample ", i), "Updated...\n")
        cat("The number of components: ", length(unique(z)), "\n")
        cat("The weights of components: ", table(z) / N1, "\n")
      }
    }
    
    
    # -------------------------------------------------------------------------- #
    ## sample each recording's latent indicator of all bird species y_i, i = 1, ..., N_1
    for (i in 1:N1) {
      y[i, ] <- -1      # set each element of ith row to -1
      n_z_k_exclude_1_list <- ((z == z[i]) & (y == 1)) %>% colSums()
      n_z_k_exclude_0_list <- ((z == z[i]) & (y == 0)) %>% colSums()
      
      ### calculate the prob terms when y_{i, k} = 1, k = 1,..., N_3 (vectorized)
      log_first_term <- T[i, , ] * as.vector(log(lambda)) + 
        (1 - T[i, , ]) * as.vector(log(1 - lambda))
      log_first_term <- log_first_term %>% colSums(na.rm = TRUE)      # in the form of a vector
      log_second_term <- mapply(beta, 
                                a_p + n_z_k_exclude_1_list + 1, 
                                b_p + n_z_k_exclude_0_list + 0) %>% log()
      
      log_prob_one <- log_first_term + log_second_term
      
      ### calculate the prob terms when y_{i, k} = 0, k = 1,..., N_3 (vectorized)
      log_first_term <- T[i, , ] * as.vector(log(psi)) + 
        (1 - T[i, , ]) * as.vector(log(1 - psi))
      log_first_term <- log_first_term %>% colSums(na.rm = TRUE)      # in the form of a vector
      log_second_term <- mapply(beta, 
                                a_p + n_z_k_exclude_1_list + 0, 
                                b_p + n_z_k_exclude_0_list + 1) %>% log()
      
      log_prob_zero <- log_first_term + log_second_term
      
      ### normalize & sample
      tmp_log_prob <- pmax(log_prob_one, log_prob_zero)
      log_prob_one <- log_prob_one - tmp_log_prob
      log_prob_zero <- log_prob_zero - tmp_log_prob
      
      prob_one <- exp(log_prob_one) / (exp(log_prob_one) + exp(log_prob_zero))
      
      pi_mat[i, ] <- prob_one
      y[i, ] <- rbinom(length(prob_one), size = 1, prob = prob_one)
    }
    

    # -------------------------------------------------------------------------- #
    ## sample \lambda_j & psi_j, j = 1, ..., N_2
    for (j in 1:N2) {
      # lambda
      a_lambda_t <- a_lambda + (y * T[, j, ]) %>% sum(na.rm = TRUE)
      b_lambda_t <- b_lambda + (y * (1 - T[, j, ])) %>% sum(na.rm = TRUE)
      lambda[j] <- rbeta(1, a_lambda_t, b_lambda_t)
      
      # psi
      a_psi_t <- a_psi + ((1 - y) * T[, j, ]) %>% sum(na.rm = TRUE)
      b_psi_t <- b_psi + ((1 - y) * (1 - T[, j, ])) %>% sum(na.rm = TRUE)
      psi[j] <- rbeta(1, a_psi_t, b_psi_t)
    }
    
    
    # -------------------------------------------------------------------------- #
    ## sample \gamma (concentration parameter in DP) via metropolis-within-gibbs
    if (t %% 50 == 0) {
      if ((accepted_gamma_dp / 50) > 0.44) {
        log_s_gamma <- log_s_gamma + min(0.01, sqrt(1 / t))
      } else {
        log_s_gamma <- log_s_gamma - min(0.01, sqrt(1 / t))
      }
      accepted_gamma_dp <- 0
    }
    
    new_gamma_dp <- gamma_dp + rnorm(1, 0, exp(log_s_gamma))
    
    if ((new_gamma_dp > 0) & (new_gamma_dp < 50^2)) {
      R <- z %>% unique() %>% length()
      alpha_gamma_dp <- (new_gamma_dp / gamma_dp)^(R + u1 - 1) * exp(-u2 * (new_gamma_dp - gamma_dp))
      alpha_gamma_dp <- alpha_gamma_dp * exp(
        sum(log(seq(0, N1 - 1) + gamma_dp)) - sum(log(seq(0, N1 - 1) + new_gamma_dp))
      )
      alpha_gamma_dp <- min(1, alpha_gamma_dp)
      
      if (runif(1) < alpha_gamma_dp) {
        gamma_dp <- new_gamma_dp
        accepted_gamma_dp <- accepted_gamma_dp + 1
      }
    }
    
    
    # -------------------------------------------------------------------------- #
    # store the posterior after burn-in
    if (t > S) {
      COMPONENT_NUM[t - S] <- length(unique(z))
      GAMMA_DP[t - S] <- gamma_dp
      
      PI_MAT[t - S, , ] <- pi_mat
      Y[t - S, , ] <- y
      
      LAMBDA[t - S, ] <- lambda
      PSI[t - S, ] <- psi
    }
    
    if (t %% print_interval == 0) {
      t2 <- Sys.time()
      cat(paste0("Iteration ", t, ": ", difftime(t2, t1, units = "secs"), " seconds\n"))
      t1 <- t2
    }
  }
  
  return(list(COMPONENT_NUM = COMPONENT_NUM, PI_MAT = PI_MAT, Y = Y,
              LAMBDA = LAMBDA, PSI = PSI, GAMMA_DP = GAMMA_DP))
}
