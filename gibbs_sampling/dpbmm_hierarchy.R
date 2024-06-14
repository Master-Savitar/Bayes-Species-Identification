dpbmm_hierarchy_gibbs <- function(T, rater_species_mat, S, a_p, b_p, mu_lambda, phi_lambda, phi_lambda_star,
                                  mu_psi, phi_psi, phi_psi_star, u1, u2, burnin = 2000, 
                                  dp_print_interval = 500, print_interval = 50) {
  # import necessary libraries
  library(dplyr)
  
  # define related functiona
  sigmoid <- function(x) {
    return(1 / (1 + exp(-x)))
  }
  
  # load data
  N1 <- dim(T)[1]
  N2 <- dim(T)[2]
  N3 <- dim(T)[3]
  
  # store the values during the Gibbs sampling
  GAMMA_DP <- vector("numeric", length = S)
  COMPONENT_NUM <- vector("numeric", length = S)
  Y <- array(rep(NaN, S*N1*N3), dim = c(S, N1, N3))
  P_MAT <- array(rep(NaN, S*N1*N3), dim = c(S, N1, N3))
  LAMBDA <- matrix(NA, nrow = S, ncol = N2)
  LMABDA_MAT <- array(rep(NaN, S*N2*N3), dim = c(S, N2, N3))
  PSI <- matrix(NA, nrow = S, ncol = N2)
  PSI_MAT <- array(rep(NaN, S*N2*N3), dim = c(S, N2, N3))
  
  # initialize model parameters
  gamma_dp <- rgamma(1, shape = u1, rate = u2)
  y <- matrix(rbinom(N1*N3, 1, a_p/(a_p+b_p)), N1, N3)
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
  z <- vector("numeric", N1)
  z <- sample(1:2, size = N1, replace = TRUE)
  accepted_gamma_dp <- 0; log_s_gamma <- 0
  
  t1 <- Sys.time()
  # Gibbs sampling
  for (t in 1:S) {
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
        # cat("The concentration parameter: ", gamma_dp, "\n")
      }
    }
    
    
    # -------------------------------------------------------------------------- #
    ## sample each recording's latent indicator of all bird species y_i, i = 1, ..., N_1
    for (i in 1:N1) {
      y[i, ] <- -1      # set each element of ith row to -1
      n_z_k_exclude_1_list <- ((z == z[i]) & (y == 1)) %>% colSums()
      n_z_k_exclude_0_list <- ((z == z[i]) & (y == 0)) %>% colSums()
      
      ### calculate the prob terms when y_{i, k} = 1, k = 1,..., N_3 (vectorized)
      log_first_term <- T[i, , ] * as.vector(log(sigmoid(lambda_mat))) + 
        (1 - T[i, , ]) * as.vector(log(1 - sigmoid(lambda_mat)))
      log_first_term <- log_first_term %>% colSums(na.rm = TRUE)      # in the form of a vector
      log_second_term <- mapply(beta, 
                                a_p + n_z_k_exclude_1_list + 1, 
                                b_p + n_z_k_exclude_0_list + 0) %>% log()
      
      log_prob_one <- log_first_term + log_second_term
      
      ### calculate the prob terms when y_{i, k} = 0, k = 1,..., N_3 (vectorized)
      log_first_term <- T[i, , ] * as.vector(log(sigmoid(psi_mat))) + 
        (1 - T[i, , ]) * as.vector(log(1 - sigmoid(psi_mat)))
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
      
      p_mat[i, ] <- prob_one
      y[i, ] <- rbinom(length(prob_one), size = 1, prob = prob_one)
    }
    

    # -------------------------------------------------------------------------- #
    ## sample \lambda_1, \lambda_2, \dots, \lambda_{N_2} together
    tau_lambda_t <- 1 / (1 / phi_lambda^2 + N3 / phi_lambda_star^2)
    mu_lambda_t <- (mu_lambda/phi_lambda^2+(rowSums(lambda_mat, na.rm = TRUE))/phi_lambda_star^2) * tau_lambda_t
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
    
    
    # store the posterior after burn-in
    if (t > burnin) {
      COMPONENT_NUM[t - burnin] <- length(unique(z))

      P_MAT[t - burnin, , ] <- p_mat
      Y[t - burnin, , ] <- y

      LAMBDA[t - burnin, ] <- lambda
      LMABDA_MAT[t - burnin, , ] <- lambda_mat

      PSI[t - burnin, ] <- psi
      PSI_MAT[t - burnin, , ] <- psi_mat

      GAMMA_DP[t - burnin] <- gamma_dp
    }
    
    
    if (t %% print_interval == 0) {
      t2 <- Sys.time()
      cat(paste0("Iteration ", t, ": ", difftime(t2, t1, units = "secs"), " seconds\n"))
      cat("The concentration parameter: ", gamma_dp, "\n")
      t1 <- t2
    }
  }
  
  return(list(COMPONENT_NUM = COMPONENT_NUM, P_MAT = P_MAT, Y = Y,
              LAMBDA = LAMBDA, LMABDA_MAT = LMABDA_MAT,
              PSI = PSI, PSI_MAT = PSI_MAT, GAMMA_DP = GAMMA_DP,
              z = z, accepted_gamma_dp = accepted_gamma_dp, log_s_gamma = log_s_gamma))
}
