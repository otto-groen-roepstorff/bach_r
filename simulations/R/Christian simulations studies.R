library(glmnet)

n_sims <- 250
simulation <- matrix(data = NA, ncol = 5, nrow = n_sims)
ba_t <- -1
tau <- 3
n <- 1000

start.time <- Sys.time()
for (j in 1:n_sims){
  data <- generate_survival_data(n, ba_t = ba_t, bx_t = log(2), bz_t = log(2),
                                 ba_c = 1, bx_c = log(2), bz_c = log(2), prop_a = 1/3, seed = sample(1:100000))
  
  #data_for_breslow <- generate_survival_data(n, ba_t = ba_t, bx_t = log(2), bz_t = log(2),
  #                                           ba_c = 1, bx_c = log(2), bz_c = log(2), seed = sample(1:1000000))
  
  #data_for_cum_haz <- generate_survival_data(n, ba_t = ba_t, bx_t = log(2), bz_t = log(2),
  #                                           ba_c = 1, bx_c = log(2), bz_c = log(2), seed = sample(1:1000000))
  
  model_surv <- oracle_model(data)
  model_cens <- oracle_cens_model(data)
  
  covar_true <- get_oracle_covar(data)
  covar_cens <- get_oracle_cens_covar(data)
  
  
  props <- propensity(data)
  
    # Function to compute HAL coefficients, survival probabilities, and cumulative baseline hazard for survival analysis
  hal_survival <- function(X, y, censor_status, lambda1, lambda2, jump_times) {
    # Standardize the features
    X_std <- scale(X)
    
    # Fit a Cox proportional hazards model with HAL penalty
    fit <- cv.glmnet(X_std, Surv(y, event = censor_status), alpha = 1, family = "cox",
                     lambda = NULL, standardize = FALSE, thresh = 1e-10)
    
    # Extract the coefficients
    beta <- as.vector(coef(fit, s = lambda1))
    
    # Apply HAL penalty
    beta <- beta / (1 + lambda2 * abs(beta))
    
    # Compute survival probabilities at jump times
    if (length(jump_times) > 0) {
      surv_prob <- rep(NA, length(jump_times))
      for (i in seq_along(jump_times)) {
        t <- jump_times[i]
        new_data <- as.matrix(cbind(1, scale(c(t, rep(0, ncol(X) - 1)), center = fit$offset, scale = fit$scale)))
        linear_predictor <- sum(new_data * beta)
        surv_prob[i] <- exp(-exp(linear_predictor))
      }
    } else {
      surv_prob <- NULL
    }
    
    # Compute cumulative baseline hazard
    baseline_hazard <- basehaz(fit)
    
    return(list(coefficients = beta, survival_probabilities = surv_prob, baseline_hazard = baseline_hazard$hazard))
  }
  
  
  hal_survival(data, data$T_true, 
               censor_status = data$Uncensored, 
               lambda1 = 1, 
               lambda2 = 1, 
               jump_times = get_jump_times(get_cum(model_surv)))
  
  # Example usage:
  # Assuming X is your feature matrix, y is the survival time, and jump_times is a vector of time points
  # result <- hal_survival(X, y, lambda1 = 0.01, lambda2 = 0.01, jump_times = c(50, 100, 150))
  # print("Optimized Coefficients:", result$coefficients)
  # print("Survival Probabilities at Jump Times:", result$survival_probabilities)
  # print("Cumulative Baseline Hazard:", result$baseline_hazard)
  
  
  
  cum_matrix_obs <- get_cum(model_surv)
  jump_times_obs <- get_jump_times(cum_matrix_obs)
  cum_bas_haz_obs <- get_cum_base_haz(cum_matrix_obs)
  no_jumps_obs <- get_row_length(cum_matrix_obs)
  beta_hat_obs <- get_param_est(model_surv)
  
  
  #breslow <- breslow_estimator(data_for_breslow, beta_hat_obs, T_corr = T)
  #cum_bas_haz_obs_breslow <- breslow$breslow
  #breslow_jump_times <- breslow$jump_times_breslow
  #no_jump_times_breslow <- length(breslow_jump_times)
  #
  #cum_bas_haz_obs = rep(NA,no_jumps_obs) 
  #for(i in 1:no_jumps_obs){
  #  cum_bas_haz_obs[i] = cum_bas_haz_obs_breslow[max((1:no_jump_times_breslow)[breslow_jump_times<=jump_times_obs[i]])]
  #}
  
  #plot(jump_times_obs, cum_bas_haz_obs, type = 'l')
  #lines(jump_times_obs, cum_bas_haz_obs_old)
  #lines(breslow_jump_times, cum_bas_haz_obs_breslow)
  
  covar_A0 <- covar_true %>% mutate(A = 0)
  covar_A1 <- covar_true %>% mutate(A = 1)
  
  #Cum haz obs
  cum_hat_obs <- predict_cox.aalen(covar = covar_true, betaHat = beta_hat_obs, cum_base_haz = cum_bas_haz_obs)
  cum_hat0_obs <- predict_cox.aalen(covar = covar_A0, betaHat = beta_hat_obs, cum_base_haz = cum_bas_haz_obs)
  cum_hat1_obs <- predict_cox.aalen(covar = covar_A1, betaHat = beta_hat_obs, cum_base_haz = cum_bas_haz_obs)
  
  #Survival functions
  S_hat_obs <- exp(-cum_hat_obs)
  S_hat0_obs <- exp(-cum_hat0_obs)
  S_hat1_obs <- exp(-cum_hat1_obs)
  
  
  #Censoring
  cum_matrix_cens <- get_cum(model_cens)
  jump_times_cens <- get_jump_times(cum_matrix_cens)
  cum_bas_haz_cens <- get_cum_base_haz(cum_matrix_cens)
  no_jumps_cens <- get_row_length(cum_matrix_cens)
  beta_hat_cens <- get_param_est(model_cens)
  
  
  #breslow_cens <- breslow_estimator_cens(data_for_breslow, beta_hat_cens, cens_corr = T)
  #cum_bas_haz_obs_breslow_cens <- breslow_cens$breslow
  #breslow_jump_times_cens <- breslow_cens$jump_times_breslow
  #no_jump_times_breslow_cens <- length(breslow_jump_times_cens)
  
  #plot(jump_times_cens, cum_bas_haz_cens, type = 'l')
  #lines(breslow_jump_times_cens, cum_bas_haz_obs_breslow_cens)
  #lines(jump_times_cens, cum_bas_haz_cens)
  
  #cum_bas_haz_cens = rep(NA,no_jumps_cens) 
  #for(i in 1:no_jumps_cens){
  #  cum_bas_haz_cens[i] = cum_bas_haz_obs_breslow_cens[max((1:no_jump_times_breslow_cens)[breslow_jump_times_cens<=jump_times_cens[i]])]
  #}
  
  
  cum_hat_cens <- predict_cox.aalen(covar = covar_cens, betaHat = beta_hat_cens, cum_base_haz = cum_bas_haz_cens)
  
  S_hat_cens <- exp(- cum_hat_cens)
  
  
  Khat_cens = matrix(NA, nrow = n, ncol = no_jumps_obs) 
  for(i in 1:no_jumps_obs){
    Khat_cens[,i] = S_hat_cens[,max((1:no_jumps_cens)[jump_times_cens<=jump_times_obs[i]])]
  }
  
  #plot(jump_times_obs, colMeans(Khat_cens), type = 'l')
  #lines(jump_times_obs, Khat_cens[1,])
  #lines(jump_times_cens, colMeans(S_hat_cens))
  
  
  #Martingales
  T_obs <- get_observed_times(data)
  
  dN <- outer(T_obs, jump_times_obs, '==')
  
  at_risk <- outer(T_obs, jump_times_obs, '>=')
  dL <- cbind(0, cum_hat_obs[,-1] - cum_hat_obs[,-ncol(cum_hat_obs)])
  
  dM <- dN - at_risk * dL
  
  par(mfrow = c(1,1))
  plot(jump_times_obs, cumsum(colSums(dM)), type = 'l', ylab = 'Kummuleret dM.')
  plot(jump_times_obs, cumsum(colSums(dN)), type = 'l', ylab = 'Kummuleret dN., dL.', col = 'red')
  lines(jump_times_obs, cumsum(colSums(at_risk*dL)), col = alpha('blue', 0.2))
  
  #plot(jump_times_obs, cumsum(dM[233,]), type = 'l', ylim = c(-3,1.1))
  #lines(jump_times_obs, cumsum(dN[233,]))
  #lines(jump_times_obs, -cumsum(at_risk[233,]*dL[233,]))
  
  #EIF
  no_times_to_keep <- sum(jump_times_obs <= tau)
  
  d_cum_base_haz <- c(0, cum_bas_haz_obs[-1]-cum_bas_haz_obs[-no_jumps_obs])
  d_cum_base_haz_tau <- c(0, cum_bas_haz_obs[-1]-cum_bas_haz_obs[-no_jumps_obs])[1:no_times_to_keep]
  
  S_hat_obs_tau <- S_hat_obs[,1:no_times_to_keep]
  S_hat0_obs_tau <- S_hat0_obs[,1:no_times_to_keep]
  S_hat1_obs_tau <- S_hat1_obs[,1:no_times_to_keep]
  
  Khat_cens_tau <- Khat_cens[,1:no_times_to_keep]
  
  dM_tau <- dM[,1:no_times_to_keep]
  dN_tau <- dN[,1:no_times_to_keep]
  dLambda_tau <- (at_risk * dL)[,1:no_times_to_keep]
  
  pi_0 <- props$propens$pi_a_0
  pi_1 <- props$propens$pi_a_1
  
  
  #Part 1
  p1 <- pi_0 * rowSums((S_hat1_obs_tau/Khat_cens_tau)*dM_tau)
  
  
  #Part 2
  prop_sum <- pi_0 + pi_1
  exp_0_multiplier <- exp(data.matrix(covar_A0) %*% beta_hat_obs)
  
  #Inner integrand
  inner_integrand <- dM_tau/(S_hat_obs_tau*Khat_cens_tau)
  inner_integral <- t(apply(inner_integrand, MARGIN = 1, FUN = cumsum))
  
  #Outer integrand
  outer_integrand <- t(t(S_hat0_obs_tau * S_hat1_obs_tau) * d_cum_base_haz_tau)
  
  p2 <- - prop_sum * exp_0_multiplier * rowSums(outer_integrand * inner_integral)
  
  #Part 3 (naive estimate)
  p3 <- exp_0_multiplier * rowSums(outer_integrand)
  
  p4 <- mean(p3)
  
  EIF_var_est <- var(p1 + p2 + p3 - p4)
  
  
  simulation[j,1] <- mean(p1)
  simulation[j,2] <- mean(p2)
  simulation[j,3] <- mean(p3)
  simulation[j,4] <- mean(p1 + p2 + p3)
  simulation[j,5] <- EIF_var_est
  print(paste0('Simulation ',j,' done'))
}
end.time <- Sys.time()
end.time - start.time
mean(simulation[,4])
mean(simulation[,3])
mean(simulation[,5]) 


#True value: 0.7215708370

#Nice misspecified survival
#0.6912407 - 300, 0.6980863 - 600, 0.7030793 - 900, 0.7076685 - 1500, 0.7089811 - 2400
#0.6891559 - 300, 0.6974173 - 600, 0.7030534 - 900, 0.7085803 - 1500, 0.7096222 - 2400

#Correctly specified survival
#0.6949226 - 300, 0.7058112 - 600, 0.7149363 - 1200, 0.7190154 - 2400, 0.7198531 - 4800
#0.6942408 - 300, 0.7049965 - 600, 0.7155028 - 1200, 0.7188588 - 2400, 0.7196702 - 4800

#Nasty misspecified survival
# 0.7509969 - 300, 0.7569597 - 600, 0.7608548 - 900, 0.7675732 - 2400 
# 0.7449133 - 300, 0.7504189 - 600, 0.7551813 - 900, 0.7619871 - 2400

#Correctly specified survival and misspecified censoring
# 0.6960505 - 300, 0.7077723 - 600, 0.7137625 - 900, 0.7182949 - 2400 
# 0.6944582 - 300, 0.7083863 - 600, 0.7135256 - 900, 0.7187663 - 2400



#simulations_df$sim12 <- c(n, 
#                          tau, 
#                          ba_t, 
#                          'False', 
#                          'True', 
#                          mean(simulation[,4]), 
#                          mean(simulation[,3]),
#                          mean(simulation[,1]) + mean(simulation[,2]))
#

#simulations_df <- data.frame('Specifikationer' = 
#                               c('Antal individer', 
#                                 'tau', 
#                                 'Beta_A', 
#                                 'Surv_corr', 
#                                 'Cens_corr', 
#                                 'Estimat',
#                                 'Naivt estimat',
#                                 'Bias correction'))
#








