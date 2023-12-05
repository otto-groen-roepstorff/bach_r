

test <- function(){
  
  n <- 1000
  data <- generate_survival_data(n, ba_t = 0, bx_t = log(2), bz_t = log(2),
                                 ba_c = 1, bx_c = log(2), bz_c = -2, seed = sample(1:1000))
  #Changing all data to observed data
  data$Uncensored <- T
  data$T_obs <- data$T_true
  
  model_surv <- oracle_model(data)
  
  props <- propensity(data)
  pi_a_0 <- props$propens$pi_a_0
    #2*(1-get_A_val(data))#props$propens$pi_a_0
  
  pi_a_1 <- props$propens$pi_a_1 #2*get_A_val(data) #props$propens$pi_a_1
  pi_sum = pi_a_0 + pi_a_1
  
  
  cum_matrix_obs <- get_cum(model_surv)
  jump_times_obs <- get_jump_times(cum_matrix_obs)
  cum_bas_haz_obs <- get_cum_base_haz(cum_matrix_obs)
  no_jumps_obs <- get_row_length(cum_matrix_obs)
  beta_hat_obs <- get_param_est(model_surv)
  max_time <- 3
  max_time_filter <- jump_times_obs<max_time
  
  
  
  
  covar_true <- get_oracle_covar(data)
  covar_A0 <- data.matrix(covar_true %>% mutate(A = 0))
  covar_A1 <- covar_true %>% mutate(A = 1)
  
  #Cum haz obs
  cum_hat_obs <- predict_cox.aalen(covar = covar_true, betaHat = beta_hat_obs, cum_base_haz = cum_bas_haz_obs)
  cum_hat0_obs <- predict_cox.aalen(covar = covar_A0, betaHat = beta_hat_obs, cum_base_haz = cum_bas_haz_obs)
  cum_hat1_obs <- predict_cox.aalen(covar = covar_A1, betaHat = beta_hat_obs, cum_base_haz = cum_bas_haz_obs)
  
    #Survival functions
  S_hat_obs <- exp(-cum_hat_obs)
  S_hat0_obs <- exp(-cum_hat0_obs)
  S_hat1_obs <- exp(-cum_hat1_obs)
  
  
  #Finding martingales
  T_obs <- get_observed_times(data)
  at_risk <- outer(X = T_obs, Y = jump_times_obs, FUN = ">=")
  dN <- outer(X = T_obs, Y = jump_times_obs, FUN = "==")
  
  alpha_A_0 <- exp(covar_A0%*%beta_hat_obs)[,1]
  
  denom <- matrix(data = NA, nrow = n, ncol = no_jumps_obs)
  
  for (i in (1:no_jumps_obs)){
    denom[,i] <- S_hat0_obs[,i]*alpha_A_0
  }  
  
  sum_integrand <- pi_a_1*(at_risk/(S_hat1_obs)-1)+pi_a_0*(dN/denom -1)
  
  
  
  #Cumulative hazard Y*alpha
  
  dH <- cbind(0, cum_hat_obs[,-1]-cum_hat_obs[,-ncol(cum_hat_obs)])
  
  p1_integrand <- sum_integrand*S_hat0_obs*S_hat1_obs*dH
  p1_integrand_pruned <- p1_integrand[,max_time_filter]
  p1 <- rowSums(p1_integrand_pruned)
  p1_alt <- rowSums((pi_a_1*(at_risk-S_hat1_obs)*S_hat0_obs*dH)[,max_time_filter])
  p2_alt <- rowSums((pi_a_0*(dN/alpha_A_0-S_hat0_obs)*S_hat1_obs*dH)[,max_time_filter])
  
  #p2 
  p_extend_surv <- P_treatment_extend_survival(model = model_surv, max_time = max_time, model_cov = covar_true)
  p3 <- p_extend_surv$res[,1]
  dim(p1_alt)

  
  return(c(mean(p1_alt), mean(p2_alt), mean(p3), mean(p1_alt+p2_alt + p3)))
}

res <- test()
res
