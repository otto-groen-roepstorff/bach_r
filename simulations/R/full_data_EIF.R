

test <- function(){
  
  n <- 1000
  data <- generate_survival_data(n, ba_t = -1, bx_t = log(2), bz_t = log(2),
                                 ba_c = 1, bx_c = log(2), bz_c = -2, seed = sample(1:1000))
  #Changing all data to observed data
  data$Uncensored <- T
  data$T_obs <- data$T_true
  
  #creating model
  model_surv <- oracle_model(data)
  
  #Calculating propensities
  props <- propensity(data)
  pi_a_0 <- props$propens$pi_a_0 #2*(1-get_A_val(data))#props$propens$pi_a_0
  pi_a_1 <- props$propens$pi_a_1 #2*get_A_val(data) #props$propens$pi_a_1
  pi_sum = pi_a_0 + pi_a_1
  
  #Getting model estimates
  cum_matrix_obs <- get_cum(model_surv)
  jump_times_obs <- get_jump_times(cum_matrix_obs)
  cum_bas_haz_obs <- get_cum_base_haz(cum_matrix_obs)
  no_jumps_obs <- get_row_length(cum_matrix_obs)
  beta_hat_obs <- get_param_est(model_surv)
  max_time <- 100
  max_time_filter <- jump_times_obs<= max_time
  
  covar_true <- get_oracle_covar(data)
  #Counterfactual data
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
  
  
  #Finding the different components
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
  dt <- cbind(0, jump_times_obs[-1] - jump_times_obs[-length(jump_times_obs)])
  
  p1_integrand <- sum_integrand*S_hat0_obs*S_hat1_obs*dH
  p1_integrand_pruned <- p1_integrand[,max_time_filter]
  p1 <- rowSums(p1_integrand_pruned)
  
  p1_alt <- rowSums((pi_a_1*(at_risk-S_hat1_obs)*S_hat0_obs*dH) %>% prune_data(max_time_filter = max_time_filter))
  
  
  p2_alt <- rowSums((pi_a_0*(dN-alpha_A_0*S_hat0_obs)*S_hat1_obs*dt) %>% prune_data(max_time_filter = max_time_filter))
  testies <- pi_a_0*(dN/(alpha_A_0*S_hat0_obs)-1)
  View(testies)
  
  #p2 
  p_extend_surv <- P_treatment_extend_survival(model = model_surv, max_time = max_time, model_cov = covar_true)
  integrandzz <- p_extend_surv$res
  
  
  p3 <- rowSums(integrandzz %>% prune_data(max_time_filter = max_time_filter)) 
  
  
  return(c(mean(p1_alt), mean(p2_alt), mean(p3), mean(p1_alt+p2_alt + p3)))
}

test2 <- function(n = 300, ba_t = -1, max_time = 3){
  
  data <- generate_survival_data(n, ba_t = ba_t, bx_t = log(2), bz_t = log(2),
                                 ba_c = 1, bx_c = log(2), bz_c = -2, seed = sample(1:100000))
  #Changing all data to observed data
  data$Uncensored <- T
  data$T_obs <- data$T_true
  T_obs <- data$T_obs
  
  #creating model
  model_surv <- non_oracle_nice_model(data)
  covar_true <- get_n_oracle_nice_covar(data)
  
  #Calculating propensities
  props <- propensity(data)
  pi_a_0 <- props$propens$pi_a_0 #2*(1-get_A_val(data))#props$propens$pi_a_0
  pi_a_1 <- props$propens$pi_a_1 #2*get_A_val(data) #props$propens$pi_a_1
  pi_sum = pi_a_0 + pi_a_1
  
  #Getting model estimates
  cum_matrix_obs <- get_cum(model_surv)
  jump_times_obs <- get_jump_times(cum_matrix_obs)
  cum_bas_haz_obs <- get_cum_base_haz(cum_matrix_obs)
  no_jumps_obs <- get_row_length(cum_matrix_obs)
  beta_hat_obs <- get_param_est(model_surv)
  max_time_filter <- jump_times_obs<= max_time
  #Counterfactual data
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
  pruned_S_hat0_obs <- S_hat0_obs %>% prune_data(max_time_filter = max_time_filter)
  
  at_risk <- outer(X = T_obs, Y = jump_times_obs, FUN = ">=")
  dN <- outer(X = T_obs, Y = jump_times_obs, FUN = "==")
  sup_observed_time <- max(jump_times_obs[jump_times_obs<= max_time])
  
  
  
  
  tau_prune <-  matrix(rep(sup_observed_time == jump_times_obs, n), byrow = T, nrow = n, ncol = no_jumps_obs)
  
  obs_filter <- T_obs<= max_time
  chosen_time <- dN*obs_filter+tau_prune*(1-obs_filter)
  
  p_extend_surv <- P_treatment_extend_survival(model = model_surv, max_time = max_time, model_cov = covar_true)
  phi_W <- p_extend_surv$res[,1]
  
  
  
  observed_S_hat0 <- S_hat0_obs*chosen_time
  observed_S_hat1 <- data.matrix(S_hat1_obs*dN)[, max_time_filter]
  rowSums(observed_S_hat0)
  
  p1 <- pi_a_1*(1-rowSums(observed_S_hat0))
  p2 <- pi_a_0*rowSums(observed_S_hat1)
  p3 <- (pi_a_0+pi_a_1-1)*phi_W
    
  p_estimate <- p1+ p2- p3
    
  output <- c(de_biased_estimate = mean(p_estimate), de_bias = mean(p_estimate-phi_W), naive_estimate = mean(phi_W))#p1 = p1, p2 = p2, p3 = p3, )
  
  return(output)
}

start_time <- Sys.time()
res <- replicate(250, test2(n =2400))
end_time <- Sys.time()
round(end_time - start_time, 2)
rowMeans(res)

#250 sim
#n = 300,   33  s
#n = 600,    1  min
#n = 1200,   4.7min
#n = 2400,  11  min 

#True value 0.721570837

#Oracle model - 250 sim
#de_biased_est  0.72792141 - 300, 0.727393269 - 600, 0.722685203 - 1200, 0.721944585 - 2400                
#naive_es       0.71754802 - 300, 0.720787283 - 600, 0.719998011 - 1200, 0.720906869 - 2400

#Non-oracle nice - 250 sim
#de_biased_est  0.71433986 - 300, 0.71708553 - 600,  0.711820790 - 1200,  - 2400                         
#naive_es       0.69848042 - 300, 0.70675960 - 600,  0.703265426 - 1200,  - 2400



#non


