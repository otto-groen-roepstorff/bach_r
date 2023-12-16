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
  
  
  props <- propensity(data, glm = F)
  
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
  
  EIF_var_est <- var(p1 + p2 + p3 - p4)/n
  
  
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


#Variance:  Correctly specified:        0.0005908318
#           Misspecified propensity:    0.0007062777

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




#Simulations for ATE

n_sims <- 1000
simulation_list <- list()
ba_t <- -1
tau <- 3
n <- 600

start.time <- Sys.time()
for (j in 1:n_sims){
  data <- generate_survival_data(n, ba_t = ba_t, bx_t = log(2), bz_t = log(2),
                                 ba_c = 1, bx_c = log(2), bz_c = log(2), prop_a = 1/2, seed = sample(1:100000))
  
  #data_for_breslow <- generate_survival_data(n, ba_t = ba_t, bx_t = log(2), bz_t = log(2),
  #                                           ba_c = 1, bx_c = log(2), bz_c = log(2), seed = sample(1:1000000))
  
  
  model_surv <- non_oracle_model(data)
  model_cens <- oracle_cens_model(data)
  
  covar_true <- get_n_oracle_covar(data)
  covar_cens <- get_oracle_cens_covar(data)
  
  
  props <- propensity(data, glm = T)
  
  
  cum_matrix_obs <- get_cum(model_surv)
  jump_times_obs <- get_jump_times(cum_matrix_obs)
  cum_bas_haz_obs <- get_cum_base_haz(cum_matrix_obs)
  no_jumps_obs <- get_row_length(cum_matrix_obs)
  beta_hat_obs <- get_param_est(model_surv)
  
  
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
  
  plot(jump_times_obs, colMeans(S_hat_obs), type = 'l')
  lines(jump_times_obs, colMeans(S_hat0_obs), col = 'red')
  lines(jump_times_obs, colMeans(S_hat1_obs), col = 'blue')
  
  
  #Censoring
  cum_matrix_cens <- get_cum(model_cens)
  jump_times_cens <- get_jump_times(cum_matrix_cens)
  cum_bas_haz_cens <- get_cum_base_haz(cum_matrix_cens)
  no_jumps_cens <- get_row_length(cum_matrix_cens)
  beta_hat_cens <- get_param_est(model_cens)
  
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
  
  #par(mfrow = c(1,1))
  #plot(jump_times_obs, cumsum(colSums(dM)), type = 'l', ylab = 'Kummuleret dM.')
  #plot(jump_times_obs, cumsum(colSums(dN)), type = 'l', ylab = 'Kummuleret dN., dL.', col = 'red')
  #lines(jump_times_obs, cumsum(colSums(at_risk*dL)), col = alpha('blue', 0.2))
  
  #plot(jump_times_obs, cumsum(dM[233,]), type = 'l', ylim = c(-3,1.1))
  #lines(jump_times_obs, cumsum(dN[233,]))
  #lines(jump_times_obs, -cumsum(at_risk[233,]*dL[233,]))
  
  #EIF
  no_times_to_keep <- jump_times_obs <= tau
  jump_times_kept <- jump_times_obs[no_times_to_keep]
  
  d_cum_base_haz <- c(0, cum_bas_haz_obs[-1]-cum_bas_haz_obs[-no_jumps_obs])
  d_cum_base_haz_tau <- c(0, cum_bas_haz_obs[-1]-cum_bas_haz_obs[-no_jumps_obs])[no_times_to_keep]
  
  S_hat_obs_tau <- S_hat_obs[,no_times_to_keep]
  S_hat0_obs_tau <- S_hat0_obs[,no_times_to_keep]
  S_hat1_obs_tau <- S_hat1_obs[,no_times_to_keep]
  
  Khat_cens_tau <- Khat_cens[,no_times_to_keep]
  
  dM_tau <- dM[,no_times_to_keep]
  dN_tau <- dN[,no_times_to_keep]
  dLambda_tau <- (at_risk * dL)[,no_times_to_keep]
  
  pi_0 <- props$propens$pi_a_0
  pi_1 <- props$propens$pi_a_1
  
  
  #One-step estimator
  
  #Inner integrand
  inner_integrand <- dM_tau/(S_hat_obs_tau*Khat_cens_tau)
  inner_integral <- t(apply(inner_integrand, MARGIN = 1, FUN = cumsum))
  
  #Part 1
  p1 <- pi_1 * inner_integral * S_hat_obs_tau
  
  #Part 2
  p2 <- S_hat1_obs_tau
  
  
  # The Whole Thing
  TWT <- p1 + p2
  
  #The estimator
  one_step_estimator <- colMeans(TWT)
  
  
  new_element <- list(jump_times = jump_times_kept, one_step_estimator = one_step_estimator)
  simulation_list[[length(simulation_list) + 1]] <- new_element
  print(paste0('Simulation ',j,' done'))
}
end.time <- Sys.time()
end.time - start.time


plot(simulation_list[[1]][[1]], 
     simulation_list[[1]][[2]], 
     xlab = 'Time', 
     ylab = 'Estimate', 
     type = 'l', 
     col = rgb(red = 0.5, green = 0, blue = 0.7, alpha = 0.3))
for (i in 24:24){
  lines(simulation_list[[i]][[1]], simulation_list[[i]][[2]],
        col = rgb(red = 0.5, green = 0, blue = 0.7, alpha = 0.3))
}
lines(truth$taus, truth$true_vals)

#True values bx = bz = log(2), ba = -1, tau = 3
taus <- seq(0,3,0.1)
true_values <- c(1, 0.9445018885, 
                 0.8926635295, 0.8442070210, 
                 0.7988779580, 0.7564432805, 
                 0.7166893300, 0.6794200945, 
                 0.6444556235, 0.6116305930, 
                 0.5807930095, 0.5518030325, 
                 0.5245319145, 0.4988610342, 
                 0.4746810230, 0.4518909722, 
                 0.4303977123, 0.4101151594, 
                 0.3909637198, 0.3728697502, 
                 0.3557650644, 0.3395864864, 
                 0.3242754414, 0.3097775838, 
                 0.2960424576, 0.2830231864, 
                 0.2706761908, 0.2589609291, 
                 0.2478396611, 0.2372772317, 
                 0.2272408724)
truth <- data.frame('taus' = taus, 'true_vals' = true_values)
lines(truth$taus, truth$true_vals)














################################################################
#                         COXPH
################################################################


################################################################
#                     Generating data
################################################################
n_sims <- 1000
ba_t <- -1
tau <- 3
n <- 300
OSEC_lst      <- list()
OSEM_lst      <- list()
OSEBM_lst     <- list()
NAIVE_lst     <- list()
NAIVE_MIS_lst <- list()

for (j in 1:n_sims){
  data <- generate_survival_data(n, ba_t = ba_t, bx_t = log(2), bz_t = log(2),
                                 ba_c = 1, bx_c = log(2), bz_c = log(2), seed = sample(1:100000))
  
  data_A0 <- data %>% mutate(A = 0)
  data_A1 <- data %>% mutate(A = 1)
  
  
  ################################################################
  #                     Fitting models
  ################################################################
  fit_coxph <- coxph(Surv(T_obs, Uncensored) ~ A + X + Z, data = data)
  fit_coxph_mis <- coxph(Surv(T_obs, Uncensored) ~ A + fake_1 + fake_2 + fake_3, data = data)
  
  cum_base_haz <- basehaz(fit_coxph, newdata = data)
  cum_base_haz_mis <- basehaz(fit_coxph_mis, newdata = data)
  
  survfit_all <- survfit(fit_coxph, newdata = data)
  survfit_all_mis <- survfit(fit_coxph_mis, newdata = data)
  
  #plot(survfit_all$time, cumsum(fit_coxph$residuals), type = 'l')
  #lines(survfit_all_mis$time, cumsum(fit_coxph_mis$residuals), type = 'l', col = 'blue')
  
  
  
  ################################################################
  #           Survival functions and cumulative hazards
  ################################################################
  
  #----------------Correctly specified--------------
  survfit <- survfit(fit_coxph, newdata = data)
  survfit_0 <- survfit(fit_coxph, newdata = data_A0)
  survfit_1 <- survfit(fit_coxph, newdata = data_A1)
  
  cum_haz <- survfit$cumhaz
  cum_haz_0 <- survfit_0$cumhaz
  cum_haz_1 <- survfit_1$cumhaz
  
  surv <- survfit$surv
  surv_0 <- survfit_0$surv
  surv_1 <- survfit_1$surv
  
  
  #plot(survfit_0$time, surv_0[,1], type = 'l', col = 'red')
  #for (i in 2:length(survfit_0$time)){
  #  lines(survfit_0$time, surv_0[,i], col = 'red')
  #  lines(survfit_1$time, surv_1[,i], col = 'blue')
  #}
  
  
  #----------------Misspecified--------------
  survfit_mis <- survfit(fit_coxph_mis, newdata = data)
  survfit_0_mis <- survfit(fit_coxph_mis, newdata = data_A0)
  survfit_1_mis <- survfit(fit_coxph_mis, newdata = data_A1)
  
  cum_haz_mis <- survfit_mis$cumhaz
  cum_haz_0_mis <- survfit_0_mis$cumhaz
  cum_haz_1_mis <- survfit_1_mis$cumhaz
  
  surv_mis <- survfit_mis$surv
  surv_0_mis <- survfit_0_mis$surv
  surv_1_mis <- survfit_1_mis$surv
  
  #plot(survfit_0$time, surv_0_mis[,1], type = 'l', col = 'red')
  #for (i in 2:length(survfit_0$time)){
  #  lines(survfit_0$time, surv_0_mis[,i], col = 'red')
  #  lines(survfit_1$time, surv_1_mis[,i], col = 'blue')
  #}
  
  
  ################################################################
  #                     Martingale errors
  ################################################################
  
  #Note these martingales are aggregated across all individuals
  dM <- fit_coxph$residuals
  dM_mis <- fit_coxph_mis$residuals
  
  
  
  #----------------Manual calculations--------------
  cum_haz <- t(survfit_all$cumhaz)
  cum_haz_mis <- t(survfit_all_mis$cumhaz)
  
  dL <- cbind(0, cum_haz[,-1] - cum_haz[,-ncol(cum_haz)])
  dL_mis <- cbind(0, cum_haz_mis[,-1] - cum_haz_mis[,-ncol(cum_haz)])
  
  dN <- outer(data$T_obs, survfit_all$time, '==') * data$Uncensored
  at_risk <- outer(data$T_obs, survfit_all$time, '>=')
  
  dM_manual <- dN - dL*at_risk
  dM_mis_manual <- dN - dL_mis*at_risk
  
  #plot(survfit_all$time, cumsum(rowSums(dN)), type = 'l')
  #lines(survfit_all$time, cumsum(rowSums(at_risk*dL)))
  #lines(survfit_all$time, cumsum(rowSums(at_risk*dL_mis)))
  
  
  
  ################################################################
  #                     Censoring distribution
  ################################################################
  fit_coxph_cens <- coxph(Surv(T_obs, Uncensored == F) ~ A + X + Z, data = data)
  censfit <- survfit(fit_coxph_cens, newdata = data)
  cens <- censfit$surv
  
  fit_coxph_cens_mis <- coxph(Surv(T_obs, Uncensored == F) ~ A + fake_1 + fake_2 + fake_3, data = data)
  censfit_mis <- survfit(fit_coxph_cens_mis, newdata = data)
  cens_mis <- censfit_mis$surv
  
  
  
  
  ################################################################
  #                           EIF
  ################################################################
  
  time_indicator <- survfit_all$time <= tau
  jump_times_kept <- survfit_all$time[time_indicator]
  
  dL_tau <- dL[time_indicator]
  dL_mis_tau <- dL_mis[time_indicator]
  
  surv_tau    <-  t(surv)[,time_indicator]
  surv_0_tau  <-  t(surv_0)[,time_indicator] 
  surv_1_tau  <-  t(surv_1)[,time_indicator]
  
  surv_mis_tau    <-  t(surv_mis)[,time_indicator]
  surv_0_mis_tau  <-  t(surv_0_mis)[,time_indicator] 
  surv_1_mis_tau  <-  t(surv_1_mis)[,time_indicator]
  
  
  Khat_cens_tau <- t(cens)[,time_indicator]
  Khat_cens_mis_tau <- t(cens_mis)[,time_indicator]
  
  dM_manual_tau       <- dM_manual[, time_indicator]
  dM_mis_manual_tau   <- dM_mis_manual[, time_indicator]
  
  
  props <- propensity(data, glm = T)
  pi_0 <- props$propens$pi_a_0
  pi_1 <- props$propens$pi_a_1
  
  
  #One-step estimator
  
  #Correctly specified
  inner_integrand <- dM_manual_tau/(surv_tau*Khat_cens_tau)
  inner_integral <- t(apply(inner_integrand, MARGIN = 1, FUN = cumsum))
  
  OSEC <- colMeans(pi_1 * inner_integral * surv_tau + surv_1_tau)
  
  #Misspecified
  inner_integrand_mis <- dM_mis_manual_tau/(surv_mis_tau*Khat_cens_tau)
  inner_integral_mis <- t(apply(inner_integrand_mis, MARGIN = 1, FUN = cumsum))
  
  OSEM <- colMeans(pi_1 * inner_integral_mis * surv_mis_tau + surv_1_mis_tau)
  
  #Both misspecified
  inner_integrand_mis <- dM_mis_manual_tau/(surv_mis_tau*Khat_cens_mis_tau)
  inner_integral_mis <- t(apply(inner_integrand_mis, MARGIN = 1, FUN = cumsum))
  
  OSEBM <- colMeans(pi_1 * inner_integral_mis * surv_mis_tau + surv_1_mis_tau)
  
  
  NAIVE <- colMeans(surv_1_tau)
  NAIVE_MIS <- colMeans(surv_1_mis_tau)
  
  
  OSEC_new_element          <- list(jump_times = jump_times_kept, one_step_estimator = OSEC)
  OSEM_new_element          <- list(jump_times = jump_times_kept, one_step_estimator = OSEM)
  OSEBM_new_element         <- list(jump_times = jump_times_kept, one_step_estimator = OSEBM)
  NAIVE_new_element         <- list(jump_times = jump_times_kept, one_step_estimator = NAIVE)
  NAIVE_MIS_new_element     <- list(jump_times = jump_times_kept, one_step_estimator = NAIVE_MIS)
  
  
  OSEC_lst[[length(OSEC_lst) + 1]] <- OSEC_new_element
  OSEM_lst[[length(OSEM_lst) + 1]] <- OSEM_new_element
  OSEBM_lst[[length(OSEBM_lst) + 1]] <- OSEBM_new_element
  NAIVE_lst[[length(NAIVE_lst) + 1]] <- NAIVE_new_element
  NAIVE_MIS_lst[[length(NAIVE_MIS_lst) + 1]] <- NAIVE_MIS_new_element
  
  
  print(paste0('Simulation ',j,' done'))
}

#True values bx = bz = log(2), ba = -1, tau = 3
taus <- seq(0,3,0.1)
true_values <- c(1, 0.9445018885, 
                 0.8926635295, 0.8442070210, 
                 0.7988779580, 0.7564432805, 
                 0.7166893300, 0.6794200945, 
                 0.6444556235, 0.6116305930, 
                 0.5807930095, 0.5518030325, 
                 0.5245319145, 0.4988610342, 
                 0.4746810230, 0.4518909722, 
                 0.4303977123, 0.4101151594, 
                 0.3909637198, 0.3728697502, 
                 0.3557650644, 0.3395864864, 
                 0.3242754414, 0.3097775838, 
                 0.2960424576, 0.2830231864, 
                 0.2706761908, 0.2589609291, 
                 0.2478396611, 0.2372772317, 
                 0.2272408724)
truth <- data.frame('taus' = taus, 'true_vals' = true_values)


plot(jump_times_kept, OSEC, type = 'l', col = 'blue')
lines(jump_times_kept, OSEM, type = 'l', col = 'red')
lines(jump_times_kept, NAIVE, type = 'l', col = 'green')
lines(jump_times_kept, NAIVE_MIS, type = 'l', col = 'yellow')
lines(truth$taus, truth$true_vals)

plot(OSEC_lst[[1]][[1]], 
     OSEC_lst[[1]][[2]], 
     main = 'Correctly specified',
     xlab = 'Time', 
     ylab = 'Estimate', 
     type = 'l', 
     col = rgb(red = 0.5, green = 0, blue = 0.7, alpha = 0.3))
for (i in 2:length(OSEC_lst)){
  lines(OSEC_lst[[i]][[1]], OSEC_lst[[i]][[2]],
        col = rgb(red = 0.5, green = 0, blue = 0.7, alpha = 0.3))
}
lines(truth$taus, truth$true_vals)


plot(OSEM_lst[[1]][[1]], 
     OSEM_lst[[1]][[2]], 
     main = 'Misspecified survival',
     xlab = 'Time', 
     ylab = 'Estimate', 
     type = 'l', 
     col = rgb(red = 0.5, green = 0, blue = 0.7, alpha = 0.3),
     ylim = c(0,1))
for (i in 2:length(OSEM_lst)){
  lines(OSEM_lst[[i]][[1]], OSEM_lst[[i]][[2]],
        col = rgb(red = 0.5, green = 0, blue = 0.7, alpha = 0.3))
}
lines(truth$taus, truth$true_vals)

plot(OSEBM_lst[[1]][[1]], 
     OSEBM_lst[[1]][[2]], 
     main = "Both censoring and survival misspecified",
     xlab = 'Time', 
     ylab = 'Estimate', 
     type = 'l', 
     col = rgb(red = 0.5, green = 0, blue = 0.7, alpha = 0.3),
     ylim = c(0,1))
for (i in 2:length(OSEBM_lst)){
  lines(OSEBM_lst[[i]][[1]], OSEBM_lst[[i]][[2]],
        col = rgb(red = 0.5, green = 0, blue = 0.7, alpha = 0.3))
}
lines(truth$taus, truth$true_vals)


plot(NAIVE_lst[[1]][[1]], 
     NAIVE_lst[[1]][[2]], 
     main = "Naive estimate",
     xlab = 'Time', 
     ylab = 'Estimate', 
     type = 'l', 
     col = rgb(red = 0.5, green = 0, blue = 0.7, alpha = 0.3),
     ylim = c(0,1))
for (i in 2:length(NAIVE_lst)){
  lines(NAIVE_lst[[i]][[1]], NAIVE_lst[[i]][[2]],
        col = rgb(red = 0.5, green = 0, blue = 0.7, alpha = 0.3))
}
lines(truth$taus, truth$true_vals)


plot(NAIVE_MIS_lst[[1]][[1]], 
     NAIVE_MIS_lst[[1]][[2]], 
     main = "Naive misspecified estimate",
     xlab = 'Time', 
     ylab = 'Estimate', 
     type = 'l', 
     col = rgb(red = 0.5, green = 0, blue = 0.7, alpha = 0.3),
     ylim = c(0,1))
for (i in 2:length(NAIVE_MIS_lst)){
  lines(NAIVE_MIS_lst[[i]][[1]], NAIVE_MIS_lst[[i]][[2]],
        col = rgb(red = 0.5, green = 0, blue = 0.7, alpha = 0.3))
}
lines(truth$taus, truth$true_vals)




