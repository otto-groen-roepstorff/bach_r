


n <- 500
data <- generate_survival_data(n, ba_t = 0, bx_t = log(2), bz_t = log(2),
                                     ba_c = 1, bx_c = log(2), bz_c = 1)
EIF_test <- EIF(data, T_corr = F, max_time = 2)
EIF_test[4]

n_sims <- 250
simulation <- matrix(data = NA, ncol = 4, nrow = n_sims)
max_time <- 2

for (i in 1:n_sims){
  n <- 300
  data <- generate_survival_data(n, ba_t = 0, bx_t = log(2), bz_t = log(2),
                                 ba_c = 1, bx_c = log(2), bz_c = 1, seed = i)
  
  EIF_sim <- EIF(data, Cens_corr = T, T_corr = F, max_time = max_time)
  simulation[i,1] <- mean(EIF_sim[[1]])
  simulation[i,2] <- mean(EIF_sim[[2]])
  simulation[i,3] <- mean(EIF_sim[[3]])
  simulation[i,4] <- EIF_sim[[4]]
  print(paste0('Simulation ',i,' done'))
}
mean(simulation[,1])




test <- function(){

n <- 1000
data <- generate_survival_data(n, ba_t = 0, bx_t = log(2), bz_t = log(2),
                               ba_c = 1, bx_c = log(2), bz_c = -2, seed = sample(1:1000))


model_surv <- non_oracle_model(data)
model_cens <- oracle_cens_model(data)


props <- propensity(data)
pi_a_0 <- 2*(1-get_A_val(data))#props$propens$pi_a_0

pi_a_1 <- 2*get_A_val(data) #props$propens$pi_a_1
pi_sum = pi_a_0 + pi_a_1


cum_matrix_obs <- get_cum(model_surv)
jump_times_obs <- get_jump_times(cum_matrix_obs)
cum_bas_haz_obs <- get_cum_base_haz(cum_matrix_obs)
no_jumps_obs <- get_row_length(cum_matrix_obs)
beta_hat_obs <- get_param_est(model_surv)
max_time <- 3
max_time_filter <- jump_times_obs<max_time



covar_true <- get_n_oracle_covar(data)
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

covar_cens <- get_oracle_cens_covar(data)
beta_hat_cens <- get_param_est(model_cens)

cum_hat_cens <- predict_cox.aalen(covar = covar_cens, betaHat = beta_hat_cens, cum_base_haz = cum_bas_haz_cens)

S_hat_cens <- exp(- cum_hat_cens)


Khat_cens = matrix(NA, nrow = n, ncol = no_jumps_obs) #Maybe this can also be vectorized??
for(i in 1:no_jumps_obs){
  Khat_cens[,i] = S_hat_cens[,max((1:no_jumps_cens)[jump_times_cens<=jump_times_obs[i]])]
}

#Finding martingales
T_obs <- get_observed_times(data)
at_risk <- outer(X = T_obs, Y = jump_times_obs, FUN = ">=")
#Cumulative hazard Y*alpha
dN <- outer(X = T_obs, Y = jump_times_obs, FUN = "==")
dH <- cbind(0, cum_hat_obs[,-1]-cum_hat_obs[,-ncol(cum_hat_obs)])
dL <- at_risk*dH
dM <- dN - dL


p1_integrand <- (S_hat1_obs/Khat_cens)*dM
p1_integrand_pruned <- p1_integrand[,max_time_filter]
p1 <- pi_a_0*rowSums(p1_integrand_pruned)
p1_fake <- pi_a_0*rowSums(p1_integrand)

#p2 
p_extend_surv <- P_treatment_extend_survival(model = model_surv, max_time = max_time, model_cov = covar_true)
p2 <- p_extend_surv$res

#p3
scalar <- p_extend_surv$multiplier

sim_M <- (dM/(S_hat_obs*Khat_cens))[, max_time_filter]
cum_sim_M <- t(apply(sim_M, MARGIN = 1, cumsum))

p3_integral <- rowSums((p_extend_surv$integrand)*cum_sim_M)
p3 <- -pi_sum*scalar*p3_integral

return(c(mean(p1_fake), mean(p1), mean(p2), mean(p3), mean(p1+p2+p3)))
}


res <- replicate(100, test())
rowMeans(res)




res



