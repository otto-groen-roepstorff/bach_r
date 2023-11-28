

test <- function(){
  
n <- 1000
#test_data <- generate_survival_data(n, ba_t = 1, bx_t = 1, bz_t = log(6), surv_is_cox = F,
#                                    ba_c = log(2), bx_c = 1, bz_c = -1, cens_is_cox = F)
train_data <- generate_survival_data(n, ba_t = -1, bx_t = log(2), bz_t = log(2), surv_is_cox = T,
                                     ba_c = 1, bx_c = log(2), bz_c = 1, cens_is_cox = T)

max_time <- 2

################################################
#Martingales
################################################
martingale_estimates <- estimate_martingales(train_data = train_data)
jump_times <- martingale_estimates$jump_times
jump_times_to_keep <- sum(jump_times <= max_time)

dM <- martingale_estimates$mod_dM[,1:jump_times_to_keep]
dL <- martingale_estimates$mod_dL[,1:jump_times_to_keep]




#Survival function estimates
S_hats <- get_S_hats(train_data)
S_hat <- S_hats$Shat[,1:jump_times_to_keep]
S_hat0 <- S_hats$Shat_0[,1:jump_times_to_keep]
S_hat1 <- S_hats$Shat_1[,1:jump_times_to_keep]
betahat <- S_hats$beta_hat


#Censoring function estimate
K_C_hat <- get_K_hat(train_data, surv_is_cox = T)[,1:jump_times_to_keep]


#Propensity scores
prop <- propensity(train_data)$propens
prop_0 <- prop[,'pi_a_0']
prop_1 <- prop[,'pi_a_1']
prop_sum <- prop_0 + prop_1



#Part 1
p1 <- prop_0 * rowSums(S_hat1/K_C_hat*dM)


#Part 2
p_treat_list <- P_treatment_extend_survival(train_data, max_time = max_time)
p2 <- p_treat_list$res

#Part 3
inner_integrand <- dM/(S_hat * K_C_hat)
inner_integral <- t(apply(inner_integrand, MARGIN = 1, FUN = cumsum))


p3 <- prop_sum * p_treat_list$multiplier * rowSums(p_treat_list$integrand[,1:jump_times_to_keep] * inner_integral)


print('Run done')

return(mean(p1+p2-p3))
}

res_T_T <- replicate(10, test())
mean(res_T_T)
beta_hat <- c(-1,log(2))
theoretical_value_11(beta_hat, 2)


test_2 <- function(){
  n <- 1000
  train_data <- generate_survival_data(n, ba_t = -1, bx_t = log(2), bz_t = log(2), surv_is_cox = T,
                                       ba_c = -1, bx_c = log(2), bz_c = 1, cens_is_cox = T)
  
  p_treat_list <- P_treatment_extend_survival(train_data)
  p2 <- p_treat_list$res
}

test_2_reps <- replicate(10, test_2())











res_T_F <- replicate(250, test())
res_F_T <- replicate(250, test())
res_F_F <- replicate(250, test())




integrand <- t(apply((dM/(S_hat*K_C_hat)), MARGIN = 1, FUN = cumsum))

p3 <- (pi_a_0 + pi_a_1) * mean(exp(betahat[2] * test_data$X) * rowSums(integrand*S_hat0*S_hat1*dL))
p3_est[i,] <- p3


estimates[i,] <- p1 + p2 + p3




par(mfrow = c(1,1))
hist(estimates)
par(mfrow = c(1,3))
hist(p1_est)
hist(p2_est)
hist(p3_est)

