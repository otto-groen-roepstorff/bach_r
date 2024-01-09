

test <- function(){
  
n <- 1000
#test_data <- generate_survival_data(n, ba_t = 1, bx_t = 1, bz_t = log(6), surv_is_cox = F,
#                                    ba_c = log(2), bx_c = 1, bz_c = -1, cens_is_cox = F)
train_data <- generate_survival_data(n, ba_t = 0, bx_t = log(6), bz_t = log(2),
                                     ba_c = 1, bx_c = log(2), bz_c = 1)

max_time <- 100

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

res <- list('p1_mean' = mean(p1), 'p2_mean' = mean(p2), 'p3_mean' = mean(p3), 'res' = mean(p1+p2-p3))

return(mean(p2))
}

res_T_T <- replicate(50, test())
mean(res_T_T)
beta_hat <- c(-1,log(2))
theoretical_value_11(beta_hat, 2)

correct_res_T_T

between(c(0.7070660870-1.96*sd(res_T_T), 0.7070660870+1.96*sd(res_T_T))


