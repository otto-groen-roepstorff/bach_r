

p1_est <- matrix(data = NA, nrow = 100, ncol = 1)
p2_est <- matrix(data = NA, nrow = 100, ncol = 1)
p3_est <- matrix(data = NA, nrow = 100, ncol = 1)
estimates <- matrix(data = NA, nrow = 100, ncol = 1)

test <- function(){
  
n <- 1000
#test_data <- generate_survival_data(n, ba_t = 1, bx_t = 1, bz_t = log(6), surv_is_cox = F,
#                                    ba_c = log(2), bx_c = 1, bz_c = -1, cens_is_cox = F)

train_data <- generate_survival_data(n, ba_t = 0, bx_t = 1, bz_t = log(2), surv_is_cox = T,
                                     ba_c = log(2), bx_c = 1, bz_c = -1, cens_is_cox = F)

#Survival function estimates
S_hats <- get_S_hats(train_data)
S_hat <- S_hats$Shat
S_hat0 <- S_hats$Shat_0
S_hat1 <- S_hats$Shat_1
betahat <- S_hats$beta_hat


#Censoring function estimate
K_C_hat <- get_K_hat(train_data, surv_is_cox = T)

################################################
#Martingales
################################################
martingale_estimates <- estimate_martingales(test_data = test_data, train_data = train_data)
jump_times <- martingale_estimates$jump_times
dM <- martingale_estimates$mod_dM
dL <- martingale_estimates$mod_dL


#Propensity scores
prop <- propensity(train_data)$propens
prop_0 <- prop[,'pi_a_0']
prop_1 <- prop[,'pi_a_1']
prop_sum <- prop_0 + prop_1



#Part 1
p1 <- prop_0 * rowSums(S_hat1/K_C_hat*dM)


#Part 2
p2 <- P_treatment_extend_survival(train_data)

#Part 3
inner_integrand <- dM/(S_hat * K_C_hat)
inner_integral <- t(apply(inner_integrand, MARGIN = 1, FUN = cumsum))


tS <- cbind(0,(S_hat0[,-1] - S_hat0[,-ncol(S_hat0)]))
#tt <- jump_times[-1] - jump_times[-length(jump_times)]
#S_hat0_deriv <- cbind(0, t(t(tS)/tt))

outer_integrand <- inner_integral * tS * S_hat1
outer_integral <- rowSums(outer_integrand)

p3 <- prop_sum*outer_integral

return(mean(p1+p2+p3))
}

res <- replicate(100, test())


integrand <- t(apply((dM/(S_hat*K_C_hat)), MARGIN = 1, FUN = cumsum))

p3 <- (pi_a_0 + pi_a_1) * mean(exp(betahat[2] * test_data$X) * rowSums(integrand*S_hat0*S_hat1*dL))
p3_est[i,] <- p3


estimates[i,] <- p1 + p2 + p3



}
par(mfrow = c(1,1))
hist(estimates)
par(mfrow = c(1,3))
hist(p1_est)
hist(p2_est)
hist(p3_est)

