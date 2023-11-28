library(tidyverse)

n <- 1000
#test_data <- generate_survival_data(n, ba_t = 1, bx_t = 1, bz_t = log(6), surv_is_cox = F,
#                                    ba_c = log(2), bx_c = 1, bz_c = -1, cens_is_cox = F)

train_data <- generate_survival_data(n, ba_t = 0, bx_t = 1, bz_t = log(2), surv_is_cox = F,
                                     ba_c = log(2), bx_c = 1, bz_c = -1, cens_is_cox = F)


df1 <- train_data[,!names(train_data) %in% c("X", "Z", "Uncensored")]
df2 <- pivot_longer(data = df1, cols = c("T_true", "C", "T_obs"), names_to = "Source", )


#x <- rgamma(n/2, shape = 5, rate = 1)
#y <- rgamma(n/2, shape = 25, rate = 2)
#z <- c(x,y)
#a_z <- c(rep(0, n/2), rep(1, n/2))
#source_z <- c(rep("x", n/2), rep("y", n/2))
#z_df <- data.frame(A= a_z, Source = source_z, value = z)
full_df <- df2 #bind_rows(z_df, df2)
full_df$A <- as.factor(full_df$A)

ggplot(full_df)+ geom_density(aes(value, fill = A), alpha = 0.4) + facet_wrap(~Source)
proportion_observed(train_data)

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