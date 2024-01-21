library(glmnet)
library(rje)
library(gplm)
library(np)
library(targeted)
library(ibr)





##################################################
#             Double robustness study            #
##################################################

n_sims <- 1000000
sim_matrix <- matrix(data = NA, ncol = 1, nrow = n_sims)

for (j in 1:n_sims){
  X <- runif(n,-2,2)
  Z <- rbinom(n, 1, 0.5)
  A <- rbinom(n, 1, 0.5)
  Y <- expit(1 + 3*X + 3*Z)
  
  sim_matrix[j,1] <- mean(Y)
  
  if (j %% 100000 == 0){
    print(paste0('Sim ',j,' done')) 
  }
}
mean(sim_matrix)
#true mean 0.6629087 


study <- function(n){
  
  
  n_sims <- 2000
  result_vector <- rep(NA, n_sims)
  bias_matrix <- matrix(data = NA, ncol = 5, nrow = n_sims)
  result_matrix <- matrix(data = NA, ncol = 5, nrow = n_sims)
  n <- n
  
  for (j in 1:n_sims){
    X <- runif(n,-2,2)
    Z <- rbinom(n, 1, 0.5)
    A <- rbinom(n, 1, 0.5)
    Y <- rbinom(n, 1, expit(A + 3*X + 3*Z))
    
    data <- data.frame('Y' = Y, 'X' = X, 'Z' = Z, 'A' = A)
    
    
    new_data <- data.frame('Y' = rbinom(n, 1, expit(rbinom(n, 1, 0.5) + 3*runif(n,-2,2) + 3*rbinom(n, 1, 0.5))), 
                           'X' = runif(n,-2,2), 
                           'Z' = rbinom(n, 1, 0.5), 'A' = rbinom(n, 1, 0.5))
    new_data1 <- new_data %>% mutate(A = 1)
    
    #y <- data$Y
    #x <- cbind(data$X,data$Z,data$A)
    
    #outreg <- loess(Y ~ X + Z + A, span = 0.5, degree = 2, data = data)
    #mean(Y)
    #mean(outreg$fitted)
    
    #fitting models
    outreg  <- glm(Y ~ A + X + Z - 1, data = data, family = binomial)
    propreg <- glm(A ~ X + Z, family = binomial)
    
    
    #predictions
    outpred <- expit(outreg$coefficients %*% t(cbind(1, new_data$X, new_data$Z)))
    outmis  <- rnorm(n, 0.75, 2)
    
    
    
    proppred <- propreg$fitted.values
    propmis  <- runif(n, 0.75, 1)
    
    
    PI      <-   (outpred) 
    OSBC    <-   (outpred + A*(Y - outpred)/proppred)
    OSORC   <-   (outpred + A*(Y - outpred)/propmis)
    OSPSC   <-   (outmis + A*(Y - outmis)/proppred)
    OSBM    <-   (outmis + A*(Y - outmis)/propmis)
    
    
    #mean(PI)
    #mean(OSBC)
    #mean(OSORC)
    #mean(OSPSC)
    #mean(OSBM)
    
    bias_matrix[j,1] <- 0.7028033 - mean(PI)
    bias_matrix[j,2] <- 0.7028033 - mean(OSBC)
    bias_matrix[j,3] <- 0.7028033 - mean(OSORC)
    bias_matrix[j,4] <- 0.7028033 - mean(OSPSC)
    bias_matrix[j,5] <- 0.7028033 - mean(OSBM)
    
    result_matrix[j,1] <- mean(PI)
    result_matrix[j,2] <- mean(OSBC)
    result_matrix[j,3] <- mean(OSORC)
    result_matrix[j,4] <- mean(OSPSC)
    result_matrix[j,5] <- mean(OSBM)
  }
  
  print(paste0('Rep done'))
  
  res <- c(mean(bias_matrix[,1]),
           mean(bias_matrix[,2]),
           mean(bias_matrix[,3]),
           mean(bias_matrix[,4]),
           mean(bias_matrix[,5]),
           var(result_matrix[,1]),
           var(result_matrix[,2]),
           var(result_matrix[,3]),
           var(result_matrix[,4]),
           var(result_matrix[,5]))
  
  
  return(res)
}


start.time <- Sys.time()
summary_data_frame <- data.frame()
#summary_vector <- c()

for (i in seq(300,9300,900)){
  set.seed(100)
  res_mat <- replicate(1, study(i))
  
  rowMeans(res_mat)
  
  summary_data_frame <- rbind(summary_data_frame, c(i,rowMeans(res_mat)))
}


names(summary_data_frame) <- c('n', 'NaiveBias',
                               'CorrSpecBias', 'PropMisBias',
                               'OutRegBias', 'MisBias',
                               'NaiveVar', 'CorrSpecVar',
                               'PropMisVar', 'OutRegVar',
                               'MisVar')
end.time <- Sys.time()
end.time - start.time

summary_data_frame

plot(summary_data_frame$n, summary_data_frame$CorrSpecBias, type = 'b', 
     pch = 1, ylim = c(-0.022,0.007), xlab = 'n', ylab = 'Bias')
points(summary_data_frame$n, summary_data_frame$PropMisBias, type = 'b', pch = 2)
points(summary_data_frame$n, summary_data_frame$OutRegBias, type = 'b', pch = 3)
points(summary_data_frame$n, summary_data_frame$MisBias, type = 'b', pch = 5)
abline(0,0)
legend("topright", 
       legend = c("Correctly specified", "Propensity misspecified", "Outcome regression misspecified", "Both misspecified"), 
       pch = c(1, 2, 3, 5),
       bty = "n", cex = 1.2)


plot(summary_data_frame$n, summary_data_frame$CorrSpecVar, type = 'b', 
     pch = 1, ylim = c(0,0.015), xlab = 'n', ylab = 'Variance')
points(summary_data_frame$n, summary_data_frame$PropMisVar, type = 'b', pch = 2)
points(summary_data_frame$n, summary_data_frame$OutRegVar, type = 'b', pch = 3)
points(summary_data_frame$n, summary_data_frame$MisVar, type = 'b', pch = 5)
legend("topright", 
       legend = c("Correctly specified", "Propensity misspecified", "Outcome regression misspecified", "Both misspecified"), 
       pch = c(1, 2, 3, 5),
       bty = "n", cex = 1.2)

####################################################################################################







##################################################
#             EIF Simulation Study               #
##################################################

n_sims <- 2000
simulation <- matrix(data = NA, ncol = 4, nrow = n_sims)
ba_t <- -1
tau <- 3
n <- 250

EIF_sim_func <- function(n, n_sims, tau, ba_t, S_true = T, K_true = T, prop_true = T, true_val){
  start.time <- Sys.time()
  simulation <- matrix(data = NA, ncol = 3, nrow = n_sims)
  
  EIF_coverage <- 0
  NAIVE_coverage <- 0
  
  for (j in 1:n_sims){
    ###################################################
    #           Naive estimate
    ##################################################
    data <- generate_survival_data(n, ba_t = ba_t, bx_t = log(2), bz_t = log(2),
                                   ba_c = 1, bx_c = log(2), bz_c = log(2), prop_a = 1/2, seed = sample(1:100000))
    
    #Model
    if (S_true == T){
      model_surv <- oracle_model(data)
      covar_true <- get_oracle_covar(data)
    } else {
      model_surv <- non_oracle_model(data)
      covar_true <- get_n_oracle_covar(data)
    }
    
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
    
    no_times_to_keep <- sum(jump_times_obs <= tau)
    d_cum_base_haz <- c(0, cum_bas_haz_obs[-1]-cum_bas_haz_obs[-no_jumps_obs])
    d_cum_base_haz_tau <- c(0, cum_bas_haz_obs[-1]-cum_bas_haz_obs[-no_jumps_obs])[1:no_times_to_keep]
    
    S_hat0_obs_tau <- S_hat0_obs[,1:no_times_to_keep]
    S_hat1_obs_tau <- S_hat1_obs[,1:no_times_to_keep]
    
    outer_integrand <- t(t(S_hat0_obs_tau * S_hat1_obs_tau) * d_cum_base_haz_tau)
    exp_0_multiplier <- exp(data.matrix(covar_A0) %*% beta_hat_obs)
    
    #Part 3 (naive estimate)
    p3_temp <- exp_0_multiplier * rowSums(outer_integrand)
    p4      <- mean(p3_temp)
    
    ##################################################
    
    
    
    
    ##################################################
    #             Estimating bias term
    ##################################################
    data <- generate_survival_data(n, ba_t = ba_t, bx_t = log(2), bz_t = log(2),
                                   ba_c = 1, bx_c = log(2), bz_c = log(2), prop_a = 1/2, seed = sample(1:100000))
    
    #data_for_breslow <- generate_survival_data(n, ba_t = ba_t, bx_t = log(2), bz_t = log(2),
    #                                           ba_c = 1, bx_c = log(2), bz_c = log(2), seed = sample(1:1000000))
    
    #data_for_cum_haz <- generate_survival_data(n, ba_t = ba_t, bx_t = log(2), bz_t = log(2),
    #                                           ba_c = 1, bx_c = log(2), bz_c = log(2), seed = sample(1:1000000))
    
    
    if (S_true == T){
      model_surv <- oracle_model(data)
      covar_true <- get_oracle_covar(data)
    } else {
      model_surv <- non_oracle_model(data)
      covar_true <- get_n_oracle_covar(data)
    }
    
    if (K_true == T){
      model_cens <- oracle_cens_model(data)
      covar_cens <- get_oracle_cens_covar(data)
    } else {
      model_cens <- non_oracle_cens_model(data)
      covar_cens <- get_n_oracle_covar(data)
    }
    
    if (prop_true == T){
      props <- propensity(data, glm = T)
    } else {
      props <- propensity(data, glm = F)
    }
    
    pi_0 <- props$propens$pi_a_0
    pi_1 <- props$propens$pi_a_1
    
    
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
      Khat_cens[,i] <- S_hat_cens[,max((1:no_jumps_cens)[jump_times_cens<=jump_times_obs[i]])]
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
    
    
    
    EIF_mean <- p4 + mean(p1 + p3 - p4 + p2)
    EIF_var_est <- mean((p1 + p3 - mean(p4) + p2)^2)/n
#    EIF_coverage <- EIF_coverage + (EIF_mean - 1.96*sqrt(EIF_var_est) < true_val & true_val < EIF_mean + 1.96*sqrt(EIF_var_est))
    
    NAIVE_mean <- mean(rbind(p3,p3_temp))
#    NAIVE_var_est <- var(rbind(p3,p3_temp))
#    NAIVE_coverage <- NAIVE_coverage + (NAIVE_mean - 1.96*sqrt(NAIVE_var_est) < true_val & true_val < NAIVE_mean + 1.96*sqrt(NAIVE_var_est))
    
    
    simulation[j,1] <- EIF_mean
    simulation[j,2] <- NAIVE_mean
    simulation[j,3] <- EIF_var_est
    
    
    if (j %% 100 == 0){
      print(paste0('Simulation ',j,' done')) 
    }
  }
  end.time <- Sys.time()
  print(end.time - start.time)
  
  return(c(true_val - mean(simulation[,1]), true_val - mean(simulation[,2]), mean(simulation[,3]), var(simulation[,2])))#, EIF_coverage/n_sims, NAIVE_coverage/n_sims))
}

res300 <- EIF_sim_func(n = 300, 
                    n_sims = 2000, 
                    tau = 3, 
                    ba_t = -1, 
                    S_true = T, 
                    K_true = T, 
                    prop_true = T,
                    true_val = 0.7215708370)

res600 <- EIF_sim_func(n = 600, 
                       n_sims = 2000, 
                       tau = 3, 
                       ba_t = -1, 
                       S_true = T, 
                       K_true = T, 
                       prop_true = T,
                       true_val = 0.7215708370)
















#ATE SIMULATION WITH OUTLIERS
study <- function(n){
  

n_sims <- 2000
bias_matrix <- matrix(data = NA, ncol = 5, nrow = n_sims)
result_matrix <- matrix(data = NA, ncol = 5, nrow = n_sims)
n <- n

for (j in 1:n_sims){
  X <- rnorm(n, 2, 1)
  Z <- runif(n, 0, 10)
  A <- rbernoulli(n, 0.5)
  error <- rnorm(n,0,1)
  Y <- X + Z + 4*A + error
  
  #Introducing outliers to induce small sample bias for clear results
  X_outliers <- rnorm(5, 5, 5)
  Z_outliers <- runif(5, -10, 5)
  A_outliers <- rbernoulli(5, 0.5)
  Y_outliers <-  5*X_outliers + 5*Z_outliers + 4*A_outliers
  
  data <- data.frame('Y' = c(Y, Y_outliers), 'X' = c(X, X_outliers), 'Z' = c(Z, Z_outliers), 'A' = c(A, A_outliers))
  data1 <- data %>% mutate(A = TRUE)
  data0 <- data %>% mutate(A = FALSE)
  
  linmod <- lm(Y ~ X + Z + A, data = data)
  logmod <- glm(c(A, A_outliers) ~ c(X, X_outliers), family = binomial(link = "logit"))
  
  mu_pred1 <- predict(linmod, newdata = data1)
  mu_mis1  <- rnorm(n+5, -2, 1)
  
  mu_pred0 <- predict(linmod, newdata = data0)
  mu_mis0  <- predict(linmod_mis, newdata = data0)
  
  
  pi_pred <- logmod$fitted.values
  pi_mis <- runif(n+5, min = 0.5, max = 1)
  
  
  PI      <-   (mu_pred1) 
  OSBC    <-   (mu_pred1 + c(A, A_outliers)*(c(Y, Y_outliers) - mu_pred1)/pi_pred) 
  OSORC   <-   (mu_pred1 + c(A, A_outliers)*(c(Y, Y_outliers) - mu_pred1)/pi_mis) 
  OSPSC   <-   (mu_mis1 + c(A, A_outliers)*(c(Y, Y_outliers) - mu_mis1)/pi_pred) 
  OSBM    <-   (mu_mis1 + c(A, A_outliers)*(c(Y, Y_outliers) - mu_mis1)/pi_mis) 
  
  bias_matrix[j,1] <- 11 - mean(PI)
  bias_matrix[j,2] <- 11 - mean(OSBC)
  bias_matrix[j,3] <- 11 - mean(OSORC)
  bias_matrix[j,4] <- 11 - mean(OSPSC)
  bias_matrix[j,5] <- 11 - mean(OSBM)
  
  result_matrix[j,1] <- mean(PI)
  result_matrix[j,2] <- mean(OSBC)
  result_matrix[j,3] <- mean(OSORC)
  result_matrix[j,4] <- mean(OSPSC)
  result_matrix[j,5] <- mean(OSBM)
}

print(paste0('Rep done'))

res <- c(mean(bias_matrix[,1]),
mean(bias_matrix[,2]),
mean(bias_matrix[,3]),
mean(bias_matrix[,4]),
mean(bias_matrix[,5]),
var(result_matrix[,1]),
var(result_matrix[,2]),
var(result_matrix[,3]),
var(result_matrix[,4]),
var(result_matrix[,5]))

res

return(res)
}


start.time <- Sys.time()
summary_data_frame <- data.frame()

for (i in seq(50,200,50)){
  set.seed(1)
  res_mat <- replicate(1, study(i))
  
  summary_data_frame <- rbind(summary_data_frame, c(i,rowMeans(res_mat)))
}


names(summary_data_frame) <- c('n', 'NaiveBias',
                                  'CorrSpecBias', 'PropMisBias',
                                  'OutRegBias', 'MisBias',
                                  'NaiveVar', 'CorrSpecVar',
                                  'PropMisVar', 'OutRegVar',
                                  'MisVar')
end.time <- Sys.time()
end.time - start.time



run_small_samples

plot(summary_data_frame$n, summary_data_frame$NaiveBias, type = 'b', pch = 1)
points(summary_data_frame$n, summary_data_frame$CorrSpecBias, type = 'b', pch = 2)
points(summary_data_frame$n, summary_data_frame$PropMisBias, type = 'b', pch = 3)
points(summary_data_frame$n, summary_data_frame$OutRegBias, type = 'b', pch = 4)
points(summary_data_frame$n, summary_data_frame$MisBias, type = 'b', pch = 5)

start.time <- Sys.time()
res_mat <- replicate(1, study(1000))
end.time <- Sys.time()
end.time - start.time
rowMeans(res_mat)


n_values    <-  c(50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 1000, 2000)#, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000)
naive       <- c(-1.5941919, -0.8775049, -0.5758256, -0.4171911, -0.3454646, -0.2629853, -0.2512203, -0.2202675, -0.1922048, -0.1822255, -0.08203459, -0.04354280)
corr_spec   <- c(-1.5366578, -0.8622949, -0.5559424, -0.3981158, -0.3391294, -0.2516576, -0.2440550, -0.2170266, -0.1845441, -0.1823530, -0.08190077, -0.04354355)
prop_mis    <- c(-1.5917011, -0.8795916, -0.5734816, -0.4185494, -0.3469696, -0.2606739, -0.2497822, -0.2230006, -0.1932641, -0.1837025, -0.08278444, -0.04334590)
out_reg_mis <- c(-1.3291300, -0.7885355, -0.5220827, -0.3809086, -0.3268552, -0.2372085, -0.2361884, -0.2116473, -0.1800100, -0.1815568, -0.08342307, -0.04371297)
mis_spec    <- c(-0.3163342, 0.3848232, 0.6750370, 0.8025472, 0.8842106, 0.9742906, 0.9852160, 1.0040802, 1.0383293, 1.0419013, 1.14359967, 1.18374498)

naive_var       <- c(10.9754104, 3.1568230, 1.4805774, 0.8868513, 0.6090755, 0.4382829, 0.2964046, 0.2526940, 0.2073874, 0.1768011, 0.04956658, 0.01657892)
corr_spec_var   <- c(11.4440419, 3.4640062, 1.5354914, 0.9192067, 0.6464731, 0.4701266, 0.3172775, 0.2768897, 0.2181831, 0.1863968, 0.05183884, 0.01683329)
prop_mis_var    <- c(11.0650669, 3.1762310, 1.4962839, 0.8964394, 0.6130945, 0.4418999, 0.2993455, 0.2583992, 0.2107911, 0.1789154, 0.05071778, 0.01687444)
out_reg_mis_var <- c(13.6607109, 4.7147113, 1.7702785, 1.0533592, 0.7653503, 0.5607253, 0.3676652, 0.3378899, 0.2513753, 0.2118299, 0.06115317, 0.02100236)
mis_spec_var    <- c(9.4532938, 2.6367472, 1.2045705, 0.7389088, 0.4953927, 0.3624715, 0.2504770, 0.2136922, 0.1750629, 0.1444112, 0.04543471, 0.01715284)



plot(n_values, naive, type = 'b', pch = 1, ylim = c(-1.5,0))
points(n_values, corr_spec, type = 'b', pch = 2)
points(n_values, prop_mis, type = 'b', pch = 3)
points(n_values, out_reg_mis, type = 'b', pch = 4)
points(n_values, mis_spec, type = 'b', pch = 5)
abline(0,0)



plot(n_values, naive_var*sqrt(n_values), type = 'b', pch = 1, ylim = c(0,85))
points(n_values, corr_spec_var*sqrt(n_values), type = 'b', pch = 2)
points(n_values, prop_mis_var*sqrt(n_values), type = 'b', pch = 3)
points(n_values, out_reg_mis_var*sqrt(n_values), type = 'b', pch = 4)
points(n_values, mis_spec_var*sqrt(n_values), type = 'b', pch = 5)
abline(0,0)


n_sims <- 2000
simulation <- matrix(data = NA, ncol = 6, nrow = n_sims)
ba_t <- -1
tau <- 3
n <- 500

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
  
  
  props <- propensity(data, glm = T)
  
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
  
  #par(mfrow = c(1,1))
  #plot(jump_times_obs, cumsum(colSums(dM)), type = 'l', ylab = 'Kummuleret dM.')
  #plot(jump_times_obs, cumsum(colSums(dN)), type = 'l', ylab = 'Kummuleret dN., dL.', col = 'red')
  #lines(jump_times_obs, cumsum(colSums(at_risk*dL)), col = alpha('blue', 0.2))
  
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
  
  
  NAIVE_var_est <- var(p3)/n
  EIF_var_est <- mean((p1 + p3 - mean(p3) + p2)^2)/n
  EIF_mean <- mean(p1 + p3 - mean(p3) + p2)
  
  
  
  #simulation[j,1] <- mean(p1)
  #simulation[j,2] <- mean(p2)
  #simulation[j,3] <- mean(p3)
  simulation[j,4] <- mean(p1 + p2 + p3)
  #simulation[j,5] <- NAIVE_var_est
  simulation[j,6] <- EIF_var_est
  simulation[j,7] <- EIF_mean
  #print(paste0('Simulation ',j,' done'))
}
end.time <- Sys.time()
end.time - start.time
mean(simulation[,3])
mean(simulation[,4])
mean(simulation[,5])
mean(simulation[,6])


#True value: 0.7215708370

#1.935071e-05

#0.6078786 0.6123083

x_axis <- c(50, 100, 200, 400, 800, 1600)#, 3200)
naive_bias <- c(0.7215708370-0.6086665, 
                0.7215708370-0.6607547, 
                0.7215708370-0.68665, 
                0.7215708370-0.6959826, 
                0.7215708370-0.7094471, 
                0.7215708370-0.7177696)#, 
                #0.7215708370-0.7191623)
eif_bias <- c(0.7215708370-0.6123555,
              0.7215708370-0.6688572,
              0.7215708370-0.6865309,
              0.7215708370-0.6946947,
              0.7215708370-0.7095038,
              0.7215708370-0.7174959)
#,
#              0.7215708370-0.7195862)
eif_var <- c(0.01112825,
             0.005346707,
             0.003019546,
             0.001465691,
             0.0007382434,
             0.0003714565,
             0.0001866196)

plot(x_axis, log(naive_bias), type = 'b', pch = 2)
points(x_axis, log(eif_bias), type = 'b', pch = 5)

plot(x_axis, eif_var, type = 'b', pch = 2)

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


#             n     50          100             200             400           800             1600          3200
#Naive estimate   0.6086665     0.6607547       0.68665         0.6959826     0.7094471       0.7177696     0.7191623
#EIF estimate     0.6123555     0.6688572       0.6865309       0.6946947     0.7095038       0.7174959     0.7195862
#EIF variance     0.01112825    0.005346707     0.003019546     0.001465691   0.0007382434    0.0003714565  0.0001866196



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



#Simulations for full data ATE
n_sims <- 50
ba_t <- -1
tau <- 2
n <- 3000

start.time <- Sys.time()
for (j in 1:n_sims){
  data <- generate_survival_data(n, ba_t = ba_t, bx_t = log(2), bz_t = log(2),
                                 ba_c = 1, bx_c = log(2), bz_c = log(2), prop_a = 1/2, seed = sample(1:100000))
  
  
  model_surv <- cox.aalen(Surv(T_true, rep(1, nrow(data))) ~ prop(A) + prop(X) + prop(Z) , data = data)
  model_surv_mis <- cox.aalen(Surv(T_true, rep(1, nrow(data))) ~ prop(A) + prop(fake_2) + prop(Z*A) + prop(X^2) , data = data)
  
  outreg      <- 1/(1 + exp(model_surv$gamma[1]))
  outreg_mis  <- 1/(1 + exp(model_surv_mis$gamma[1]))
  prop_true   <- 0.5
  prop_mis    <- rnorm(n, 1, 1)
  
  NAIVE <- outreg
  OS <- mean(outreg + data$A/prop_true*(data$T_true - outreg))
  
  
  
  
  
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
  
  #plot(jump_times_obs, colMeans(S_hat_obs), type = 'l')
  #lines(jump_times_obs, colMeans(S_hat0_obs), col = 'red')
  #lines(jump_times_obs, colMeans(S_hat1_obs), col = 'blue')
  
  
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
  #
  #plot(jump_times_obs, cumsum(dM[232,]), type = 'l', ylim = c(-3,1.1))
  #lines(jump_times_obs, cumsum(dN[232,]))
  #lines(jump_times_obs, -cumsum(at_risk[232,]*dL[232,]))
  
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
  
  OS <- S_hat1_obs_tau - (pi_1)*inner_integral
  
  #The estimator
  one_step_estimator <- colMeans(OS)
  
  
  new_element <- list(jump_times = jump_times_kept, one_step_estimator = one_step_estimator)
  simulation_list[[length(simulation_list) + 1]] <- new_element
  print(paste0('Simulation ',j,' done'))
}
end.time <- Sys.time()
end.time - start.time




#Simulations for ATE

n_sims <- 50
simulation_list <- list()
ba_t <- -1
tau <- 2
n <- 3000

start.time <- Sys.time()
for (j in 1:n_sims){
  data <- generate_survival_data(n, ba_t = ba_t, bx_t = log(2), bz_t = log(2),
                                 ba_c = 1, bx_c = log(2), bz_c = log(2), prop_a = 1/2, seed = sample(1:100000))
  
  
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
  
  #plot(jump_times_obs, colMeans(S_hat_obs), type = 'l')
  #lines(jump_times_obs, colMeans(S_hat0_obs), col = 'red')
  #lines(jump_times_obs, colMeans(S_hat1_obs), col = 'blue')
  
  
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
  #
  #plot(jump_times_obs, cumsum(dM[232,]), type = 'l', ylim = c(-3,1.1))
  #lines(jump_times_obs, cumsum(dN[232,]))
  #lines(jump_times_obs, -cumsum(at_risk[232,]*dL[232,]))
  
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
  
  OS <- S_hat1_obs_tau - (pi_1)*inner_integral
  
  #The estimator
  one_step_estimator <- colMeans(OS)
  
  
  new_element <- list(jump_times = jump_times_kept, one_step_estimator = one_step_estimator)
  simulation_list[[length(simulation_list) + 1]] <- new_element
  print(paste0('Simulation ',j,' done'))
}
end.time <- Sys.time()
end.time - start.time




#True values bx = bz = log(2), ba = -1, tau = 3
taus <- seq(0,2,0.05)
true_values <- c(1, 0.9717751250, 0.9445018885, 0.9181431865, 0.8926635295, 0.8680289690, 0.8442070210, 0.8211666020, 0.7988779580, 0.7773126080, 0.7564432805, 0.7362438595, 0.7166893300, 0.6977557285, 0.6794200945, 0.6616604230, 0.6444556235, 0.6277854760, 0.6116305930, 0.5959723825, 0.5807930095, 0.5660753650, 0.5518030325, 0.5379602570, 0.5245319145, 0.5115034875, 0.4988610342, 0.4865911662, 0.4746810230, 0.4631182490, 0.4518909722, 0.4409877824, 0.4303977123, 0.4201102176, 0.4101151594, 0.4004027866, 0.3909637198, 0.3817889354, 0.3728697502, 0.3641978078, 0.3557650644)
truth <- data.frame('taus' = taus, 'true_vals' = true_values)

plot(simulation_list[[1]][[1]], 
     simulation_list[[1]][[2]], 
     xlab = 'Time', 
     ylab = 'Estimate', 
     type = 'l', 
     col = rgb(red = 0.5, green = 0, blue = 0.7, alpha = 0.3))
for (i in 1:n_sims){
  lines(simulation_list[[i]][[1]], simulation_list[[i]][[2]],
        col = rgb(red = 0.5, green = 0, blue = 0.7, alpha = 0.3))
}
lines(truth$taus, truth$true_vals)














################################################################
#                         COXPH
################################################################


################################################################
#               Full data one-step estimator of ATE
################################################################

n_sims <- 250
ba_t <- -1
tau <- 3
n <- 900

NAIVE_lst   <- list()
OSECOR_lst  <- list()
OSEMOR_lst  <- list()
OSEMPS_lst  <- list()


for (j in 1:n_sims){
  data <- generate_survival_data(n, ba_t = ba_t, bx_t = log(2), bz_t = log(2),
                                 ba_c = 1, bx_c = log(2), bz_c = log(2), seed = sample(1:100000))
  data_A0 <- data %>% mutate(A = 0)
  data_A1 <- data %>% mutate(A = 1)
  status <- rep(1, nrow(data))
  
  ######################################################
  
  ######################################################
  fit_coxph_full_data     <- coxph(Surv(T_true, status) ~ A + X + Z, data = data)
  fit_coxph_full_data_mis <- lm(T_true ~ A + X + Z - 1, data = data)
  ######################################################
  
  ######################################################
  surv_fit   <- survfit(fit_coxph_full_data, newdata = data)
  surv_fit_0 <- survfit(fit_coxph_full_data, newdata = data_A0)
  surv_fit_1 <- survfit(fit_coxph_full_data, newdata = data_A1)
  
  surv    <- t(surv_fit$surv)
  surv_0  <- t(surv_fit_0$surv)
  surv_1  <- t(surv_fit_1$surv)
  ######################################################
  
  ######################################################
  surv_mis    <- predict(fit_coxph_full_data_mis, newdata = data)
  surv_0_mis  <- predict(fit_coxph_full_data_mis, newdata = data_A0)
  surv_1_mis  <- predict(fit_coxph_full_data_mis, newdata = data_A1)
  ######################################################
  
  ######################################################
  props <- propensity(data, glm = T)
  pi_0 <- props$propens$pi_a_0
  pi_1 <- props$propens$pi_a_1
  
  props_mis <- propensity(data, glm = F)
  pi_1_mis <- props_mis$propens$pi_a_1
  
  ######################################################
  
  at_risk <- outer(data$T_true, survfit$time, '>=')
  
  NAIVE   <- colMeans(surv_1)
  OSECOR  <- colMeans(surv_1 + t(pi_1*t(at_risk - surv_1)))
  OSEMOR  <- colMeans(surv_1_mis + t(pi_1*t((at_risk - surv_1_mis))))
  OSEMPS  <- colMeans(surv_1 + t(pi_1_mis*t((at_risk - surv_1))))
  
  
  NAIVE_new_element         <- list(jump_times = survfit$time, one_step_estimator = NAIVE)
  OSECOR_new_element        <- list(jump_times = survfit$time, one_step_estimator = OSECOR)
  OSEMOR_new_element        <- list(jump_times = survfit$time, one_step_estimator = OSEMOR)
  OSEMPS_new_element        <- list(jump_times = survfit$time, one_step_estimator = OSEMPS)
  
  
  NAIVE_lst[[length(NAIVE_lst) + 1]]    <- NAIVE_new_element
  OSECOR_lst[[length(OSECOR_lst) + 1]]  <- OSECOR_new_element
  OSEMOR_lst[[length(OSEMOR_lst) + 1]]  <- OSEMOR_new_element
  OSEMPS_lst[[length(OSEMPS_lst) + 1]]  <- OSEMPS_new_element
  
  print(paste0('Simulation ',j,' done'))
}

taus <- seq(0,10,0.1)
true_values <- c(1, 0.9445018885, 0.8926635295, 0.8442070210, 0.7988779580, 0.7564432805, 0.7166893300, 0.6794200945, 0.6444556235, 0.6116305930, 0.5807930095, 0.5518030325, 0.5245319145, 0.4988610342, 0.4746810230, 0.4518909722, 0.4303977123, 0.4101151594, 0.3909637198, 0.3728697502, 0.3557650644, 0.3395864864, 0.3242754414, 0.3097775838, 0.2960424576, 0.2830231864, 0.2706761908, 0.2589609291, 0.2478396611, 0.2372772317, 0.2272408724, 0.2177000200, 0.2086261504, 0.1999926251, 0.1917745521, 0.1839486566, 0.1764931627, 0.1693876852, 0.1626131289, 0.1561515970, 0.1499863059, 0.1441015072, 0.1384824156, 0.1331151419, 0.1279866322, 0.1230846106, 0.1183975268, 0.1139145076, 0.1096253120, 0.1055202892, 0.1015903404, 0.09782688290, 0.09422181685, 0.09076749430, 0.08745669090, 0.08428257910, 0.08123870330, 0.07831895710, 0.07551756175, 0.07282904615, 0.07024822825, 0.06777019785, 0.06539030030, 0.06310412150, 0.06090747370, 0.05879638235, 0.05676707385, 0.05481596400, 0.05293964720, 0.05113488645, 0.04939860376, 0.04772787140, 0.04611990356, 0.04457204858, 0.04308178164, 0.04164669786, 0.04026450592, 0.03893302196, 0.03765016386, 0.03641394598, 0.03522247400, 0.03407394028, 0.03296661932, 0.03189886360, 0.03086909958, 0.02987582396, 0.02891760016, 0.02799305502, 0.02710087560, 0.02623980623, 0.02540864572, 0.02460624472, 0.02383150314, 0.02308336787, 0.02236083049, 0.02166292516, 0.02098872656, 0.02033734808, 0.01970793992, 0.01909968745, 0.01851180954)
truth <- data.frame('taus' = taus, 'true_vals' = true_values)

plot(OSECOR_lst[[1]][[1]], 
     OSECOR_lst[[1]][[2]], 
     main = 'Correctly specified',
     xlab = 'Time', 
     ylab = 'Estimate', 
     type = 'l', 
     col = rgb(red = 0.5, green = 0, blue = 0.7, alpha = 0.3))
for (i in 2:length(OSECOR_lst)){
  lines(OSECOR_lst[[i]][[1]], OSECOR_lst[[i]][[2]],
        col = rgb(red = 0.5, green = 0, blue = 0.7, alpha = 0.3))
}
lines(truth$taus, truth$true_vals)


plot(OSEMOR_lst[[1]][[1]], 
     OSEMOR_lst[[1]][[2]], 
     main = 'Misspecified survival',
     xlab = 'Time', 
     ylab = 'Estimate', 
     type = 'l', 
     col = rgb(red = 0.5, green = 0, blue = 0.7, alpha = 0.3),
     ylim = c(0,1))
for (i in 2:length(OSEMOR_lst)){
  lines(OSEMOR_lst[[i]][[1]], OSEMOR_lst[[i]][[2]],
        col = rgb(red = 0.5, green = 0, blue = 0.7, alpha = 0.3))
}
lines(truth$taus, truth$true_vals)

plot(OSEMPS_lst[[1]][[1]], 
     OSEMPS_lst[[1]][[2]], 
     main = "Misspecified propensity",
     xlab = 'Time', 
     ylab = 'Estimate', 
     type = 'l', 
     col = rgb(red = 0.5, green = 0, blue = 0.7, alpha = 0.3),
     ylim = c(0,1))
for (i in 2:length(OSEMPS_lst)){
  lines(OSEMPS_lst[[i]][[1]], OSEMPS_lst[[i]][[2]],
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
  
  data_cens <- generate_survival_data(n, ba_t = ba_t, bx_t = log(2), bz_t = log(2),
                                      ba_c = 1, bx_c = log(2), bz_c = log(2), seed = sample(1:100000))
  
  data_0 <- data %>% mutate(A = 0, X = 0, Z = 0)
  data_A0 <- data %>% mutate(A = 0)
  data_A1 <- data %>% mutate(A = 1)
  
  
  ################################################################
  #                     Fitting models
  ################################################################
  fit_coxph <- coxph(Surv(T_obs, Uncensored) ~ A + X + Z, data = data)
  fit_coxph_mis <- coxph(Surv(T_obs, Uncensored) ~ A + fake_1 + fake_2 + fake_3, data = data)
  
  cum_base_haz <- basehaz(fit_coxph, newdata = data_0)
  cum_base_haz_mis <- basehaz(fit_coxph_mis, newdata = data_0)
  
  survfit_all <- survfit(fit_coxph, newdata = data)
  survfit_all_mis <- survfit(fit_coxph_mis, newdata = data)
  
  #plot(survfit_all$time, cumsum(fit_coxph$residuals), type = 'l', col = 'red')
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
  fit_coxph_cens <- coxph(Surv(T_obs, Uncensored == F) ~ A + X + Z, data = data_cens)
  censfit <- survfit(fit_coxph_cens, newdata = data)
  cens <- censfit$surv
  
  fit_coxph_cens_mis <- coxph(Surv(T_obs, Uncensored == F) ~ A + fake_1 + fake_2 + fake_3, data = data_cens)
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





#OUR EIF WITH COXPH
n_sims <- 1000
ba_t <- -1
tau <- 3
n <- 900


#Data generation
data <- generate_survival_data(n, ba_t = ba_t, bx_t = log(2), bz_t = log(2),
                               ba_c = 1, bx_c = log(2), bz_c = log(2), seed = sample(1:100000))

data_0 <- data %>% mutate(A = 0, X = 0, Z = 0)
data_A0 <- data %>% mutate(A = 0)
data_A1 <- data %>% mutate(A = 1)


################################################################
#                     Fitting models
################################################################
fit_coxph <- coxph(Surv(T_obs, Uncensored) ~ A + X + Z, data = data)

cum_base_haz <- t(basehaz(fit_coxph, newdata = data_0))[-nrow(data) - 1,]




################################################################
#           Survival functions and cumulative hazards
################################################################
survfit     <-  survfit(fit_coxph, newdata = data)
survfit_0   <-  survfit(fit_coxph, newdata = data_A0)
survfit_1   <-  survfit(fit_coxph, newdata = data_A1)

surv        <- t(survfit$surv)
surv_0      <- t(survfit_0$surv)
surv_1      <- t(survfit_1$surv)

cum_haz     <- t(survfit$cumhaz)
cum_haz_0   <- t(survfit_0$cumhaz)
cum_haz_1   <- t(survfit_1$cumhaz)


#plot(survfit_0$time, surv_0[1,], type = 'l', col = 'red')
#for (i in 2:length(survfit_0$time)){
#  lines(survfit_0$time, surv_0[i,], col = 'red')
#  lines(survfit_1$time, surv_1[i,], col = 'blue')
#}

################################################################
#                     Martingale errors
################################################################

#Note these martingales are aggregated across all individuals
dM <- fit_coxph$residuals

#----------------Manual calculations--------------
#dL <- cbind(0, cum_haz[,-1] - cum_haz[,-ncol(cum_haz)])
#
#dN <- outer(data$T_obs, survfit$time, '==') * data$Uncensored
#
#at_risk <- outer(data$T_obs, survfit$time, '>=')
#dM_manual <- t(dN - dL*at_risk)
#
#plot(survfit$time, cumsum(colSums(dN)), type = 'l')
#lines(survfit$time, cumsum(colSums(at_risk*dL)))
#
#plot(survfit$time, cumsum(dM), type = 'l')
#lines(survfit$time, cumsum(colSums(dM_manual)))


################################################################
#                     Censoring distribution
################################################################
fit_coxph_cens <- coxph(Surv(T_obs, Uncensored == F) ~ A + X + Z, data = data_cens)
censfit <- survfit(fit_coxph_cens, newdata = data)
cens <- t(censfit$surv)




################################################################
#                           EIF
################################################################
time_indicator <- survfit_all$time <= tau
jump_times_kept <- survfit_all$time[time_indicator]


multiplier <- exp(fit_coxph$coefficients %*% t(cbind(data$A, data$X, data$Z)))

dH0 <- cbind(0, cum_base_haz[,-1] - cum_base_haz[,-ncol(cum_base_haz)])
dH0_tau <- dH0[,time_indicator]

surv_tau    <-  surv[,time_indicator]
surv_0_tau  <-  surv_0[,time_indicator] 
surv_1_tau  <-  surv_1[,time_indicator]

Khat_cens_tau <- t(cens)[,time_indicator]

dM_tau        <- dM[time_indicator]
dM_manual_tau <- dM_manual[, time_indicator]


props <- propensity(data, glm = T)
pi_0 <- props$propens$pi_a_0
pi_1 <- props$propens$pi_a_1


#Naive estimate
mean(multiplier * rowSums(surv_0_tau*surv_1_tau*dH0_tau))

#One-step estimator



#Correctly specified
inner_integrand <- dM_manual_tau/(surv_tau*Khat_cens_tau)
inner_integral <- t(apply(inner_integrand, MARGIN = 1, FUN = cumsum))



















##############LEGACY##############


# ba_t <- -1
# n <- 1000
# data <- generate_survival_data(n, ba_t = ba_t, bx_t = log(2), bz_t = log(2),
#                                ba_c = 1, bx_c = log(2), bz_c = log(2), seed = sample(1:100000))
# 
# 
# model <- cox.aalen(Surv(T_true, rep(1, nrow(data))) ~ prop(A) + prop(X^2) + prop(Z^2), data)
# covar <- get_oracle_covar(data)
# 
# props <- propensity(data)
# 
# cum_matrix <- get_cum(model)
# jump_times <- get_jump_times(cum_matrix)
# cum_bas_haz <- get_cum_base_haz(cum_matrix)
# no_jumps <- get_row_length(cum_matrix)
# beta_hat <- get_param_est(model)
# 
# 
# covar_A0 <- covar %>% mutate(A = 0)
# covar_A1 <- covar %>% mutate(A = 1)
# 
# #Cum haz obs
# cum_hat_obs <- predict_cox.aalen(covar = covar, betaHat = beta_hat, cum_base_haz = cum_bas_haz)
# cum_hat0_obs <- predict_cox.aalen(covar = covar_A0, betaHat = beta_hat, cum_base_haz = cum_bas_haz)
# cum_hat1_obs <- predict_cox.aalen(covar = covar_A1, betaHat = beta_hat, cum_base_haz = cum_bas_haz)
# 
# #Survival functions
# S_hat <- exp(-cum_hat_obs)
# S_hat0 <- exp(-cum_hat0_obs)
# S_hat1 <- exp(-cum_hat1_obs)
# at_risk <- outer(data$T_true, jump_times, '>=')
# 
# mu_hat <- exp(-cum_hat1_obs)
# pi_hat <- 1/2  #props$probs[,1]
# at_risk <- outer(data$T_true, jump_times, '>=')
# A <- data$A
# 
# frac <- (A/pi_hat) * (at_risk - mu_hat)
# 
# mean(mu_hat + frac)
# mean(mu_hat)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ba_t <- log(2)
# n <- 1000
# data <- generate_survival_data(n, ba_t = ba_t, bx_t = log(2), bz_t = log(2),
#                                ba_c = 1, bx_c = log(2), bz_c = log(2), seed = sample(1:100000))
# 
# 
# model <- cox.aalen(Surv(T_true, rep(1, nrow(data))) ~ prop(A) + prop(X) + prop(Z), data)
# covar <- get_oracle_covar(data)
# 
# props <- propensity(data)
# 
# cum_matrix <- get_cum(model)
# jump_times <- get_jump_times(cum_matrix)
# cum_bas_haz <- get_cum_base_haz(cum_matrix)
# no_jumps <- get_row_length(cum_matrix)
# beta_hat <- get_param_est(model)
# 
# 
# covar_A0 <- covar %>% mutate(A = 0)
# covar_A1 <- covar %>% mutate(A = 1)
# 
# #Cum haz obs
# cum_hat_obs <- predict_cox.aalen(covar = covar, betaHat = beta_hat, cum_base_haz = cum_bas_haz)
# cum_hat0_obs <- predict_cox.aalen(covar = covar_A0, betaHat = beta_hat, cum_base_haz = cum_bas_haz)
# cum_hat1_obs <- predict_cox.aalen(covar = covar_A1, betaHat = beta_hat, cum_base_haz = cum_bas_haz)
# 
# #Survival functions
# S_hat <- exp(-cum_hat_obs)
# S_hat0 <- exp(-cum_hat0_obs)
# S_hat1 <- exp(-cum_hat1_obs)
# at_risk <- outer(data$T_true, jump_times, '>=')
# dN <- outer(data$T_true, jump_times, '==')
# exp_0_multiplier <- exp(data.matrix(covar_A0) %*% beta_hat)
# d_cum_base_haz <- c(0, cum_bas_haz[-1]-cum_bas_haz[-no_jumps])
# pi_0 <- props$propens$pi_a_0
# pi_1 <- props$propens$pi_a_1
# 
# time_ind <- jump_times < 2
# 
# 
# part1 <- pi_0 * exp_0_multiplier * rowSums(t(t(at_risk*S_hat0) * d_cum_base_haz)[,time_ind])
# part2 <- pi_1 * rowSums((S_hat1 * dN)[,time_ind])
# part3 <- - (pi_0 + pi_1)*P_treatment_extend_survival(model = model, max_time = 3, model_cov = covar)$res
# part4 <- P_treatment_extend_survival(model = model, max_time = 3, model_cov = covar)$res
# 
# 1 - mean(part1 + part2 + part3 + part4)
# mean(part4)
# 
# 
# 
# 
# 
# 
# 
# 
# n_sims <- 2000
# simulation <- matrix(data = NA, ncol = 2, nrow = n_sims)
# ba_t <- -1
# tau <- 3
# n <- 900
# 
# for (j in 1:n_sims){
#   data <- generate_survival_data(n, ba_t = ba_t, bx_t = log(2), bz_t = log(2),
#                                  ba_c = 1, bx_c = log(2), bz_c = log(2), seed = sample(1:100000))
#   
#   
#   model_surv <- oracle_model(data)
#   model_cens <- oracle_cens_model(data)
#   
#   covar_true <- get_oracle_covar(data)
#   covar_cens <- get_oracle_cens_covar(data)
#   
#   
#   props <- propensity(data)
#   
#   
#   
#   cum_matrix_obs <- get_cum(model_surv)
#   jump_times_obs <- get_jump_times(cum_matrix_obs)
#   cum_bas_haz_obs <- get_cum_base_haz(cum_matrix_obs)
#   no_jumps_obs <- get_row_length(cum_matrix_obs)
#   beta_hat_obs <- get_param_est(model_surv)
#   
#   
#   #breslow <- breslow_estimator(data_for_breslow, beta_hat_obs, T_corr = F)
#   #cum_bas_haz_obs_breslow <- breslow$breslow
#   #breslow_jump_times <- breslow$jump_times_breslow
#   #no_jump_times_breslow <- length(breslow_jump_times)
#   #
#   #cum_bas_haz_obs = rep(NA,no_jumps_obs) 
#   #for(i in 1:no_jumps_obs){
#   #  cum_bas_haz_obs[i] = cum_bas_haz_obs_breslow[max((1:no_jump_times_breslow)[breslow_jump_times<=jump_times_obs[i]])]
#   #}
#   
#   #plot(jump_times_obs, cum_bas_haz_obs, type = 'l')
#   #lines(jump_times_obs, cum_bas_haz_obs_old)
#   #lines(breslow_jump_times, cum_bas_haz_obs_breslow)
#   
#   covar_A0 <- covar_true %>% mutate(A = 0)
#   covar_A1 <- covar_true %>% mutate(A = 1)
#   
#   #Cum haz obs
#   cum_hat_obs <- predict_cox.aalen(covar = covar_true, betaHat = beta_hat_obs, cum_base_haz = cum_bas_haz_obs)
#   cum_hat0_obs <- predict_cox.aalen(covar = covar_A0, betaHat = beta_hat_obs, cum_base_haz = cum_bas_haz_obs)
#   cum_hat1_obs <- predict_cox.aalen(covar = covar_A1, betaHat = beta_hat_obs, cum_base_haz = cum_bas_haz_obs)
#   
#   #Survival functions
#   S_hat_obs <- exp(-cum_hat_obs)
#   S_hat0_obs <- exp(-cum_hat0_obs)
#   S_hat1_obs <- exp(-cum_hat1_obs)
#   
#   
#   #Censoring
#   cum_matrix_cens <- get_cum(model_cens)
#   jump_times_cens <- get_jump_times(cum_matrix_cens)
#   cum_bas_haz_cens <- get_cum_base_haz(cum_matrix_cens)
#   no_jumps_cens <- get_row_length(cum_matrix_cens)
#   beta_hat_cens <- get_param_est(model_cens)
#   
#   
#   #breslow_cens <- breslow_estimator_cens(data_for_breslow, beta_hat_cens, cens_corr = T)
#   #cum_bas_haz_obs_breslow_cens <- breslow_cens$breslow
#   #breslow_jump_times_cens <- breslow_cens$jump_times_breslow
#   #no_jump_times_breslow_cens <- length(breslow_jump_times_cens)
#   
#   #plot(jump_times_cens, cum_bas_haz_cens, type = 'l')
#   #lines(breslow_jump_times_cens, cum_bas_haz_obs_breslow_cens)
#   #lines(jump_times_cens, cum_bas_haz_cens)
#   
#   #cum_bas_haz_cens = rep(NA,no_jumps_cens) 
#   #for(i in 1:no_jumps_cens){
#   #  cum_bas_haz_cens[i] = cum_bas_haz_obs_breslow_cens[max((1:no_jump_times_breslow_cens)[breslow_jump_times_cens<=jump_times_cens[i]])]
#   #}
#   
#   
#   cum_hat_cens <- predict_cox.aalen(covar = covar_cens, betaHat = beta_hat_cens, cum_base_haz = cum_bas_haz_cens)
#   
#   S_hat_cens <- exp(- cum_hat_cens)
#   
#   
#   Khat_cens = matrix(NA, nrow = n, ncol = no_jumps_obs) 
#   for(i in 1:no_jumps_obs){
#     Khat_cens[,i] = S_hat_cens[,max((1:no_jumps_cens)[jump_times_cens<=jump_times_obs[i]])]
#   }
#   
#   #plot(jump_times_obs, colMeans(Khat_cens), type = 'l')
#   #lines(jump_times_obs, Khat_cens[1,])
#   #lines(jump_times_cens, colMeans(S_hat_cens))
#   
#   
#   #Martingales
#   T_obs <- get_observed_times(data)
#   
#   dN <- outer(T_obs, jump_times_obs, '==')
#   
#   at_risk <- outer(T_obs, jump_times_obs, '>=')
#   dL <- cbind(0, cum_hat_obs[,-1] - cum_hat_obs[,-ncol(cum_hat_obs)])
#   
#   dM <- dN - at_risk * dL
#   
#   #par(mfrow = c(1,1))
#   #plot(jump_times_obs, cumsum(colSums(dM)), type = 'l', ylab = 'Kummuleret dM.')
#   #plot(jump_times_obs, cumsum(colSums(dN)), type = 'l', ylab = 'Kummuleret dN., dL.', col = 'red')
#   #lines(jump_times_obs, cumsum(colSums(at_risk*dL)), col = alpha('blue', 0.2))
#   
#   #plot(jump_times_obs, cumsum(dM[233,]), type = 'l', ylim = c(-3,1.1))
#   #lines(jump_times_obs, cumsum(dN[233,]))
#   #lines(jump_times_obs, -cumsum(at_risk[233,]*dL[233,]))
#   
#   #EIF
#   
#   tau <- 3
#   
#   d_cum_base_haz <- c(0, cum_bas_haz_obs[-1]-cum_bas_haz_obs[-no_jumps_obs])
#   pi_0 <- props$propens$pi_a_0
#   pi_1 <- props$propens$pi_a_1
#   
#   S_hat1_obs_tau <- S_hat1_obs[,jump_times_obs <= tau]
#   S_hat0_obs_tau <- S_hat0_obs[,jump_times_obs <= tau]
#   S_hat_obs_tau <- S_hat_obs[,jump_times_obs <= tau]
#   Khat_cens_tau <- Khat_cens[,jump_times_obs <= tau]
#   dM_tau <- dM[,jump_times_obs <= tau]
#   
#   
#   delta <- S_hat1_obs_tau - S_hat0_obs_tau
#   
#   g <- pi_1 - pi_0
#   
#   
#   inner_integrand <- 1/(S_hat_obs_tau * Khat_cens_tau) * dM_tau
#   integral <- t(apply(inner_integrand, MARGIN = 1, FUN = cumsum))
#   
#   h <- g*S_hat_obs_tau * integral
#   EIF <- delta + h
#   
#   
#   simulation[j,1] <- mean(EIF)
#   simulation[j,2] <- mean(delta)
#   print(paste0('Simulation ',j,' done'))
# }
# mean(simulation[,1])
# mean(simulation[,2])
# 
# 
# 
# 
# 
# #FF EIF: 0.2328531
# #FF naive: 0.2350929
# 
# #TT EIF: 0.1807568
# #TT naive: 0.1798884
# 
# 
# #TF EIF: 0.2287061
# #TF naive: 0.2350452
# 
# 
# 
# 
# 
# 
