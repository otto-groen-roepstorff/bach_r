####Packages----
library(targets)
library(tidyverse)
library(survival)
library(timereg)
library(dplyr)
library(data.table)
library(survival)


#Get functions----------------------
#Optimized for generation method----

get_row_length <- function(x){
  return(nrow(x))
}


get_X_val <- function(data){
  return(data$X)
}

get_A_val <- function(data){
  return(data$A)
}

get_Z_val <- function(data){
  return(data$Z)
}

get_censor_status <- function(data){
  return(data$Uncensored)
}

get_observed_times <- function(data){
  return(data$T_obs)
}

get_C_times <- function(data){
  return(data$C)
}

get_nA_covariates_val <- function(data){
  return(data[, !names(data) %in% c("T_true", "C", "T_obs", "Uncensored", "A")])
}


get_cum <- function(model){
  return(model$cum)
}

get_jump_times <- function(cum_haz_matrix){
  return(cum_haz_matrix[,"time"])
}

get_cum_base_haz <- function(cum_haz_matrix){
  return(cum_haz_matrix[,"(Intercept)"])
}

get_param_est <- function(cox.aalen_model){
  return(cox.aalen_model$gamma)
}

get_mg_dM <- function(martingale_list){
  return(martingale_list$mod_dM)
}

get_mg_dL <- function(martingale_list){
  return(martingale_list$mod_dL)
}

#Generate survival data------
generate_survival_data <- function(n, 
                                   ba_t = -log(2), bx_t = log(2), bz_t = -log(2),
                                   ba_c = log(2), bx_c = log(2), bz_c = -log(2),
                                   x_range = c(-1,1), x_cont = T, x_vals = c(0,1), 
                                   prop_a = 1/2, seed = 1){
  set.seed(seed)
  #Treatment indicators
  a <- rbinom(n, 1, prop_a) 
  
  #Generate list of fake variables to increase variance and bias in estimates
  z <- runif(n = n, min = 0, max = 1)
  fake_1 <- rbinom(n, 1, prop_a*5/6) 
  fake_2 <- rbinom(n, 1, prop_a*1/3)
  fake_3 <- runif(n = n, min = 0, max = 10)
  
  if(x_cont){
    x <- runif(n, min = min(x_range), max = max(x_range))
  }
  else{
    x <- sample(x_vals, n, replace = T)
  }
  
  #Generating survival times
  surv_t <- rexp(n, 1)/exp(ba_t*a + bx_t*x + bz_t*z) 
  
  #Generating censoring times
  cens_t <- rexp(n, 1)/exp(bz_c*z + bx_c*x)
  
  #Observed values
  t_obs <- pmin(surv_t, cens_t)
  
  #Status
  delta <- surv_t < cens_t
  
  return(data.frame(T_true = surv_t, C = cens_t, T_obs = t_obs, Uncensored = delta, A = a, X = x, Z = z, fake_1 = fake_1, fake_2 = fake_2, fake_3 = fake_3))
}


#Summary statistics ----------------

proportion_observed <- function(data){
  mean(get_censor_status(data))
}

visualize_data <- function(data){
  df1 <- pivot_longer(data = data, cols = c("T_true", "C", "T_obs"), names_to = "Source", )
  df1$A <- as.factor(df1$A)
  ggplot(df1)+ geom_density(aes(value, fill = A), alpha = 0.4) + facet_wrap(~Source)
}

predict_cox.aalen <- function(covar,betaHat, cum_base_haz){ #W is the rest of the covariates, A is treatment
  covar <- data.matrix(covar)
  val <- exp(covar %*% betaHat) %*% cum_base_haz
  return(val)
}

#Fit cox models----------------------
oracle_model <- function(data){
  mod <- cox.aalen(formula = Surv(T_obs, Uncensored) ~  prop(A) + prop(X) + prop(Z), data = data)
  return(mod)
}
get_oracle_covar <- function(data){
  res <- data[,names(data) %in% c("A", "X", "Z")]
  return(res)
}


non_oracle_model <- function(data){
  mod <- cox.aalen(formula = Surv(T_obs, Uncensored) ~  prop(A) + prop(fake_2) + prop(Z*A) + prop(X^2), data = data)
  return(mod)
}

get_n_oracle_covar <- function(data){
  dat <- data %>% mutate(ZA = Z*A, X2 = X^2)
  res <- dat[,names(dat) %in% c("A", "fake_2", "ZA", "X2")]
  return(res)
}

oracle_cens_model <- function(data){
  mod <- cox.aalen(Surv(get_observed_times(data), get_censor_status(data) == F) ~ 
              prop(X) + prop(Z), data = data)
}

get_oracle_cens_covar <- function(data){
  res <- data[,names(data) %in% c("X", "Z")]
  return(res)
}

non_oracle_cens_model <- function(data){
  mod <- cox.aalen(Surv(get_observed_times(data), get_censor_status(data) == F)  ~ prop(A) + prop(fake_2) + prop(Z*A) + prop(X^2), data = data)
return(mod)
}



#EIF components--------------
#Propensity
propensity <- function(data){
  mod <- glm(A ~ X + Z + fake_1 + fake_2 + fake_3, family = binomial(link = "logit"), data = data)
  
  new_data <- get_nA_covariates_val(data) %>% data.frame()
  
  predictions <- predict(object = mod, newdata = new_data) 
  
  #probabilities
  p_a_1 <- plogis(predictions)
  p_a_0 <- 1 - p_a_1
  probs <- cbind(p_a_0, p_a_1)
  colnames(probs) <- c("p_a_0", "p_a_1")
  #propensities
  pi_a_1 <- 1/p_a_1 * (get_A_val(data) == 1)
  pi_a_0 <- 1/p_a_0 * (get_A_val(data) == 0)
  propens <- list("pi_a_0" = pi_a_0, "pi_a_1" = pi_a_1)
  
  #Coefficients
  coeff <- plogis(mod$coefficients)
  
  output <- list('probs' = probs, 'propens' = propens, 'coeff' = coeff)
  return(output)
}




#Martingale
estimate_martingales <- function(model, full_data, model_cov){
  working_model <- model
  mod_cum <- get_cum(working_model)
  mod_jump_times <- get_jump_times(cum_haz_matrix = mod_cum)
  mod_cum_baseline <- get_cum_base_haz(mod_cum)
  n_jump_times <- length(mod_jump_times)
  beta_hat <- get_param_est(working_model)
  
  
  
  
  mod_est_cum_haz <- predict_cox.aalen(covar = model_cov, betaHat = beta_hat, cum_base_haz = mod_cum_baseline)

  #List of observed times
  T_obs <- get_observed_times(full_data)
  
  #Generating at risk matrix. Notice the equality!
  at_risk <- outer(X = T_obs, Y = mod_jump_times, FUN = ">=")
  
  #Change in cumulative hazard per individual (rows) per stopping time (columns).
  mod_dL <- cbind(0, mod_est_cum_haz[,-1] - mod_est_cum_haz[,-ncol(mod_est_cum_haz)])
  
  #Estimating the change in counting process for each individual i dN_i(t) = I(T_i = t, Delta = 1):
  indicator_jump_time <- outer(T_obs, mod_jump_times, FUN = "==")
  
  mod_dN <- indicator_jump_time
  
  mod_dM <- mod_dN - at_risk*mod_dL
  
  
  res <- list(mod_dM = mod_dM, mod_dN = mod_dN, mod_dL = mod_dL, jump_times = mod_jump_times, at_risk = at_risk)
  return(res)
}




#Survival function estimation
get_S_hats <- function(mod, model_cov){
  #Fitting model
  #Returning relevant figures from the model to estimate survival functions
  n <- get_row_length(model_cov)
  cum_haz_matrix <- mod$cum
  n_jumps <- get_row_length(cum_haz_matrix)
  jump_times <- get_jump_times(cum_haz_matrix)
  cum_haz <- get_cum_base_haz(cum_haz_matrix)
  beta_hat <- get_param_est(mod)
 
  covar_true <- model_cov
  
  #creating counterfactual datasets
  covar_A1 <- model_cov %>% mutate(A = 1)
  covar_A0 <- model_cov %>% mutate(A = 0)
  
  
  cumhaz_hat <- predict_cox.aalen(covar = covar_true, betaHat = beta_hat, cum_base_haz = cum_haz)
  cumhaz1_hat <- predict_cox.aalen(covar = covar_A1, betaHat = beta_hat, cum_base_haz = cum_haz)
  cumhaz0_hat <- predict_cox.aalen(covar = covar_A0, betaHat = beta_hat, cum_base_haz = cum_haz)
  
  
  #Estimating survival function
  S_hat <- exp(-cumhaz_hat)
  S_hat_1 <- exp(-cumhaz1_hat)
  S_hat_0 <- exp(-cumhaz0_hat)
  
  return(list('S_hat' = S_hat, 'S_hat_0' = S_hat_0, 'S_hat_1' = S_hat_1, 'beta_hat' = beta_hat))
}


P_treatment_extend_survival <- function(model, max_time, model_cov){
  n <- get_row_length(data)
  cum_haz_matrix <- get_cum(model)
  n_jumps <- get_row_length(cum_haz_matrix)
  jump_times <- get_jump_times(cum_haz_matrix)
  jump_times_to_keep <- sum(jump_times <= max_time)
  
  cum_haz <- get_cum_base_haz(cum_haz_matrix)
  beta_hat <- get_param_est(model)
 
  covar_true <- model_cov
  
  #creating counterfactual datasets
  covar_A1 <- model_cov %>% mutate(A = 1)
  covar_A0 <- model_cov %>% mutate(A = 0)
  
  
  cumhaz_hat <- predict_cox.aalen(covar = covar_true, betaHat = beta_hat, cum_base_haz = cum_haz)[,1:jump_times_to_keep]
  cumhaz1_hat <- predict_cox.aalen(covar = covar_A1, betaHat = beta_hat, cum_base_haz = cum_haz)[,1:jump_times_to_keep]
  cumhaz0_hat <- predict_cox.aalen(covar = covar_A0, betaHat = beta_hat, cum_base_haz = cum_haz)[,1:jump_times_to_keep]
      
  
  #Estimating survival function
  Shat_1 <- exp(-cumhaz1_hat)
  Shat_0 <- exp(-cumhaz0_hat)
  
  
  delta_cumbasehaz <- c(0, cum_haz[-1] - cum_haz[-n_jumps])[1:jump_times_to_keep]
  
  integrand <- t(t(Shat_1 * Shat_0) * delta_cumbasehaz)
  multiplier <- exp(data.matrix(covar_A0) %*% beta_hat)
  
  res <- multiplier * rowSums(integrand)
  
  res_list <- list('res' = res, 'delta_cumbasehaz' = delta_cumbasehaz, 'integrand' = integrand, 'multiplier' = multiplier)
  
  return(res_list)
}


#Censoring distribution estimation
get_K_hat <- function(data, T_model, corr_cens_model = T){
  
  n <- get_row_length(data)
  
  ######################################################################
  #                         Getting surv jump times
  ######################################################################
  
  cum_obs <- get_cum(T_model)
  no_jumps_obs <- get_row_length(cum_obs)                    # Number of jumps T > t
  tau_obs <- get_jump_times(cum_obs)                          # Jump times tau1, tau2,...,tau_no_jumps
  
  ######################################################################
  
  
  ######################################################################
  #                         Getting censoring distribution
  ######################################################################
  #Fitting model for censoring distribution
  if(corr_cens_model){
    mod_cens <- oracle_cens_model(data)
    covar <- get_oracle_cens_covar(data)
  } else {
    mod_cens <- non_oracle_cens_model(data)
    covar <- get_n_oracle_covar(data)
  }
  
  #Jump times and cumulative baseline hazard:
  cum_cens <- get_cum(mod_cens)
  no_jumps_cens <- get_row_length(cum_cens)                   #Number of jumps T > t
  tau_cens <- get_jump_times(cum_cens)                        #Jump times tau_cens1, tau_cens2,...,tau_cens_no_jumps
  cum_basehaz_cens <- get_cum_base_haz(cum_cens)              #Cumulative baseline hazard   Lambda0(tau0), Lambda0(tau1),...,Lambda0(tau_no_jumps)
  
  #Betahat
  betahat_cens <- get_param_est(mod_cens)
  
  #Estimated cumulative hazard assuming the cox-model
  cumhaz_cens_hat <- predict_cox.aalen(covar = covar, betaHat = betahat_cens, cum_base_haz = cum_basehaz_cens)
  
  #Estimating censoring function K_C(t | A, X)
  Khat_temp <- exp(-cumhaz_cens_hat)
  
  #The problem now is that the censoring times and failure times do not jump at the same time and we wish to look at the jump times
  #of the survival function. That is we wish to return Khat not in the jump times tau_cens but in tau. Since Khat is constant in between
  #jump times, we can simply take the value of Khat corresponding to the largest value of tau_cens that is lower than tau
  
  
  Khat = matrix(NA, nrow = n, ncol = no_jumps_obs) #Maybe this can also be vectorized??
  for(i in 1:no_jumps_obs){
    Khat[,i] = Khat_temp[,max((1:no_jumps_cens)[tau_cens<=tau_obs[i]])]
  }
  
  return(list(Khat = Khat, mod_cens = mod_cens))
}




data <- generate_survival_data(100)
T_corr = T
Cens_corr = T
max_time = 5
EIF <- function(data, T_corr = T, Cens_corr = T, max_time = 5){
  #Building models of T and Cens
  if(T_corr){
    T_model <- oracle_model(data)
    T_model_cov <- get_oracle_covar(data)}
  else{
    T_model <- non_oracle_model(data)
    T_model_cov <- get_n_oracle_covar(data)}
  
  
  #martingales
  martingale_estimates <- estimate_martingales(T_model,  full_data = data, model_cov = T_model_cov)
  #martingale_estimates <- estimate_martingales_old(data)
  #mean(martingale_estimates$mod_dL == martingale_estimates_1$mod_d)
  #View(martingale_estimates_1)
  jump_times <- martingale_estimates$jump_times
  jump_times_to_keep <- sum(jump_times <= max_time)
  
  #Only sum op to tau
  dM <- get_mg_dM(martingale_estimates)[,1:jump_times_to_keep]
  dL <- get_mg_dL(martingale_estimates)[,1:jump_times_to_keep]
  
  #Survival functions
  S_hats <- get_S_hats(mod = T_model, model_cov = T_model_cov)
  S_hat <- S_hats$S_hat[,1:jump_times_to_keep]
  S_hat0 <- S_hats$S_hat_0[,1:jump_times_to_keep]
  S_hat1 <- S_hats$S_hat_1[,1:jump_times_to_keep]
  betahat <- S_hats$beta_hat
  
  
  #Censoring function estimate
  K_C_hat <- get_K_hat(data = data, T_model = T_model, corr_cens_model = Cens_corr)$Khat[,1:jump_times_to_keep]
  
  
  #Propensity scores
  prop <- propensity(data)$propens
  prop_0 <- prop$pi_a_0
  prop_1 <- prop$pi_a_1
  prop_sum <- prop_0 + prop_1
  
  
  #Her mangler vi at fjerne rækker efter tau!!!
  #Part 1
  p1 <- prop_0 * rowSums((S_hat1/K_C_hat)*dM)
  
  
  #Part 2
  p_treat_list <- P_treatment_extend_survival(model = T_model,model_cov = T_model_cov, max_time = max_time)
  p2 <- p_treat_list$res
  
  #Her mangler vi at fjerne rækker efter tau!!!
  #Part 3
  inner_integrand <- dM/(S_hat * K_C_hat)
  inner_integral <- t(apply(inner_integrand, MARGIN = 1, FUN = cumsum))
  
  p3 <- - prop_sum * p_treat_list$multiplier * rowSums(p_treat_list$integrand * inner_integral)
  
  output <- list(p1, p2, p3, est = mean(p1+p2+p3))
  return(output)
}

check_EIF <- function(n, T_corr = T, Cens_corr = T, max_time = 5, ba_t = -1){
  dat <- generate_survival_data(n, ba_t = ba_t)
  proportion_observed(dat)
  visualize_data(dat)
  return(EIF(dat, T_corr = T_corr, Cens_corr = Cens_corr, max_time = max_time))
}

#check_EIF(1000, ba_t = 1)

