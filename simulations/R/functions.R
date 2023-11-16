####Packages----
library(targets)
library(tidyverse)
library(survival)
library(timereg)
library(dplyr)
library(data.table)




#Get functions----------------------
#Optimized for generation method----




get_X_values <- function(data){
  return(data$X)
}

get_A_values <- function(data){
  return(data$A)
}
get_T_values <- function(data){
  return(data$T_true)
}

get_row_length <- function(data){
  return(nrow(data))
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

get_A_estimate <- function(model){
  return(model$coefficients["A"])
}

get_X_estimate <- function(model){
  return(model$coefficients["X"])
}

get_se_A_estimate <- function(model){
  return(sqrt(model$var[1,1]))
}

get_se_X_estimate <- function(model){
  return(sqrt(model$var[2,2]))
}

get_jump_times <- function(cum_hazard_matrix){
  return(cum_hazard_matrix[,"time"])
}

get_cum_baseline_hazard_times <- function(cum_hazard_matrix){
  return(cum_hazard_matrix[,"(Intercept)"])
}

get_parameter_estimates <- function(cox.aalen_model){
  return(cox.aalen_model$gamma)
}

#Filtering functions ---------------
filter_on_a <- function(data, val){
  return(data %>% filter(a == val))
}

filter_on_x <- function(data, val){
  return(data %>% filter(x == val))
}



#GENERATING METHODS-----------------


generate_survival_data <- function(n, ba_t = log(2), bx_t = log(2), bz_t = log(2),
                                   ba_c = log(2), bx_c = log(2), bz_c = log(2),
                                   surv_is_cox = T, cens_is_cox = T, 
                                   prop_a = 1/2, 
                                   x_range = c(-1,1), x_cont = T, x_vals = c(0,1)){
  #Generate list of treatment indicators
  
  
  z <- runif(n = n, min = 0, max = 1)
  
  if(!surv_is_cox){
    prop_a <- z
  }
  
  a <- rbinom(n, 1, prop_a) 
  
  if(x_cont){
    x <- runif(n, min = min(x_range), max = max(x_range))
  }
  else{
    x <- sample(x_vals, n, replace = T)
  }
  
  #Generating survival times
  if(surv_is_cox){denom_surv <- exp(ba_t*a + bx_t*x)} 
    else {denom_surv <- exp(ba_t*a + bx_t*(x^2)+bz_t*z)}
  
  surv_t <- rexp(n, 1)/denom_surv
  
  #Generating censoring times
  if(cens_is_cox){denom_cens <- exp(bx_c*x)} 
    else {denom_cens <- exp(bx_c*x^2+bz_c*z)}
  
  cens_t <- rexp(n, 1)/denom_cens
  #Observed values
  t_obs <- pmin(surv_t, cens_t)
  
  #Status
  delta <- surv_t < cens_t
  
  return(data.frame(T_true = surv_t, C = cens_t, T_obs = t_obs, Uncensored = delta, X = x, A = a, Z = z))
}




#Summary statistics ----------------

proportion_censored <- function(data){
  mean(get_censor_status(data))
}



#Cox models-------------------------

#Fit Cox-model
cox_naive_model <- function(data){
  coxph(formula = Surv(T_obs, Uncensored) ~  A + X, data = data)
}

cox.aalen_naive_model <- function(data){
  cox.aalen(formula = Surv(T_obs, Uncensored) ~  prop(A) + prop(X), data = data)
}

#Fit Oracle Cox_model
cox_oracle_model <- function(data, surv_is_cox = T){
  if(surv_is_cox){
    mod <- cox_naive_model(data)
  }
  else{
    mod <- coxph(formula = Surv(T_obs, Uncensored) ~  A + X2 + Z, data = data)
  }
  return(mod)
}


#Helper functions------------------


#Empirical distribution
empirical_dist <- function(data, a, x){
  
  filtered_data <- data %>% filter_on_a(a) %>% filter_on_x(x)
  
  cond_empirical_dist <- function(t){
    
    cal <- function(t){
      return(mean(get_observed_times(filtered_data) > t))
    }
    
    output <- apply(matrix(t), 1, FUN = cal)
    
    return(output)
  }
  
  return(cond_empirical_dist)
}


#Predict based on estimates
predict_cox.aalen <- function(A, W, betaHat, cum_base_haz, n_jumps){ #W is the rest of the covariates, A is treatment
  val <- (exp(cbind(A, W) %*% betaHat) %*% cum_base_haz)[,2:(n_jumps)] #hvorfor skubbe med 1?
  return(val)
}
cbind(A,X)%*%beta_hat%*%cum_haz

#Confidence intervals---------------

#Create CI for A based on Cox-regression
confidence_interval_A_cox <- function(model, conf_level = 0.05){
  CI <- confint(model, parm = "A", conf_level = conf_level)
  return(CI)
}

#Create CI for X based on Cox regression
confidence_interval_X_cox <- function(model, conf_level = 0.05){
  CI <- confint(model, parm = "X", conf_level = conf_level)
  return(CI)
}


# A function for checking if the model correctly includes the true value in the CI
check_val_in_CI <- function(val, CI){
return(between(x = val, lower = min(CI), upper = max(CI)))
}


check_estimate_A <- function(model, ba_t){
  return(check_val_in_CI(ba_t, confidence_interval_A_cox(model)))
}

check_estimate_X <- function(model, bx_t){
  return(check_val_in_CI(bx_t, confidence_interval_X_cox(model = model)))
}


#Troubleshooting code--------------

#A function that generates data with n observations and variable if Cox-generated 
#RETURNS: T/F whether model included the true value of treatment effect
check_mod <- function(n = 400, surv_is_cox = F, treatment_effect =  log(2), cens_is_cox = T, x_vals = c(0,1)){
  #Generate data
  dat <- generate_survival_data(n = n, ba_t = treatment_effect, surv_is_cox = surv_is_cox, ba_c = 1, cens_is_cox = cens_is_cox, x_vals = x_vals)
  #Fit Cox model
  mod <- cox_naive_model(dat)
  return(check_estimate_A(mod, treatment_effect))
}


#Create a model comparison function
#INPUT: Data (df), surv_is_cox (Boolean) - decides how Oracle is defined, 
#OUTPUT: Boolean on whether the CI included true value
model_comparison <- function(data, surv_is_cox, treatment_effect){
  #If the data has not been generated as a Cox-model, we need to find X^2
  if(!surv_is_cox){
    data$X2 <- get_X_values(data)^2
  }
  
  #Fit naive model
  mod_naive <- cox_naive_model(data)
  #Fit Oracle model Should never perform worse than naive model
  mod_oracle <- cox_oracle_model(data, surv_is_cox)
  
  #Checking if the models have included the true value of ba
  a_naive <- get_A_estimate(mod_naive)
  a_oracle <- get_A_estimate(mod_oracle)
  naive_correct <- check_estimate_A(model = mod_naive, ba_t = treatment_effect)
  oracle_correct <- check_estimate_A(model = mod_oracle,ba_t = treatment_effect )
  
  output <- c(naive_correct, oracle_correct, a_naive, a_oracle)
  names(output ) <- c("naive", "oracle", "naive", "oracle")
  return(output)
}


#INPUT: n - number of data points, Treatment_effect, x_vals - the range for the x_values
#Creating function that showcases the problem. Assume that the treatment effect is 0
showcase_oracle_not_better <- function(n = 400, treatment_effect = log(2), 
                                       x_range = c(0,4), surv_is_cox = F, cens_is_cox =F){
  dat <- generate_survival_data(n, ba_t = treatment_effect,  
                                surv_is_cox = surv_is_cox, cens_is_cox = cens_is_cox,
                                x_range = x_range)
  res <- model_comparison(dat, surv_is_cox = surv_is_cox, treatment_effect = treatment_effect)
  cens <- proportion_censored(dat)
  names(cens) <- c("prop not-censored")
  return(c(res, cens))
}
#showcase_oracle_not_better()
#
#showcase_oracle_not_better(treatment_effect = )
#res <- replicate(1000, showcase_oracle_not_better(n = 1000, treatment_effect = log(2), x_range = c(0,3), surv_is_cox = F, cens_is_cox = T))




#PROBLEMS-----------------------------
#The code below is programmed for reproducibility

###
#COX MODEL HAS TOO HIGH PRECISION
###


#Do the experiment l times
#l <- 1000
#res1 <- replicate(l, check_mod(x_vals = c(0,5), surv_is_cox = T))
##Find the average number of times the interval included the true effect.
##Expected: 95%
##Reason: We construct a 95% CI around the point estimate, and the function is correctly specified
#sum(res1);mean(res1)
##Result: 100% inclusion <- something is wrong in the procedure
#
##Do the same for non-cox generated data
#res2 <- replicate(l, check_mod(surv_is_cox = F))
##Find the average number of times the interval included the true effect.
##Expected: less than 95%
##Reason: The data has not been generated as a Cox_model, so we are using the wrong model to analyze the data
#sum(res2);mean(res2)


###
#ORACLE MODEL DOES NOT PERFORM AS GOOD AS EXPECTED
###


#res <- replicate(n = 1000, check_mod())
#res <- replicate(n = 10000, showcase_oracle_not_better(n = 1000))

#rowMeans(res)
#
#
#
##replicating study l times
#l <- 1000
##result <- replicate(l, showcase_problem(surv_is_cox = T))
#result_2 <- replicate(l, showcase_oracle_not_better(surv_is_cox = F, x_range = c(0,4)))
####Expected: Oracle should perform better than normal cox
####Reason: data has not been created in Cox-manner
##rowSums(result);rowMeans(result)
#rowSums(result_2); rowMeans(result_2)
##result: The oracle and the cox model both predict correctly with 100% accuracy, which is very problemati
#
#save_data <- matrix(c(0,0,0,0), nrow = 1, ncol = 4)
#save_data <- rbind(save_data, c(0,1,2))
#save_data
#for (i in (1:10)){
#  res <- replicate(l, showcase_oracle_not_better(surv_is_cox = F, x_range = c(0,i)))
#  results <- c(i, rowMeans(res))
#  save_data <- rbind(save_data, results)
#}
#




#m <- cox.aalen(formula = Surv(t_observed, has_been_censored) ~ prop(a)+prop(x), data = dat)
#plot(m)
#
#?coxph





########Plots for checking if stuff works
#med_test <- median(get_observed_times(dat))
#f <- empirical_dist(dat, 1, 0)
#g <- empirical_dist(dat, 1, 1)
#h <- empirical_dist(dat, 1, 2)
#f(med_test)

#plot(f, from = 0, to = 1)
#plot(g, from = 0, to = 1, add = T, col = "red")
#plot(h, from = 0, to = 1, add = T, col = "green")


#EIF components----------------------


#Jeg forsøgt at implementere det, men der nogle problmer med dimensioner, som jeg ikke helt kan forstå
#Estimating P(T1 > T0 | W)
P_treatment_extend_survival <- function(data){
  mod <- cox.aalen_naive_model(data)
  n <- get_row_length(data)
  cum_haz_matrix <- mod1$cum
  n_jumps <- get_row_length(cum_haz_matrix)
  jump_times <- get_jump_times(cum_haz_matrix)
  cum_haz <- get_cum_baseline_hazard_times(cum_haz_matrix)
  beta_hat <- get_parameter_estimates(mod)
  s
  A <- get_A_values(data)
  X <- get_X_values(data)
  #Z <- get_Z_values(data)
  cumhaz_hat <- predict_cox.aalen(A = A, W = X, betaHat = beta_hat, cum_base_haz = cum_haz, n_jumps = n_jumps)
  cumhaz1_hat <- predict_cox.aalen(A = 1, W = X, betaHat = beta_hat, cum_base_haz = cum_haz, n_jumps = n_jumps)
  cumhaz0_hat <- predict_cox.aalen(A = 0, W = X, betaHat = beta_hat, cum_base_haz = cum_haz, n_jumps = n_jumps)
  
  #Estimating survival function
  Shat_1 <- exp(-cumhaz1_hat)
  Shat_0 <- exp(-cumhaz0_hat)
  
  
  delta_cumbasehaz <- matrix(data = cum_haz[-1], nrow = n, ncol = n_jumps-1, byrow = T) -
    matrix(data = cum_haz[-(n_jumps)], nrow = n, ncol = n_jumps-1, byrow = T)
  
  delta_cumbasehaz2 <- c(delta_cumbasehaz, delta_cumbasehaz[,230])
  
  
  integrand <- Shat_1 * Shat_0 * delta_cumbasehaz2
  
  return(mean(exp(beta_hat["prop(X)"] * X) * rowSums(integrand))
)
  
}

data <- generate_survival_data(100)
P_treatment_extend_survival(dat)

#----------------------------------------------
###########################
#LEGACY CODE
#Redundant methods made obsolete by Martinussen's method
###########################
#-------------------------------

##############
#Get functions
##############
#get_uni_a0 <- function(uni_data){
#  return(uni_data$A0)
#}
#get_uni_a1 <- function(uni_data){
#  return(uni_data$A1)
#}


###################
#Quantile functions
###################

#Assume that baseline hazard is 1

##Cox-generated survival times
#quantile_surv_generate_cox <- function(q,b1 = 1,b2 = 1, a, x){                 #quantile function for h0
#  (-log(1-q))/(exp(b1*a+b2*x))
#}
#
##Not Cox-generated survival times
#quantile_surv_generate_not_cox <- function(q,b1 = 1,b2 = 1, a, x){                 #quantile function for h0
#  (-log(1-q))/(exp(b1*a+b2*x^2))
#}
#
##Censoring times
#conditional_censoring_cox <- function(q,b1 = 1,b2 = 1, a, x){                 #quantile function for h0
#  (-log(1-q))/(exp(b1*a+b2*x))
#}
#  
#
#conditional_censoring_not_cox <- function(q,b1 = 1,b2 = 1, a, x){                 #quantile function for h0
#  (-log(1-q))/(exp(b1*a+b2*x^2))
#}
#

#####################
#Generating methods
#####################
#generate_random <- function(n, prop = 1/2){
#  uniformT <- runif(n, 0, 1)
#  uniformT0 <- uniformT[1:(floor(n*prop))]                    #data that follows cox-model
#  uniformT1 <- uniformT[(floor(n*prop+1)):n]                  #data that does not follow cox-model
#  return(list(A0 = uniformT0, A1 = uniformT1))
#}



#generate_survival_times <- function(n, x_vals, b1 = 1, b2 = 1, is_cox = T){
#  len_x <- length(x_vals)
#  n_corrected <- n-n%%len_x
#  data_frame_storage <- data.frame()
#  
#  for(x in x_vals){
#    randomness <- generate_random(n_corrected/len_x)
#    q_a0 <- get_uni_a0(randomness)
#    q_a1 <- get_uni_a1(randomness)
#    
#    if(is_cox){
#      a0_times <- quantile_surv_generate_cox(q_a0, 
#                                             b1 = b1, 
#                                             b2 = b2, 
#                                             a = 0, 
#                                             x = x)
#      
#      a1_times <- quantile_surv_generate_cox(q_a1, 
#                                             b1 = b1, 
#                                             b2 = b2, 
#                                             a = 1, 
#                                             x = x)
#      
#    }
#    else{
#      a0_times <- quantile_surv_generate_not_cox(q_a0, 
#                                                 b1 = b1, 
#                                                 b2 = b2, 
#                                                 a = 0, 
#                                                 x = x)
#      a1_times <- quantile_surv_generate_not_cox(q_a1, 
#                                                 b1 = b1, 
#                                                 b2 = b2, 
#                                                 a = 1, 
#                                                 x = x)
#    }
#    
#    temp_df0 <- data.frame(t_true = a0_times, a = rep(0, length(a0_times)), x = x)
#    temp_df1 <- data.frame(t_true = a1_times, a = rep(1, length(a1_times)), x = x)
#    data_frame_storage <- rbind(data_frame_storage, temp_df1, temp_df0)
#  }
#  
#  
#  
#  return(data_frame_storage)
#}
#
#generate_censoring_times <- function(n, x_vals, b1 = 1, b2 = 1, is_cox = T){
#  len_x <- length(x_vals)
#  n_corrected <- n-n%%len_x
#  data_frame_storage <- data.frame()
#  
#  for(x in x_vals){
#    randomness <- generate_random(n_corrected/len_x)
#    q_a0 <- get_uni_a0(randomness)
#    q_a1 <- get_uni_a1(randomness)
#    
#    if(is_cox){
#      a0_times <- conditional_censoring_cox(q_a0, 
#                                            b1 = b1, 
#                                            b2 = b2, 
#                                            a = 0, 
#                                            x = x)
#      
#      a1_times <- conditional_censoring_cox(q_a1, 
#                                            b1 = b1, 
#                                            b2 = b2, 
#                                            a = 1, 
#                                            x = x)
#    }
#    else{
#      a0_times <- conditional_censoring_not_cox(q_a0, 
#                                                b1 = b1, 
#                                                b2 = b2, 
#                                                a = 0, 
#                                                x = x)
#      a1_times <- conditional_censoring_not_cox(q_a1, 
#                                                b1 = b1, 
#                                                b2 = b2, 
#                                                a = 1, 
#                                                x = x)
#    }
#    
#    temp_df0 <- data.frame(t_true = a0_times, a = rep(0, length(a0_times)), x = x)
#    temp_df1 <- data.frame(t_true = a1_times, a = rep(1, length(a1_times)), x = x)
#    data_frame_storage <- rbind(data_frame_storage, temp_df1, temp_df0)
#  }
#  
#  return(data_frame_storage)
#}
#
#
#generate_survival_data <- function(n, x_vals, b1 = 1, b2 = 1, surv_is_cox = T, cens_is_cox = T){
#  survival_times <- generate_survival_times(n = n,b1 = b1, b2 = b2, x_vals = x_vals, is_cox = surv_is_cox)
#  censoring_times <- get_t_values(generate_censoring_times(n = n,b1 = b1, b2 = b2, x_vals = x_vals, is_cox = cens_is_cox))
#  has_not_been_censored <- get_t_values(survival_times) < censoring_times
#  t_observed <- pmin(get_t_values(survival_times), censoring_times)
#  full_data_set <- cbind(survival_times, censoring_times, t_observed, has_not_been_censored)
#  return(full_data_set)
#}
#

# CODE FOR COMPARING MODELS
#
#sim_data <- function(n, ba_t = 1, bx_t = 1, ba_c = 1, bx_c = 1, surv_is_cox = T, 
#                     cens_is_cox = T, prop_a = 1/2, 
#                     x_range = c(-1,1), x_cont = T, 
#                     x_vals = c(0,1)){
#  dat <- generate_survival_data(n, ba_t, )
#  if(!surv_is_cox){
#    dat$X2 <- get_X_values(dat)^2
#  }
#  mod_naive <- cox_naive_model(dat)
#  mod_oracle <- cox_oracle_model(dat, surv_is_cox)
#  print(dat)
#  print(proportion_censored(dat))
#  naive_correct <- check_estimate_A(model = mod_naive, ba_t = ba_t)
#  oracle_correct <- check_estimate_A(model = mod_oracle,ba_t = ba_t )
#  print(confidence_interval_A_cox(model = mod_naive))
#  print(confidence_interval_A_cox(model = mod_oracle))
#  cat("True value is ", ba_t, "\n")
#  return(c(naive_correct, oracle_correct))
#  
#}

##########################
#LEGACY CONFIDENCE INTERVAL
###########################

#Create CI for A based on Cox-regression
#confidence_interval_A_cox <- function(model, conf_level = 0.05){
   #se_A <- get_se_A_estimate(model)
  #A <- get_A_estimate(model)
  #q <- qnorm(conf_level/2)
  #lower <- A + q*se_A
  #upper <- A - q*se_A
  #CI <- c(lower, upper)
  #names(CI) <- c("lower", "upper")
#  return(CI)
#}

##Create CI for X based on Cox regression
#confidence_interval_X_cox <- function(model, conf_level = 0.05){
#  se_X <- get_se_X_estimate(model)
#  X <- get_X_estimate(model)
#  q <- qnorm(conf_level/2)
#  lower <- X + q*se_X
#  upper <- X - q*se_X
#  CI <- c(lower, upper)
#  names(CI) <- c("lower", "upper")
#  return(CI)
#}
