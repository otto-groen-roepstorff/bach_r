###Packages
#library(targets)
#library(tidyverse)
#library(survival)
#library(timereg)
#library(dplyr)
#library(data.table)



##############
#Get functions
#Optimized for generation method
##############


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
  return(sqrt(mod$var[1,1]))
}

get_se_X_estimate <- function(model){
  return(sqrt(mod$var[2,2]))
}

####################
#Filtering functions
####################

filter_on_a <- function(data, val){
  return(data %>% filter(a == val))
}

filter_on_x <- function(data, val){
  return(data %>% filter(x == val))
}


#####################
#GENERATING METHODS
#####################


generate_survival_data <- function(n, ba_t = 1, bx_t = 1, ba_c = 1, bx_c = 1, surv_is_cox = T, 
                                        cens_is_cox = T, prop_a = 1/2, 
                                        x_range = c(-1,1), x_cont = T, 
                                        x_vals = c(0,1)){
  #Generate list of treatment indicators
  a <- rbinom(n, 1, prop_a) 
  
  if(x_cont){
    x <- runif(n, min = min(x_range), max = max(x_range))
  }
  else{
    x <- sample(x_vals, n, replace = T)
  }
  #Generating survival times
  if(surv_is_cox){
    denom_surv <- exp(ba_t*a + bx_t*x)
  } else {
    denom_surv <- exp(ba_t*a + bx_t*(x^2))
  }
  surv_t <- rexp(n, 1)/denom_surv
  
  #Generating censoring times
  if(cens_is_cox){
    denom_cens <- exp(ba_c*a + bx_c*x)
  } else {
    denom_cens <- exp(ba_c*a + bx_c*x^2)
  }
  cens_t <- rexp(n, 1)/denom_cens
  #Observed values
  t_obs <- pmin(surv_t, cens_t)
  
  #Status
  delta <- surv_t < cens_t
  
  return(data.frame(T_true = surv_t, C = cens_t, T_obs = t_obs, Uncensored = delta, X = x, A = a))
}
  
  
###################
#Summary statistics
###################

proportion_censored <- function(data){
  mean(get_censor_status(data))
}

#Run Cox-model
cox_naive_model <- function(data){
  coxph(formula = Surv(T_obs, Uncensored) ~  A + X, data = data)
}

#Run Oracle Cox_model

#################
#Helper functions
#################

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


################################
#Confidence intervals
#####################
#---------------------


#Create CI for A based on Cox-regression
confidence_interval_A_cox <- function(model, conf_level = 0.05){
  se_A <- get_se_A_estimate(model)
  A <- get_A_estimate(model)
  q <- qnorm(conf_level/2)
  lower <- A + q*se_A
  upper <- A - q*se_A
  CI <- c(lower, upper)
  names(CI) <- c("lower", "upper")
  return(CI)
}

#Create CI for X based on Cox regression
confidence_interval_X_cox <- function(model, conf_level = 0.05){
  se_X <- get_se_X_estimate(model)
  X <- get_X_estimate(model)
  q <- qnorm(conf_level/2)
  lower <- X + q*se_X
  upper <- X - q*se_X
  CI <- c(lower, upper)
  names(CI) <- c("lower", "upper")
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

#############################################
#PROBLEMS
#---------------
#The code below is programmed for reproducibility
#############################

#################################
#COX MODEL HAS TOO HIGH PRECISION
#################################

#A function that generates data with n observations and variable if Cox-generated 
#RETURNS: T/F whether model included the true value of treatment effect

check_mod <- function(n = 10000, surv_is_cox = T, treatment_effect =  0, cens_is_cox = T){
      #Generate data
      dat <- generate_survival_data(n = n, ba_t = treatment_effect, surv_is_cox = surv_is_cox, ba_c = 1, cens_is_cox = cens_is_cox)
      #Fit Cox model
      mod <- cox_naive_model(dat)
      return(check_estimate_A(mod, 0))
}

#Do the experiment l times
l <- 1000
res1 <- replicate(l, check_mod(cens_is_cox = F))
#Find the average number of times the interval included the true effect.
#Expected: 95%
#Reason: We construct a 95% CI around the point estimate, and the function is correctly specified
sum(res1);mean(res1)
#Result: 100% inclusion <- something is wrong in the procedure

#Do the same for non-cox generated data
res2 <- replicate(l, check_mod(surv_is_cox = F))
#Find the average number of times the interval included the true effect.
#Expected: less than 95%
#Reason: The data has not been generated as a Cox_model, so we are using the wrong model to analyze the data
sum(res2);mean(res2)


##################################################
#ORACLE MODEL DOES NOT PERFORM AS GOOD AS EXPECTED
##################################################

#Define Oracle model
cox_oracle_model <- function(data, surv_is_cox = T){
  if(surv_is_cox){
    mod <- coxph(formula = Surv(T_obs, Uncensored) ~  A + X, data = data)
  }
  else{
    mod <- coxph(formula = Surv(T_obs, Uncensored) ~  A + X2, data = data)
    }
  return(mod)
}

#Create a model comparison function
model_comparison <- function(data, surv_is_cox, ba_T){
  #If the data has not been generated as a Cox-model, we need to find X^2
  if(!surv_is_cox){
    data$X2 <- get_X_values(data)^2
  }
  
  #Fit naive model
  mod_naive <- cox_naive_model(data)
  #Fit Oracle model Should never perform worse than naive model
  mod_oracle <- cox_oracle_model(data, surv_is_cox)
  
  #Checking if the models have included the true value of ba
  naive_correct <- check_estimate_A(model = mod_naive, ba_t = ba_T)
  oracle_correct <- check_estimate_A(model = mod_oracle,ba_t = ba_T )
  
  output <- c(naive_correct, oracle_correct)
  names(output ) <- c("naive", "oracle")
  #returning a tuple with Boolean if the models guessed correctly
  return(output)
}

#Creating function that showcases the problem. Assume that the treatment effect is 0
showcase_problem <- function(n = 10000, treatment_effect = 0){
  dat <- generate_survival_data(n, surv_is_cox = F, ba_t = treatment_effect)
  return(model_comparison(dat, surv_is_cox = F, ba_T = treatment_effect))
}

#replicating study l times
l <- 1000
result <- replicate(l, showcase_problem())
#Expected: Oracle should perform better than normal cox
#Reason: data has not been created in Cox-manner
rowSums(result)
rowMeans(result)
#result: The oracle and the cox model both predict correctly with 100% accuracy, which is very problemati






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

###################
#EIF components
###############

#extended_survival_prop <- function()


#generate_censoring_times(10,x_vals = 1, is_cox = F)

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

