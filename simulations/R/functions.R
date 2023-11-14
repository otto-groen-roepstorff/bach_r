###Packages
library(tidyverse)
library(survival)
library(timereg)



#####################
#Get functions
##############

get_uni_a0 <- function(uni_data){
  return(uni_data$A0)
}
get_uni_a1 <- function(uni_data){
  return(uni_data$A1)
}

get_x_values <- function(data){
  return(data$x)
}

get_a_values <- function(data){
  return(data$a)
}
get_t_values <- function(data){
  return(data$t_true)
}

get_row_length <- function(data){
  return(nrow(data))
}

get_censor_status <- function(data){
  data$has_not_been_censored
}

get_observed_times <- function(data){
  data$t_observed
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




###################
#Quantile functions
###################

#Assume that baseline hazard is 1

#########################################################
#       What to do with the at-risk indicator?
#########################################################

#Cox-generated survival times
quantile_surv_generate_cox <- function(q,b1 = 1,b2 = 1, a, x){                 #quantile function for h0
  (-log(1-q))/(exp(b1*a+b2*x))
}

#Not Cox-generated survival times
quantile_surv_generate_not_cox <- function(q,b1 = 1,b2 = 1, a, x){                 #quantile function for h0
  (-log(1-q))/(exp(b1*a+b2*x^2))
}

#Censoring times
conditional_censoring_cox <- function(q,b1 = 1,b2 = 1, a, x){                 #quantile function for h0
  (-log(1-q))/(exp(b1*a+b2*x))
}
  

conditional_censoring_not_cox <- function(q,b1 = 1,b2 = 1, a, x){                 #quantile function for h0
  (-log(1-q))/(exp(b1*a+b2*x^2))
}

#####################
#GENERATING METHODS
#####################

generate_random <- function(n, prop = 1/2){
  uniformT <- runif(n, 0, 1)
  uniformT0 <- uniformT[1:(floor(n*prop))]                    #data that follows cox-model
  uniformT1 <- uniformT[(floor(n*prop+1)):n]                  #data that does not follow cox-model
  return(list(A0 = uniformT0, A1 = uniformT1))
}



generate_survival_times <- function(n, x_vals, b1 = 1, b2 = 1, is_cox = T){
  len_x <- length(x_vals)
  n_corrected <- n-n%%len_x
  data_frame_storage <- data.frame()
  
  for(x in x_vals){
    randomness <- generate_random(n_corrected/len_x)
    q_a0 <- get_uni_a0(randomness)
    q_a1 <- get_uni_a1(randomness)
    
    if(is_cox){
      a0_times <- quantile_surv_generate_cox(q_a0, 
                                         b1 = b1, 
                                         b2 = b2, 
                                         a = 0, 
                                         x = x)
      
      a1_times <- quantile_surv_generate_cox(q_a1, 
                                         b1 = b1, 
                                         b2 = b2, 
                                         a = 1, 
                                         x = x)
                                         
    }
    else{
      a0_times <- quantile_surv_generate_not_cox(q_a0, 
                                             b1 = b1, 
                                             b2 = b2, 
                                             a = 0, 
                                             x = x)
      a1_times <- quantile_surv_generate_not_cox(q_a1, 
                                             b1 = b1, 
                                             b2 = b2, 
                                             a = 1, 
                                             x = x)
    }
  
    temp_df0 <- data.frame(t_true = a0_times, a = rep(0, length(a0_times)), x = x)
    temp_df1 <- data.frame(t_true = a1_times, a = rep(1, length(a1_times)), x = x)
    data_frame_storage <- rbind(data_frame_storage, temp_df1, temp_df0)
  }
  
    
  
  return(data_frame_storage)
}

generate_censoring_times <- function(n, x_vals, b1 = 1, b2 = 1, is_cox = T){
  len_x <- length(x_vals)
  n_corrected <- n-n%%len_x
  data_frame_storage <- data.frame()
  
  for(x in x_vals){
    randomness <- generate_random(n_corrected/len_x)
    q_a0 <- get_uni_a0(randomness)
    q_a1 <- get_uni_a1(randomness)
    
    if(is_cox){
      a0_times <- conditional_censoring_cox(q_a0, 
                                             b1 = b1, 
                                             b2 = b2, 
                                             a = 0, 
                                             x = x)
      
      a1_times <- conditional_censoring_cox(q_a1, 
                                             b1 = b1, 
                                             b2 = b2, 
                                             a = 1, 
                                             x = x)
    }
    else{
      a0_times <- conditional_censoring_not_cox(q_a0, 
                                                 b1 = b1, 
                                                 b2 = b2, 
                                                 a = 0, 
                                                 x = x)
      a1_times <- conditional_censoring_not_cox(q_a1, 
                                                 b1 = b1, 
                                                 b2 = b2, 
                                                 a = 1, 
                                                 x = x)
    }
    
    temp_df0 <- data.frame(t_true = a0_times, a = rep(0, length(a0_times)), x = x)
    temp_df1 <- data.frame(t_true = a1_times, a = rep(1, length(a1_times)), x = x)
    data_frame_storage <- rbind(data_frame_storage, temp_df1, temp_df0)
  }
  
  return(data_frame_storage)
}


generate_survival_data <- function(n, x_vals, b1 = 1, b2 = 1, surv_is_cox = T, cens_is_cox = T){
  survival_times <- generate_survival_times(n = n,b1 = b1, b2 = b2, x_vals = x_vals, is_cox = surv_is_cox)
  censoring_times <- get_t_values(generate_censoring_times(n = n,b1 = b1, b2 = b2, x_vals = x_vals, is_cox = cens_is_cox))
    has_not_been_censored <- get_t_values(survival_times) < censoring_times
  t_observed <- pmin(get_t_values(survival_times), censoring_times)
  full_data_set <- cbind(survival_times, censoring_times, t_observed, has_not_been_censored)
  return(full_data_set)
}

surv <- generate_survival_times(100, x_vals = c(0,1,2))
generate_survival_data(100, c(0,1,2))

cens <- generate_survival_times(100, x_vals = c(0,1,2))
cbind(surv, cens)
mean(get_t_values(surv) > get_t_values(cens))

###################
#Summary statistics
###################

proportion_censored <- function(data){
  mean(get_censor_status(data))
}




#################
#Helper functions
#################


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


cox_naive_model <- function(data){
  coxph(formula = Surv(t_observed, has_not_been_censored) ~  a + x, data = data)
}

df1 <- generate_censoring_times(1000, b1 = 1, b2 = 2, is_cox = T, x_vals = (0:1))
df2 <- generate_survival_times(1000, b1 = 1, b2 = 2, is_cox = T, x_vals = (0:1))


hist(df1$t_true)
hist(df2$t_true)

dat <- generate_survival_data(1000000, (0:2), surv_is_cox = T, b1 = -1, b2 = log(2), cens_is_cox = T)
cox_naive_model(dat)
dat
proportion_censored(dat)
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



#LEGACY CODE


#generate_censoring_times <- function(data, is_cox = T, b1, b2){
#  n <- get_row_length(data)
#  a <- get_a_values(data)
#  x <- get_x_values(data)
#  if(is_cox){
#    censoring_times <- conditional_censoring_cox(n = n, a = a, x = x)
#  }
#  else{
#    censoring_times <- conditional_censoring_not_cox()
#  }
#  return(censoring_times)
#}
#


#  conditional_censoring_cox <- function(n,a,x, b1 = 1, b2 = 1){
#  censoring_times <- rexp(n, rate = 1 + a*b1 + x*b2)
#  return(censoring_times)
#}

