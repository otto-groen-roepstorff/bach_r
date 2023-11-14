

#Simulate data
full_survival_data <- function(n, dist){
  #Randomness
  #Splitting the data into two halves
  
  
  treatment_times <- dist()
  
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
  
  
#  function(n,a,x, b1 = 1, b2 = 1){
#  censoring_times <- rexp(n, rate = 1 + a*b1 + x*b2)
#  return(censoring_times)
#}

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


generate_survival_times <- function(n,b1, b2, x_vals, is_cox = T){
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
                                         x)
    }
    else{
      a0_times <- quantile_surv_generate_not_cox(q_a0, 
                                             b1 = b1, 
                                             b2 = b2, 
                                             a = 0, 
                                             x)
      a1_times <- quantile_surv_generate_not_cox(q_a1, 
                                             b1 = b1, 
                                             b2 = b2, 
                                             a = 1, 
                                             x)
    }
  
    temp_df0 <- data.frame(t_true = a0_times, a = rep(0, length(a0_times)), x = x)
    temp_df1 <- data.frame(t_true = a1_times, a = rep(1, length(a1_times)), x = x)
    data_frame_storage <- rbind(data_frame_storage, temp_df1, temp_df0)
  }
  
    
  
  return(data_frame_storage)
}


generate_censoring_times <- function(data, is_cox = T, b1, b2){
  n <- get_row_length(data)
  a <- get_a_values(data)
  x <- get_x_values(data)
  if(is_cox){
    censoring_times <- conditional_censoring_cox(n = n, a = a, x = x)
  }
  else{
    censoring_times <- conditional_censoring_not_cox()
  }
  return(censoring_times)
}

generate_censoring_times(data)

generate_survival_data <- function(n, x_vals, b1 = 1, b2 = 1, surv_is_cox = T, cens_is_cox = T){
  survival_data <- generate_survival_times(n,b1 = b1, b2 = b2, x_vals = x_vals, is_cox = surv_is_cox)
  censoring_times <- generate_censoring_times(data = survival_data, is_cox = cens_is_cox, b1 = b1, b2 = b2)
    has_been_censored <- get_t_values(survival_data) > censoring_times
  t_observed <- pmin(get_t_values(survival_data), censoring_times)
  full_data_set <- cbind(survival_data, censoring_times, t_observed, has_been_censored)
  return(full_data_set)
}

data <- generate_survival_data(100,x_vals = c(0,1,2,6, 7, 6), is_cox = T)


