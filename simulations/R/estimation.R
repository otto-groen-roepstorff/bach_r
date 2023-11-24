

################################################################
#Generating data for testing estimators
################################################################
n <- 1000
test_data <- generate_survival_data(n, ba_t = 1, bx_t = -1, bz_t = log(2), surv_is_cox = F,
                                       ba_c = log(2), bx_c = 1, bz_c = 1, cens_is_cox = T)

test_model <- cox.aalen(Surv(T_obs, Uncensored) ~ prop(A) + prop(X), data = test_data)
cum <- test_model$cum
no_jumps <- dim(cum)[1] - 1     # Number of jumps T > t
tau <- cum[,1][-1]              # Jump times tau1, tau2,...,tau_no_jumps
cumbasehaz <- cum[,2]           # Cumulative baseline hazard Lambda0(tau0), Lambda0(tau1),...,Lambda0(tau_no_jumps)


test_model_cens <- cox.aalen(Surv(T_obs, Uncensored == F) ~ prop(A) + prop(X), data = test_data)
cum_cens <- test_model_cens$cum

#Jump times and cumulative baseline hazard:
no_jumps_cens <- dim(cum_cens)[1] - 1     #Number of jumps T > t
tau_cens <- cum_cens[,1][-1]              #Jump times tau_cens1, tau_cens2,...,tau_cens_no_jumps
cumbasehaz_cens <- cum_cens[,2]           #Cumulative baseline hazard   Lambda0(tau0), Lambda0(tau1),...,Lambda0(tau_no_jumps)

#Betahat
betahat_cens <- test_model_cens$gamma

#Estimated cumulative hazard assuming the cox-model
cumhaz_cens_hat <- (exp(cbind(test_data$A, test_data$X) %*% betahat_cens) %*% cumbasehaz_cens)        # cumhaz_cens(t | A, X)

#Estimating survival function S(t | A = 1, X) and S(t | A = 0, X)
Khat_temp <- exp(-cumhaz_cens_hat[,-1])

#The problem now is that the censoring times and failure times do not jump at the same time and we wish to look at the jump times
#of the survival function. That is we wish to return Khat not in the jump times tau_cens but in tau. Since Khat is constant in between
#jump times, we can simply take the value of Khat corresponding to the largest value of tau_cens that is lower than tau
Khat = matrix(NA, nrow = n, ncol = no_jumps) #Maybe this can also be vectorized??
for(i in 1:no_jumps){
  Khat[,i] = Khat_temp[,max((1:no_jumps_cens)[tau_cens<=tau[i]])]
}

View(Khat)


get_K_hat <- function(data, surv_is_cox = T){
  
  n <- get_row_length(data)
  
  ######################################################################
  #                         Getting surv jump times
  ######################################################################
  
  surv_model <- cox.aalen(Surv(T_obs, Uncensored) ~ prop(A) + prop(X), data = data)
  cum <- get_cum_hazard(surv_model)
  no_jumps <- get_row_length(cum)-1                   # Number of jumps T > t
  tau <- get_jump_times(cum)                          # Jump times tau1, tau2,...,tau_no_jumps
  
  ######################################################################
  
  
  ######################################################################
  #                         Getting censoring distribution
  ######################################################################
  #Fitting model for censoring distribution
  if(surv_is_cox){
    mod_cens <- cox.aalen(Surv(get_observed_times(data), get_censor_status(data) == F) ~ 
                            prop(get_A_values(data)) + prop(get_X_values(data)), data = data)
  } else {
    mod_cens <- cox.aalen(Surv(get_observed_times(data), get_censor_status(data) == F) 
                          ~ prop(get_A_values(data)) + prop(get_X_values(data)^2), data = data)
  }
  
  #Jump times and cumulative baseline hazard:
  cum_cens <- get_cum_hazard(mod_cens)
  no_jumps_cens <- get_row_length(cum_cens) - 1               #Number of jumps T > t
  tau_cens <- get_jump_times(cum_cens)                        #Jump times tau_cens1, tau_cens2,...,tau_cens_no_jumps
  cumbasehaz_cens <- get_cum_baseline_hazard_times(cum_cens)  #Cumulative baseline hazard   Lambda0(tau0), Lambda0(tau1),...,Lambda0(tau_no_jumps)
  
  #Betahat
  betahat_cens <- get_parameter_estimates(mod_cens)
  
  #Estimated cumulative hazard assuming the cox-model
  #cumhaz_cens_hat <- (exp(cbind(data$A, data$X) %*% betahat_cens) %*% cumbasehaz_cens)        # cumhaz_cens(t | A, X)
  cumhaz_cens_hat <- predict_cox.aalen(A = get_A_values(data), 
                                       W = get_X_values(data), 
                                       betaHat = betahat_cens, 
                                       cum_base_haz = cumbasehaz_cens, 
                                       n_jumps = no_jumps_cens)
  
  #Estimating survival function S(t | A = 1, X) and S(t | A = 0, X)
  Khat_temp <- exp(-cumhaz_cens_hat)
  
  #The problem now is that the censoring times and failure times do not jump at the same time and we wish to look at the jump times
  #of the survival function. That is we wish to return Khat not in the jump times tau_cens but in tau. Since Khat is constant in between
  #jump times, we can simply take the value of Khat corresponding to the largest value of tau_cens that is lower than tau
  Khat = matrix(NA, nrow = n, ncol = no_jumps) #Maybe this can also be vectorized??
  for(i in 1:no_jumps){
    Khat[,i] = Khat_temp[,max((1:no_jumps_cens)[tau_cens<=tau[i]])]
  }
  
  Khat[is.na(Khat)] <- 1
  
  return(Khat)
}




  
################################################################
#Generating data for testing estimators
################################################################
n <- 1000


test_data <- generate_survival_data(n, ba_t = -1, bx_t = -1, bz_t = log(6), surv_is_cox = T,
                                    ba_c = log(2), bx_c = 1, bz_c = -1, cens_is_cox = T)


P_treatment_extend_survival(test_data)
S_hats <- get_S_hats(test_data)
K_hat <- get_K_hat(test_data)





mod_cens <- cox.aalen(Surv(get_observed_times(test_data), get_censor_status(test_data) == F) 
                        ~ prop(get_A_values(test_data)) + prop(get_X_values(test_data)), data = test_data)

n <- get_row_length(test_data)
cum_haz_matrix_cens <- mod_cens$cum
n_jumps_cens <- get_row_length(cum_haz_matrix_cens)-1
jump_times_cens <- get_jump_times(cum_haz_matrix_cens)
cum_haz <- get_cum_baseline_hazard_times(cum_haz_matrix_cens)
beta_hat <- get_parameter_estimates(mod_cens)
A <- get_A_values(test_data)
X <- get_X_values(test_data)
#Z <- get_Z_values(data)

cumhaz_hat_cens <- predict_cox.aalen(A = A, W = X, betaHat = beta_hat, cum_base_haz = cum_haz, n_jumps = n_jumps_cens)
Khat_temp <- exp(-cumhaz_hat_cens)

######################################################
#Getting real number of jumps and real jump times
model <- cox.aalen(Surv(get_observed_times(test_data), get_censor_status(test_data)) 
                   ~ prop(get_A_values(test_data)) + prop(get_X_values(test_data)), data = test_data)

cum_haz_matrix <- model$cum
n_jumps <- get_row_length(cum_haz_matrix)-1
jump_times <- get_jump_times(cum_haz_matrix)
######################################################


Khat = matrix(NA, nrow = n, ncol = n_jumps) #Maybe this can also be vectorized??
for(i in 1:n_jumps){
  Khat[,i] = Khat_temp[,max((1:n_jumps_cens)[jump_times_cens<=jump_times[i]])]
}










################################################################
#Estimating P(T1 > T0 | W)
################################################################

test_model <- cox.aalen(Surv(T_obs, Uncensored) ~ prop(A) + prop(X), data = test_data)
cum <- test_model$cum

################################
#Estimating P(T* > t | A = 1, W)
################################


#Jump times and cumulative baseline hazard:
no_jumps <- dim(cum)[1] - 1     # Number of jumps T > t
tau <- cum[,1][-1]              # Jump times tau1, tau2,...,tau_no_jumps
cumbasehaz <- cum[,2]           # Cumulative baseline hazard Lambda0(tau0), Lambda0(tau1),...,Lambda0(tau_no_jumps)

#Betahat
betahat <- test_model$gamma


#Estimated cumulative hazard assuming the cox-model
cumhaz_hat <- (exp(cbind(test_data$A, test_data$X) %*% betahat) %*% cumbasehaz)[,-1]          # cumhaz(t | A, X)
cumhaz1_hat <- (exp(cbind(1, test_data$X) %*% betahat) %*% cumbasehaz)[,-1]            # cumhaz(t | A = 1, X)
cumhaz0_hat <- (exp(cbind(0, test_data$X) %*% betahat) %*% cumbasehaz)[,-1]            # cumhaz(t | A = 0, X)

#Estimating survival function S(t | A = 1, X) and S(t | A = 0, X)
Shat_1 <- exp(-cumhaz1_hat)
Shat_0 <- exp(-cumhaz0_hat)


#Estimating the change in the cumulative baseline hazard (dA_0)
delta_cumbasehaz <- matrix(data = cumbasehaz[-1], nrow = n, ncol = no_jumps, byrow = T) -
  matrix(data = cumbasehaz[-(no_jumps+1)], nrow = n, ncol = no_jumps, byrow = T)

#Defining the integrand
integrand <- Shat_1 * Shat_0 * delta_cumbasehaz

#Estimating P(T1 > T0 | W)
exp(betahat[2] * test_data$X) * rowSums(integrand)

#Taking the mean over each individual
mean(exp(betahat[2] * test_data$X) * rowSums(integrand))






################################################################
#               Estimating change in the martingale
################################################################

#Determining at-risk indicators
T_obs_times_test <- matrix(data = test_data$T_obs, nrow = n, ncol = no_jumps_test)
jump_times_test <- matrix(data = tau_test, nrow = n, ncol = no_jumps_test, byrow = T)

T_obs_times_train <- matrix(data = train_data$T_obs, nrow = n, ncol = no_jumps)
jump_times_train <- matrix(data = tau, nrow = n, ncol = no_jumps, byrow = T)

#Matrix [i,j] displaying if individual i is at risk at jump time tau_j
at_risk <- matrix(data = as.numeric(T_obs_times_train > jump_times_train), nrow = n, ncol = no_jumps)

#Estimating the change in the cumulative hazard for each individual i dL_i(t| A_i, W_i) = L_i(t| A_i, W_i) - L_i(t-1 | A_i, W_i):
dL <- cumhaz_hat[,-1] - cumhaz_hat[,-(no_jumps+1)]

#Estimating the change in counting process for each individual i dN_i(t) = I(T_i = t, Delta = 1):
dN <- matrix(data = as.numeric(T_obs_times_train == jump_times_train), nrow = n, ncol = no_jumps) * matrix(data = train_data$Uncensored, nrow = n, ncol = no_jumps)

#Estimating change in the martingale for individual i dM_i(t| A, X) = dN_i(t) - I(T_i > t) dL_i(t | A_i, W_i):
dM <- dN - at_risk*dL

plot(tau, cumsum(colSums(dN)))
points(tau, cumsum(colSums(at_risk*dL)), col = 'red')


plot(tau,colSums(dM))


################################################################
#               Estimating martingale integral
################################################################

#The martingale integral sum_{T_j} 1/(S(T_j| A , X)*K_C(T_j | A, X)) dM(T_j | A, X).

#We have two integrals one which depend on u and one which depend on t. How to handle this? We are first integrating from 0 to t and
#then from 0 to tau. I assume that the first integral should be of the same dimension as S(t | A, X), but not sure.

#First estimating the survival function in jump times S(T_j | A, X):
Shat <- exp(-cumhaz_hat)


############################
#     Censoring distribution
############################
#In addition we need the censoring distribution. We proceed exactly like with the survival times, but this time with Uncensored == F
test_model_cens <- cox.aalen(Surv(T_obs, Uncensored == F) ~ prop(A) + prop(X), data = test_data)
cum_cens <- test_model$cum

#Jump times and cumulative baseline hazard:
no_jumps_cens <- dim(cum_cens)[1] - 1     #Number of jumps T > t
tau_cens <- cum_cens[,1][-1]              #Jump times tau_cens1, tau_cens2,...,tau_cens_no_jumps
cumbasehaz_cens <- cum_cens[,2]           #Cumulative baseline hazard   Lambda0(tau0), Lambda0(tau1),...,Lambda0(tau_no_jumps)

#Betahat
betahat_cens <- test_model_cens$gamma

#Estimated cumulative hazard assuming the cox-model
cumhaz_cens_hat <- (exp(cbind(test_data$A, test_data$X) %*% betahat_cens) %*% cumbasehaz_cens)        # cumhaz_cens(t | A, X)

#Estimating survival function S(t | A = 1, X) and S(t | A = 0, X)
Khat_temp <- exp(-cumhaz_cens_hat)


test_model <- cox.aalen(Surv(T_obs, Uncensored) ~ prop(A) + prop(X), data = test_data)
cum <- test_model$cum

#Jump times and cumulative baseline hazard:
no_jumps <- dim(cum)[1] - 1     # Number of jumps T > t
tau <- cum[,1][-1]              # Jump times tau1, tau2,...,tau_no_jumps
cumbasehaz <- cum[,2]           # Cumulative baseline hazard Lambda0(tau0), Lambda0(tau1),...,Lambda0(tau_no_jumps)

#The problem now is that the censoring times and failure times do not jump at the same time and we wish to look at the jump times
#of the survival function. That is we wish to return Khat not in the jump times tau_cens but in tau. Since Khat is constant in between
#jump times, we can simply take the value of Khat corresponding to the largest value of tau_cens that is lower than tau
Khat = matrix(NA, nrow = n, ncol = no_jumps) #Maybe this can also be vectorized??
for(i in 1:no_jumps){
  Khat[,i] = Khat_temp[,max((1:no_jumps_cens)[tau_cens<=tau[i]])]
}


#Finally the martingale integral becomes
rowSums(dM/(Khat*Shat))











