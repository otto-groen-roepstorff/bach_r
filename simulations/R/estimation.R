

################################################################
#Generating data for testing estimators
################################################################
n <- 10000
test_data <- generate_survival_data(n, ba_t = -2, bx_t = -2, bz_t = 1, surv_is_cox = T,
                                       ba_c = log(2), bx_c = 1, bz_c = 1, cens_is_cox = T)






################################################################
#Estimating P(T1 > T0 | W)
################################################################

test_model <- cox.aalen(Surv(T_obs, Uncensored) ~ prop(A) + prop(X), data = test_data)
cum <- test_model$cum

################################
#Estimating P(T* > t | A = 1, W)
################################


#Jump times and cumulative baseline hazard:
no_jumps <- dim(cum)[1] - 1    #Number of jumps T > t
tau <- cum[,1][-1]      #Jump times tau1, tau2,...,tau_no_jumps
cumbasehaz <- cum[,2]   #Cumulative baseline hazard   Lambda0(tau0), Lambda0(tau1),...,Lambda0(tau_no_jumps)


#Betahat
betahat <- test_model$gamma


#Estimated cumulative hazard assuming the cox-model
cumhaz_hat <- (exp(cbind(test_data$A, test_data$X) %*% betahat) %*% cumbasehaz)             # cumhaz(t | A, X)
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
T_obs_times <- matrix(data = test_data$T_obs, nrow = n, ncol = no_jumps)
jump_times <- matrix(data = tau, nrow = n, ncol = no_jumps, byrow = T)

#Matrix displaying if individual i is at risk at jump time tau_j
at_risk <- as.numeric(t(t(T_obs_times)) > jump_times)

#Estimating the change in the cumulative hazard for individual i dL_i(t| A_i, W_i) = L_i(t| A_i, W_i) - L_i(t-1 | A_i, W_i):
dL <- cumhaz_hat[,-1] - cumhaz_hat[,-(no_jumps+1)]

#Estimating the change in counting process for individual i dN_i(t) = I(T_i = t, Delta = 1):
dN <- as.numeric(t(t(T_obs_times)) == jump_times) * matrix(data = test_data$Uncensored, nrow = n, ncol = no_jumps)

#Estimating change in the martingale for individual i dM_i(t| A, X) = dN_i(t) - I(T_i > t) L_i(t | A_i, W_i):
dM <- dN - at_risk*dL







################################################################
#               Estimating martingale integral
################################################################

#The martingale integral sum_{T_j} 1/(S(T_j| A , X)*K_C(T_j | A, X)) dM(T_j | A, X)

#First estimating the survival function in jump times S(T_j | A, X):
Shat <- exp(-cumhaz_hat[,-1])


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



#The problem now is that the censoring times and failure times do not jump at the same time and we wish to look at the jump times
#of the survival function. That is we wish to return Khat not in the jump times tau_cens but in tau. Since Khat is constant in between
#jump times, we can simply take the value of Khat corresponding to the largest value of tau_cens that is lower than tau
Khat = matrix(NA, nrow = n, ncol = no_jumps) #Maybe this can also be vectorized??
for(i in 1:no_jumps){
  Khat[,i] = Khat_temp[,max((1:no_jumps_cens)[tau_cens<=tau[i]])]
}



#Finally the martingale integral becomes
rowSums(dM/(Khat*Shat))











