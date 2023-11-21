

################################################################
#Generating data for testing estimators
################################################################
n <- 100
train_data <- generate_survival_data(n, ba_t = 1, bx_t = -1, bz_t = log(6), surv_is_cox = T,
                                     ba_c = log(6), bx_c = 1, bz_c = -1, cens_is_cox = F)
test_data <- generate_survival_data(n, ba_t = 1, bx_t = -1, bz_t = log(6), surv_is_cox = T,
                                    ba_c = log(6), bx_c = 1, bz_c = -1, cens_is_cox = F)



################################################################
#Fitting models
################################################################

train_model <- cox.aalen(Surv(T_obs, Uncensored) ~ prop(A) + prop(X), data = train_data)
test_model <- cox.aalen(Surv(T_obs, Uncensored) ~ prop(A) + prop(X), data = test_data)


################################################################
#Jump times and cumulative baseline hazard:
################################################################

#Train data
cum_train <- train_model$cum
no_jumps_train <- dim(cum_train)[1] - 1           # Number of jumps T > t
tau_train <- cum_train[,1][-1]                    # Jump times tau1, tau2,...,tau_no_jumps


#Test data
cum_test <- test_model$cum
no_jumps_test <- dim(cum_test)[1] - 1
tau_test <- cum_test[,1][-1]



#Betahat from training model
betahat <- train_model$gamma

#test_model_coxph <- coxph(Surv(T_obs, Uncensored) ~ A + X, data = test_data)
#breslow(test_model_coxph)

#Estimated cumulative hazard assuming the cox-model based on training data
cumbasehaz_train <- cum_train[,2]                                                                 # Cumulative baseline hazard Lambda0(tau0), Lambda0(tau1),...,Lambda0(tau_no_jumps)
cumhaz_hat_train <- exp(cbind(train_data$A, train_data$X) %*% betahat) %*% cumbasehaz_train       # cumhaz(t | A, X)


################################################################
#            Estimating change in martingale on train data
################################################################

#Determining at-risk indicators
T_obs_times_train <- matrix(data = train_data$T_obs, nrow = n, ncol = no_jumps_train)
jump_times_train <- matrix(data = tau_train, nrow = n, ncol = no_jumps_train, byrow = T)

#Matrix [i,j] displaying if individual i is at risk at jump time tau_j
at_risk_train <- matrix(data = as.numeric(T_obs_times_train >= jump_times_train), nrow = n, ncol = no_jumps_train)

#Estimating the change in the cumulative hazard for each individual i dL_i(t| A_i, W_i) = L_i(t| A_i, W_i) - L_i(t-1 | A_i, W_i):
dL <- cumhaz_hat_train[,-1] - cumhaz_hat_train[,-(no_jumps_train+1)]

#Estimating the change in counting process for each individual i dN_i(t) = I(T_i = t, Delta = 1):
dN_train <- matrix(data = as.numeric(T_obs_times_train == jump_times_train), nrow = n, ncol = no_jumps_train) * matrix(data = train_data$Uncensored, nrow = n, ncol = no_jumps_train)


plot(tau_train, cumsum(colSums(dN_train)))
lines(tau_train, cumsum(colSums(at_risk_train*dL)), col = 'red')

dM <- dN_train - at_risk_train*dL
plot(tau_train, colSums(dM), type = 'l')

################################################################
#           Estimating change in the martingale on test data
################################################################

#Determining at-risk indicators
T_obs_times_test <- matrix(data = test_data$T_obs, nrow = n, ncol = no_jumps_test)
jump_times_test <- matrix(data = tau_test, nrow = n, ncol = no_jumps_test, byrow = T)

#Matrix [i,j] displaying if individual i is at risk at jump time tau_j
at_risk_test <- matrix(data = as.numeric(T_obs_times_test > jump_times_test), nrow = n, ncol = no_jumps_test)

#Estimating the change in the cumulative hazard for each individual i dL_i(t| A_i, W_i) = L_i(t| A_i, W_i) - L_i(t-1 | A_i, W_i):
#dL <- cumhaz_hat[,-1] - cumhaz_hat[,-(no_jumps+1)]

#Estimating the change in counting process for each individual i dN_i(t) = I(T_i = t, Delta = 1):
dN_test <- matrix(data = as.numeric(T_obs_times_test == jump_times_test), nrow = n, ncol = no_jumps_test) * matrix(data = test_data$Uncensored, nrow = n, ncol = no_jumps_test)

dM <- dN_test - at_risk_train*dL
plot(tau_train, colSums(dM), type = 'l')

plot(tau_test, cumsum(colSums(dN_test)))
lines(tau_train, cumsum(colSums(at_risk_train*dL)), col = 'red')





