


n <- 500
data <- generate_survival_data(n, ba_t = -1, bx_t = log(2), bz_t = log(2),
                                     ba_c = 1, bx_c = log(2), bz_c = 1)
EIF_test <- EIF(data, T_corr = F, max_time = 2)
EIF_test[4]

n_sims <- 10
simulation <- matrix(data = NA, ncol = 4, nrow = n_sims)
max_time <- 5

for (i in 1:n_sims){
  n <- 10000
  data <- generate_survival_data(n, ba_t = -1, bx_t = log(2), bz_t = log(2),
                                 ba_c = 1, bx_c = log(2), bz_c = 1)
  
  EIF_sim <- EIF(data, T_corr = F, max_time = max_time)
  simulation[i,1] <- mean(EIF_sim[[1]])
  simulation[i,2] <- mean(EIF_sim[[2]])
  simulation[i,3] <- mean(EIF_sim[[3]])
  simulation[i,4] <- EIF_sim[[4]]
  print(paste0('Simulation ',i,' done'))
}
mean(simulation[,4])


sum(between(simulation[,4], 0.7310462715-1.96*sd(simulation[,4]),0.7310462715+1.96*sd(simulation[,4])))/100

weak_eq
strong_eq

######################################################
#       10000 individuals, 100 simulations
#       ba_t = -1
#       ba_c = 1
#       bx_t = log(2)
#       bz_t = log(2)
#       bx_c = log(2)
#       bz_c = 1
#       max_time = 10
#       z unif(0,1) and x unif(-1,1)
#       Both censoring and survival correctly specified
sim_10000 <- simulation
######################################################


#Truth 0.7310585710
0.7076703 #n_individuals = 500
0.7121777 #n_individuals = 1000
0.722216  #n_individuals = 1500
0.7268907 #n_individuals = 2000
0.7240742 #n_individuals = 2500
0.7269909 #n_individuals = 3000
0.7265684 #n_individuals = 3500
0.7267914 #n_individuals = 4000

simulate_prop <- matrix(data = NA, nrow = 500, ncol = 1)

for (i in 1:500){
n <- 1000
data <- generate_survival_data(n, ba_t = 0, bx_t = log(2), bz_t = log(2),
                               ba_c = 1, bx_c = log(2), bz_c = 1)

model <- oracle_model(data)
model_cov <- get_oracle_covar(data)
max_time <- 1000

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
simulate_prop[i] <- mean(res)
print(paste0('Simulation ',i,' done'))
}
mean(simulate_prop)










P10 <- P_treatment_extend_survival(oracle_model(data), max_time = 2, get_oracle_covar(data))

EIF_test <- EIF(data, T_corr = T, Cens_corr = T, max_time = 2)
mean(EIF_test[[4]])




























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






########################################
#Martingale estimation again again
########################################
n <- 1000
test_data <- generate_survival_data(n, ba_t = 1, bx_t = log(2), bz_t = log(2), surv_is_cox = T,
                                    ba_c = log(2), bx_c = 1, bz_c = 1, cens_is_cox = T)


#Fitting model and returning beta estimates
model <- cox.aalen(Surv(T_obs, Uncensored) ~ prop(A) + prop(X), data = test_data)
beta_hat <- model$gamma


#Returning jump times tau and cumulative baseline hazard cum_basehaz
cum <- model$cum
tau <- cum[,1]
cum_basehaz <- cum[,2]
K <- length(tau)


#Calculating cumulative hazard
cum_haz <- exp(cbind(test_data$A, test_data$X) %*% beta_hat) %*% cum_basehaz


#Estimating survival function
surv_func <- exp(-cum_haz)

plot(tau,surv_func[1,], type = 'l', ylim = c(0,1))
for (i in 1:49){
  lines(tau, surv_func[i+1, ])
}


#Estimating change in cumulative hazard
dL <- cbind(0,cum_haz[,-1] - cum_haz[,- (K)])

#Determing at_risk
at_risk <- outer(X = test_data$T_obs, Y = tau, FUN = ">=")


#Determining dN
indicator_jump_time <- outer(test_data$T_obs, tau, FUN = "==")
indicator_uncensored <- test_data$Uncensored

dN <- indicator_jump_time*indicator_uncensored


#Determining dM
dM <- dN - at_risk * dL


surv_func <- exp(-cum_haz)

plot(tau,cumsum(dM[1,]), type = 'l', ylim = c(-3,1.1))
for (i in 2:n){
  lines(tau, cumsum(dM[i, ]))
}

dM_cumsum <- dM %>% apply(FUN = cumsum, MARGIN = 1) %>% t()
colnames(dM_cumsum) <- tau

set.seed(1)
n <- 100
test_data <- generate_survival_data(n, ba_t = 1, bx_t = log(2), bz_t = log(2), surv_is_cox = T,
                                    ba_c = log(2), bx_c = 1, bz_c = 1, cens_is_cox = T)

martingales <- estimate_martingales(train_data = test_data)
at_risk <- martingales$at_risk
dN <- martingales$mod_dN
dL <- martingales$mod_dL
dM <- martingales$mod_dM
tau <- martingales$jump_times



plot(tau,cumsum(dM[1,]), type = 'l', ylim = c(-3,1), ylab = 'dM')
for (i in 1:99){
  lines(tau, cumsum(dM[i, ]))
}

dM_cumsum <- dM %>% apply(FUN = cumsum, MARGIN = 1) %>% t()
mean(dM_cumsum[,ncol(dM_cumsum)])



S_hats <- get_S_hats(test_data)
dim(S_hats$Shat)
View(S_hats$Shat)
View(get_K_hat(test_data))




set.seed(1)
n <- 1000
train_data <- generate_survival_data(n, ba_t = 1, bx_t = log(2), bz_t = log(2), surv_is_cox = T,
                                    ba_c = log(2), bx_c = 1, bz_c = 1, cens_is_cox = T)


working_model <- cox_naive_model(train_data)
mod_cum <- get_cum_hazard(working_model)
mod_jump_times <- get_jump_times(mod_cum)
mod_cum_baseline <- get_cum_baseline_hazard_times(mod_cum)
n_jump_times <- get_row_length(mod_cum)
beta_hat <- get_parameter_estimates(working_model)

train_A <- get_A_values(train_data)
train_X <- get_X_values(train_data)
#Jeg synes, at dette er ret suspekt - vi fitter baseline hazard på noget vi har estimeret det ud fra. Selvfølgelig bliver det meget pænt!!
mod_est_cum_haz <- exp(cbind(train_A, train_X) %*% beta_hat) %*% mod_cum_baseline

#List of observed times
T_obs <- get_observed_times(train_data)

#Generating at risk matrix. Notice the equality!
at_risk <- outer(X = T_obs, Y = mod_jump_times, FUN = ">=")

#Change in cumulative hazard per individual (rows) per stopping time (columns).
mod_dL <- cbind(0, mod_est_cum_haz[,-1] - mod_est_cum_haz[,-ncol(mod_est_cum_haz)])

#Estimating the change in counting process for each individual i dN_i(t) = I(T_i = t, Delta = 1):
indicator_jump_time <- outer(T_obs, mod_jump_times, FUN = "==")
indicator_uncensored <- get_censor_status(train_data)

mod_dN <- indicator_jump_time*indicator_uncensored

cumulative_count_obs <- cumsum(colSums(mod_dN))

cumulative_risk <-cumsum(colSums(at_risk*mod_dL))

mod_dM <- mod_dN - at_risk*mod_dL

counting_process_plot <- 'plot'














############################################################
#                 Plotting functions
############################################################
set.seed(1)
n <- 1000
data <- generate_survival_data(n, ba_t = 0, bx_t = log(2), bz_t = log(2), surv_is_cox = T,
                                     ba_c = 1, bx_c = 0, bz_c = 1, cens_is_cox = T, x_cont = F, x_vals = c(0,1,2,3))

par(mfrow = c(1,4))
df1 <- data[,!names(data) %in% c( "Z", "Uncensored")]
df2 <- pivot_longer(data = df1, cols = c("T_true", "C", "T_obs"), names_to = "Source", )
df2$A <- as.factor(df2$A)
df3 <- df2 %>% filter(Source == "T_true")
df3
ggplot(df3)+ geom_density(aes(value, fill = A), alpha = 0.4) + facet_wrap(~X)


P_T0_T1 <- P_treatment_extend_survival(data)$res
hist(P_T0_T1)
par(mfrow = c(3,2))

props <- propensity(data)


martingales <- estimate_martingales(train_data = data)
at_risk <- martingales$at_risk
dN <- martingales$mod_dN
dL <- martingales$mod_dL
dM <- martingales$mod_dM
tau <- martingales$jump_times


plot(tau,cumsum(dN[1,]), type = 'l', ylab = 'Counting processes')
for (i in 2:100){
  lines(tau, cumsum(dN[i, ]))
}


plot(tau,cumsum((at_risk*dL)[1,]), type = 'l', ylim = c(0,3), ylab = 'Cumulative hazards')
for (i in 2:100){
  lines(tau, cumsum((at_risk*dL)[i, ]))
}


plot(tau,cumsum(dM[1,]), type = 'l', ylim = c(-2,1), ylab = 'Martingales')
for (i in 2:100){
  lines(tau, cumsum(dM[i, ]))
}



S_hats <- get_S_hats(data)
S_hat <- S_hats$Shat
S_hat0 <- S_hats$Shat_0
S_hat1 <- S_hats$Shat_1

plot(tau, S_hat[1,], type = 'l', ylim = c(0,1), ylab = 'Survival function')
for (i in 2:100){
  lines(tau, S_hat[i,])
}

plot(tau, S_hat1[1,], type = 'l', ylim = c(0,1), ylab = 'Survival function')
lines(tau, S_hat0[1,], col = 'red')
for (i in 2:100){
  lines(tau, S_hat1[i,])
  lines(tau, S_hat0[i,], col = 'red')
}


K_C_hat <- get_K_hat(data)
plot(tau, K_C_hat[1,], type = 'l', ylab = 'K_C')
for (i in 2:100){
  lines(tau, K_C_hat[i,])
}












