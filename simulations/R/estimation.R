

#Estimating P(T1 > T0 | W)
n <- 10000
test_data <- torben_generate_survival_times(n, ba = 0, bx = 1)
test_model <- cox.aalen(Surv(T_obs, Status) ~ prop(A) + prop(X), data = test_data)

cum <- test_model$cum


#Estimating P(T* > t | A = 1, W)


#Cumulative jump times:
K <- dim(cum)[1] - 1    #Number of jumps T > t
tau <- cum[,1][-1]      #Jump times tau1, tau2,...,tauK
cumbasehaz <- cum[,2]   #Cumulative baseline hazard   Lambda0(tau0), Lambda0(tau1),...,Lambda0(tauK)


#Betahat
betahat <- test_model$gamma


#Estimated cumulative hazard assuming the cox-model
cumhaz_hat <- (exp(cbind(test_data$A, test_data$X) %*% betahat) %*% cumbasehaz)[,2:(K+1)]   # cumhaz(t | A, X)
cumhaz1_hat <- (exp(cbind(1, test_data$X) %*% betahat) %*% cumbasehaz)[,2:(K+1)]            # cumhaz(t | A = 1, X)
cumhaz0_hat <- (exp(cbind(0, test_data$X) %*% betahat) %*% cumbasehaz)[,2:(K+1)]            # cumhaz(t | A = 0, X)

#Estimating survival function S(t | A = 1, X)
Shat_1 <- exp(-cumhaz1_hat)
Shat_0 <- exp(-cumhaz0_hat)

delta_cumbasehaz <- matrix(data = cumbasehaz[-1], nrow = n, ncol = K, byrow = T) -
                    matrix(data = cumbasehaz[-(K+1)], nrow = n, ncol = K, byrow = T)

integrand <- Shat_1 * Shat_0 * delta_cumbasehaz

mean(exp(betahat[2] * test_data$X) * rowSums(integrand))


