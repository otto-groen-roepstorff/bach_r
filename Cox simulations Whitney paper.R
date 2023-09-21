library(tidyverse)
library(ggplot2)
library(survival)
library(survminer)



#Genererer data
t <- runif(2000, 0, 1)

t1 <- t[1:1000]
t2 <- t[1001:2000]




a_seq <- seq(1,2.5,0.5)







#Genererer censoreringer
C_0 <- rexp(1000, 2)

C_1 <- rexp(1000, 2)




cdf <- data.frame(c_data)
colnames(cdf) <- c('Time', 'Censored', 'Distribution')

cox <- coxph(Surv(Time,Censored)~Distribution, data = cdf)


c_data <- cbind(apply(cbind(full_data[,4], C_full), MARGIN = 1, FUN = min), full_data[,4] < C_full, c(rep(0,1000),rep(1,1000)))


cdf <- data.frame(c_data)
colnames(cdf) <- c('Time', 'Censored', 'Distribution')

cox <- coxph(Surv(Time,Censored)~Distribution, data = cdf)
cox$coefficients


#Samlet analyse
n_alpha <- 4
n_gamma <- 3
alph <- seq(1,2.5, lengt=n_alpha)
gam <- c(0.2, 1.1, 2)

n <- 100000
t <- runif(n, 0, 1)


analyse_matrix <- matrix(data = NA, nrow = n_alpha*n_gamma^2, 4)

#Vi arbejder med samme t hver gang
t1 <- t[1:(n/2)]
t2 <- t[(n/2+1):n]

#Kvantilfunktioner
F_0 <- function(t){
  -log(1-t)
}

F_1 <- function(t,a){
  (-log(1-t))^(1/a)
}

F_0_sim <- F_0(t1)
F_0_sim_full <- cbind(F_0_sim, F_0_sim, F_0_sim, F_0_sim)

j = 0

#
for (gamma0 in gam){
  for (gamma1 in gam){
    C_0 <- rexp(n/2, gamma0) #Censoreringer
    C_1 <- rexp(n/2, gamma1)
    C_full <- c(C_0, C_1)
    
    #Genererer data
    F_1_sim <- matrix(data = NA, nrow = length(t2), ncol = length(alph))
    
    for (i in alph){
      F_1_sim <- F_1(t2, i)  #Generer F_1 baseret på alpha og gamma0 og gamma1
      
      full_data <- c(F_0_sim, F_1_sim)
      
      censored_data <- cbind(
        apply(cbind(full_data, C_full), MARGIN = 1, FUN = min), #Censoreret data
                      full_data < C_full, 
                      c(rep(0,n/2),rep(1,n/2))
        )
      censored_df <- data.frame(censored_data)
      colnames(censored_df) <- c('Time', 'Censored', 'Distribution')
      
      cox <- coxph(Surv(Time,Censored)~Distribution, data = censored_df)
      estimated_beta <- cox$coefficients
      
      values <- c(estimated_beta, i, gamma0, gamma1)
      
      
      j = j + 1
      analyse_matrix[j, ] <- values
    }
  }
}


analyse_matrix_df <- data.frame(analyse_matrix)
colnames(analyse_matrix_df) <- c('theta', 'alpha', 'gamma0', 'gamma1')

analyse_matrix_df$gamma_values <-  paste0('(',analyse_matrix_df$gamma0,',',analyse_matrix_df$gamma1,')')

ggplot(analyse_matrix_df, aes(x = alpha, y = theta, color = gamma_values)) + geom_line()








