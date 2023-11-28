
library(pracma)

theoretical_value_01 <- function(beta_hat, upper_time){
  beta_a <- beta_hat[1]
  beta_w <- beta_hat[2]
  
  nom_1 <- -expint_E1(upper_time * (exp(beta_a) + 1))
  nom_2 <- expint_E1(upper_time * (exp(beta_a) + 1) * exp(beta_w)) + beta_w
  nom <- nom_1 + nom_2
  
  denom <- (exp(beta_a) + 1) * beta_w
 
  res <- nom/denom
  
  return(res)
}
beta_hat <- c(-1,-log(2))
theoretical_value_01(beta_hat, 10)

theo_a <- function(x){
  return(theoretical_value(c(x, log(1.0001)), 1))
}


x <- data.frame(seq(-log(6), log(6), length = 100))

y <- apply(x, FUN = theo_a, MARGIN = 1)
plot(x$seq..log.6...log.6...length...100., y)
abline(h = 1-1/7)





theoretical_value_11 <- function(beta_hat, upper_time){
  beta_a <- beta_hat[1]
  beta_w <- beta_hat[2]
  
  nom_1 <- -expint_E1(upper_time * (exp(beta_a) + 1) * exp(-beta_w)) + beta_w
  nom_2 <- expint_E1(upper_time * (exp(beta_a) + 1) * exp(beta_w)) + beta_w
  nom <- nom_1 + nom_2
  
  denom <- 2 * beta_w * (exp(beta_a) + 1)
  
  res <- nom/denom
  
  return(res)
}
beta_hat <- c(0,log(2))
theoretical_value_11(beta_hat, 5)








beta <- log(6)
A <- 1
tau_func <- function(tau){
  expint_E1(tau * exp(beta) * (exp(A) + 1)) - expint_E1(tau * exp(-beta) * (exp(A) + 1))  
}

x <- data.frame(seq(0.0000001, 5, length = 100))
y <- apply(x, FUN = tau_func, MARGIN = 1)
plot(x[,1],y)
sum(y < -0.5)

y
