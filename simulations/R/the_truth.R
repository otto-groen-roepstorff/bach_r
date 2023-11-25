
library(pracma)

theoretical_value <- function(beta_hat, upper_time){
  beta_a <- beta_hat[1]
  beta_w <- beta_hat[2]
  
  nom_1 <- -expint_Ei(-exp(beta_w) * (1 + exp(beta_a)) * upper_time)
  nom_2 <- expint_Ei(-(1 + exp(beta_a))*upper_time) + beta_w
  nom <- nom_1 + nom_2
  
  denom <- (exp(beta_a) + 1) * beta_w
 
  res <- nom/denom
  
  return(as.numeric(res))
}

beta_hat <- c(-1,log(2))
theoretical_value(beta_hat, 4)
