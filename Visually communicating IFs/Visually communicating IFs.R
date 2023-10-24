##################
#Packages
##################
library(tidyverse)  #For tidy data
library(ggplot2)    #For plotting
library(e1071)      #For moment calculations
#-----------------


############################################################
#        Visually communicating intuition for IFs
############################################################

#We follow the example of the paper using T(P) = E(p(z)^2) for a 1-dimensional r.v.



############################################################
#     Defining function for making mixed distributions
############################################################

Fmixture_model <- function(seed, n, n_epsilon, empdist, theodist, ...){
  
  #Seed: Choose seed for reproducability
  #n: Choose number of simulations
  #n_epsilon: Choose number of epsilon-submodels with epsilon in [0,1]
  #empdist: Choose empirical distribution (rexp, rnorm, runif,...)
  #theodist: Choose theoretical distribution (dexp, dnorm, dunif,...)
  #...: Choose parameters for your chosen distribution e.g. 0,1 for standard normal
  
  
  #Set seed
  set.seed(seed)
  
  #Empirical distribution
  sim_data <- empdist(n, ...)
  density_est <- density(sim_data, n = 1024)
  emp_dist <- density_est$y
  
  #Theoretical distribution
  theo_dist <- theodist(density_est$x, ...)
  
  
  #Mixture distribution
  epsilon <- seq(0,1,1/n_epsilon)
  
  epsilon_dist <- sapply(epsilon, function(w) {
    w * emp_dist +
      (1 - w) * theo_dist
  }, simplify = FALSE)
  
  #Collecting models into dataframes
  submodels <- data.frame("X.axis" = density_est$x)
  
  for (i in (1:length(epsilon))){
    submodels[[paste0("Epsilon_", epsilon[i])]] <- epsilon_dist[[i]]
  }
  
  #Creating a long dataframe used for plots
  submodels_long <- reshape2::melt(data = submodels, id.vars = "X.axis")
  submodels_long$epsilon <- rep(epsilon, each = length(density_est$x))
  
  all_models <- data.frame(submodels, "Emp.dist" = emp_dist, "Theoretical.dist" = theo_dist)
  
  return(list("all_models" = all_models, "submodels_long" = submodels_long))
}



############################################################
#                 Simulating mixed-model data
############################################################

#Function for simulating data based on the mixed densities
Fmixed_data_sim <- function(seed, Fmixture_model_all_models, n_sims){
  
  
  set.seed(seed) #Set seed for reproduceability
  
  simulated_data <- sapply(Fmixture_model_all_models[,-1], 
         FUN = sample, x = Fmixture_model_all_models[[1]], size = n_sims, replace = TRUE) #Simulate data
  
  return(simulated_data)
}


############################################################
#                 Functional calculations
############################################################

#Function for calculating functionals
Fmixfunctionals <- function(Fmixture_model_all_models, functional, ...){
  Fmixture_model_all_models <- Fmixture_model_all_models[,2:(ncol(Fmixture_model_all_models)-2)] #Select only epsilon-distributions
  est_functionals <- sapply(Fmixture_model_all_models, functional, ...) #Applying functional to each distribution
  
  
  functionals_long <- reshape2::melt(data = est_functionals, id.vars = "Index") #Reshaping to get functionals as vector
  
  est_functionals <- data.frame("epsilon" = seq(0,1,1/(length(functionals_long$value)-1)), 
                                'functionals' = functionals_long$value) #Collecting functionals and corresponding epsilon value
  
  return(est_functionals)
}







############################################################
#             Using functions to make plots and calculations
############################################################
epsilon_model <- Fmixture_model(1, 1000, 5, rexp, dexp, 1)
long_submodels <- epsilon_model$submodels_long

#Plotting
ggplot(long_submodels, aes(x = X.axis, color = epsilon, fill = variable)) +
  # Add multiple geom_line layers for each y-axis column
  geom_line(aes(y = value), linetype = "solid") +
  scale_color_gradient(low = "red", high = "steelblue", name = "Epsilon value") +  # Adjust low and high colors as needed
  theme_classic() +
  labs(title = "A single path formed from convex combinations of distributions",
       x = "Z",
       y = "Density")+ xlim(c(0, 2.5)) + ylim(c(0.035, 1))



dists <- epsilon_model_functional_plot$all_models[-c(1,103,104)]
funcs <- est_functionals$functionals

gateau_deriv_1 <- -mean((dists$Epsilon_1 - funcs[101]))
y_intercept <- funcs[101] - gateau_deriv_1

ifz <- 2*(g(z) - T(G))

slope <- (est_functionals[nrow(est_functionals),2] - est_functionals[nrow(est_functionals) - 1,2])/
  (est_functionals[nrow(est_functionals),1] - est_functionals[nrow(est_functionals) - 1,1])

directional_derivative <- est_functionals$epsilon * slope

#Simulating same data but with more epsilons for plotting functional curve
epsilon_model_functional_plot <- Fmixture_model(1, 10000, 100, rexp, dexp, 0.1)

est_functionals <- Fmixfunctionals(epsilon_model_functional_plot$all_models, functional = moment, 2)

ggplot(est_functionals, aes(x = epsilon, y = functionals)) +
  geom_line(linetype = "solid") + theme_classic() +
  labs(title = "Functional values along the path",
       x = "Epsilon",
       y = "Functional") + ylim(c(0.00044, 0.00052)) + geom_abline(intercept = funcs[101] - slope, slope = slope)






