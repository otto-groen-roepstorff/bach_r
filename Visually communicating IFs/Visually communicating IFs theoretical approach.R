##################
#Packages
##################
library(tidyverse)  #For tidy data
library(ggplot2)    #For plotting
library(e1071)      #For moment calculations
library(Pareto)
#-----------------


##################
#Generating distributions for plotting
##################
#Set seed for reproducability
set.seed(1)

#Simulate observations for the empirical distribution for plotting
sim_data <- seq(0,3,0.02)

#Empirical distribution for plotting
emp_dist <- dexp(sim_data, 1.5)

#Theoretical distribution for plotting
theo_dist <- dexp(sim_data, 1)


#Mixture distribution
epsilon <- seq(0,1,1/5)

epsilon_dist <- sapply(epsilon, function(w) {
  w * emp_dist +
    (1 - w) * theo_dist
}, simplify = FALSE)

#Collecting models into dataframes
submodels <- data.frame("X.axis" = sim_data) #Column of observations
for (i in (1:length(epsilon))){
  submodels[[paste0("Epsilon_", epsilon[i])]] <- epsilon_dist[[i]] #columns of probabilities
}

#Creating a long dataframe used for plots
submodels_long <- reshape2::melt(data = submodels, id.vars = "X.axis")
submodels_long$epsilon <- rep(epsilon, each = length(sim_data))


#Plotting densities stratified on epsilon value
ggplot(submodels_long, aes(x = X.axis, color = epsilon, fill = variable)) +
  geom_line(aes(y = value), linetype = "solid") +
  scale_color_gradient(low = "red", high = "steelblue", name = "Epsilon value") +  # Adjust low and high colors as needed
  theme_bw() +
  labs(title = "A single path formed from convex combinations of distributions",
       x = "Z",
       y = "Density") + xlim(c(0,3))





#The theoretical expected density of the epsilon submodel 3/4*epsilon^2 + (1-epsilon)^2 + epsilon(1-epsilon)*6/5
exp_dens <- function(x){
  3/4*x^2 + (1-x)^2 * 1/2 + x*(1-x)*6/5
}


#Defining epsilon values for which to calculate the functional and calculating the functionals
epsilon <- seq(0,1,1/100)
functionals <- sapply(epsilon, FUN = exp_dens)

#Collectiong functionals and epsilons into dataframe for plotting
df_functionals <- data.frame('epsilon' = epsilon, 'functional' = functionals)


#Slope of functional at epsilon = 1 is calculated to 3/10
slope <- 3/10


#Plotting 
ggplot(df_functionals, aes(x = epsilon, y = functionals)) +
  geom_line(linetype = "solid") + 
  geom_segment(aes(x = 0, xend = 1, y = functionals[101] - slope, yend = functionals[101]), linetype = 'dashed') +
  geom_text(aes(x = 0.1, y = 0.44, label = '1-step estimator')) +
  geom_segment(aes(x = 0, xend = 0, y = functionals[1] - 1/20, yend = functionals[1]), linetype = 'dotted') +
  geom_text(aes(x = 0.023, y = 0.48, label = 'R2')) +
  theme_classic() +
  labs(title = "Functional values along the path",
       x = "Epsilon",
       y = "Functional") + ylim(c(0.4,0.8))




