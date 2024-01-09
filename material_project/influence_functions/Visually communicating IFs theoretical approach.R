##################
#Packages
##################
library(tidyverse)  #For tidy data
library(ggplot2)    #For plotting
library(e1071)      #For moment calculations
library(Pareto)
library(gridExtra)
#-----------------


##################
#Generating distributions for plotting simple case ( n = 100 from later code)
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
density_plot <- ggplot(submodels_long, aes(x = X.axis, color = epsilon, fill = variable)) +
  geom_line(aes(y = value), linetype = "solid") +
  scale_color_gradient(low = "lightgrey", high = "black", name = "Epsilon value") +  # Adjust low and high colors as needed
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    text=element_text(size=16,  family="serif"),
    legend.position = c(0.8, 0.8),
    legend.box = "horizontal"
  ) +
  labs(title = "(A) Single path from convex combinations of distributions",
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
functional_plot <- ggplot(df_functionals, aes(x = epsilon, y = functionals)) +
  geom_line(linetype = "solid") + 
  geom_text(aes(x = 0.07, y = 0.5, label = 'true value')) +
  geom_segment(x = 0, xend = 1, y = 0.75, yend = 0.75, linetype = 'dotted') + #Plug-in line
  geom_text(aes(x = 0.07, y = 0.76, label = 'plug-in estimate')) +
  geom_segment(aes(x = 0, xend = 1, y = functionals[101] - slope, yend = functionals[101]), linetype = 'dashed') + #One-step line
  geom_text(aes(x = 0.07, y = 0.44, label = '1-step estimate')) +
  geom_segment(aes(x = 0, xend = 0, y = functionals[1] - 1/20, yend = functionals[1]), linetype = 4) + #R2 term
  geom_text(aes(x = 0.023, y = 0.48, label = 'R2')) +
  geom_point(x = 0, y = 0.5, size = 2.5) + geom_point(x = 0, y = 0.75, shape = 15, size = 2.5) + geom_point(x = 0, y = 0.45, shape = 17, size = 2.5) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    text=element_text(size=16,  family="serif"),
    legend.position = c(0.8, 0.8),
    legend.box = "horizontal"
  ) +
  labs(title = "(B) Functional values along the path",
       x = "Epsilon",
       y = "Functional") + ylim(c(0.4,0.8))

#window plot
grid.arrange(density_plot, functional_plot, ncol = 2)


#Saving plots
wd <- getwd()
save_wd <- paste0(wd,'/Visually communicating IFs')

ggsave(paste0(save_wd,'/density_plot.png'), plot = density_plot, width = 6, height = 4, dpi = 300)
ggsave(paste0(save_wd,'/functional_plot.png'), plot = functional_plot, width = 6, height = 4, dpi = 300)


#Extension to multiple models with n dependency-----------------
#Consistent models

##################
#Generating distributions for plotting
##################
#Set seed for reproducibility
set.seed(1)

#Simulate observations for the empirical distribution for plotting
sim_data <- seq(0,3,0.02)

#Empirical distribution for plotting with different n
n <- c(10, 20, 40, 80, 160, 320)#, 640, 1000, 2000, 4000, 8000, 16000, 32000, 64000, 100000, 200000, 400000, 1000000)

cons_empd_dist <- list() #create list for densities for the different n values
for(i in 1:length(n)){ 
  emp_dist_i <- list(dexp(sim_data, 1 + 5/sqrt(n[i])))
  cons_empd_dist <- c(cons_empd_dist, emp_dist_i )
}

#Theoretical distribution for plotting
theo_dist <- dexp(sim_data, 1)


#Mixture distribution
epsilon <- seq(0,1,1/5) 
epsilon_dist <- list()
for (i in 1:length(n)) {
  epsilon_dist[[i]] <- lapply(epsilon, function(w) {
    w * cons_empd_dist[[i]] + (1 - w) * theo_dist
  })
}

#Collecting models into dataframes
submodels <- data.frame("X.axis" = sim_data) #Column of x-values for densities

for (j in 1:length(n)) {
  col_names <- paste0("Epsilon_", epsilon, ", n = ", n[j])
  submodels[col_names] <- do.call(cbind, epsilon_dist[[j]])
}


#Creating a long dataframe used for plots
submodels_long <- reshape2::melt(data = submodels, id.vars = "X.axis")
submodels_long$epsilon <- rep(epsilon, each = length(sim_data))
submodels_long$n <- rep(n, each = length(sim_data)*length(epsilon))




#Plotting densities stratified on epsilon value
density_plot <- ggplot(submodels_long, aes(x = X.axis, color = epsilon, fill = variable)) +
  geom_line(aes(y = value), linetype = "solid") +
  scale_color_gradient(low = "lightgrey", high = "black", name = "Epsilon value") +  # Adjust low and high colors as needed
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    text=element_text(size=16,  family="serif"),
    legend.position = c(0.8, 0.9),
    legend.box = "horizontal"
  ) +
  labs(title = "A single path formed from convex combinations of distributions",
       x = "Z",
       y = "Density") + xlim(c(0,3))+
  facet_wrap(.~n  )


#The product densities
fp_pure <- function(n){ #product of two f_n
  return((1+5/sqrt(n))/2)
}

fp_mix <- function(n){ #proudct of f_n and f
  (1+5/sqrt(n))/(2+5/sqrt(n))
}

fp_true <- function(n){ #product of two f
  return(1/2)
}

#The theoretical expected density of the epsilon submodel 3/4*epsilon^2 + (1-epsilon)^2*(1+5/sqrt(n))/2 + 2*epsilon(1-epsilon)*(1+5/sqrt(n))/(2+5/sqrt(n))
exp_dens <- function(x,n){
  fp_pure(n)*x^2 + (1-x)^2*fp_true(n) + 2*x*(1-x)*fp_mix(n)
}


#Defining epsilon grid for which to calculate the functional and calculating the functionals
epsilon <- seq(0,1,1/100)
functionals <- data.frame(epsilon = epsilon)
for(i in n){
  clnam <- paste("n =", i)
  res <- sapply(epsilon, FUN = function(x) exp_dens(x = x, n = i))
  functionals[[clnam]] <- res
}

#Collectiong functionals and epsilons into dataframe for plotting
df_functionals2 <- pivot_longer(functionals, !epsilon, names_to = c("sample_size"))
#finding slopes
slopes <- 2*fp_pure(n)-2*fp_mix(n)


error_terms <- 2*fp_mix(n)-fp_pure(n)-fp_true(n)

one_step_end <- as.vector(t(functionals[nrow(functionals),-1]))
one_step_start <- one_step_end -slopes

#slope data
segment_data <- data.frame(
  slope = slopes,
  one_step_estimate = one_step_start,
  plug_in_estimate = one_step_end,
  sample_size = factor(x = c(colnames(functionals[,-1])), levels = c(colnames(functionals[,-1])))
)

segment_data_long <- pivot_longer(segment_data, cols = c("one_step_estimate", "plug_in_estimate"), names_to = c("estimate_type"))

#Plotting 
functional_plot <- ggplot(df_functionals2, aes(x = epsilon, y = value)) +
  geom_abline(intercept = 0.5, slope = 0, linetype = "dotdash")+
  geom_line(aes(color = sample_size), size = 1) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    text=element_text(size=16,  family="serif"),
    legend.position = c(0.25, 0.8),
    legend.box = "horizontal"
  ) +
  labs(title = "Functional values along the path",
       x = "Epsilon",
       y = "Functional") + 
  ylim(c(min(segment_data$one_step_estimate),max(segment_data$plug_in_estimate)))+
  geom_segment(data = segment_data, aes(x = 0, xend = 1, y = one_step_start, yend = one_step_end, color = factor(sample_size)), linetype = "dashed", size = 1)+
  geom_segment(data = segment_data, aes(x = 0, xend = 1, y = one_step_end, yend = one_step_end, color = factor(sample_size)), linetype = "dotted", linewidth  =1)+
  geom_point(data = segment_data_long, aes(color = sample_size, shape = estimate_type, y = value), x = 0)+
  scale_shape_manual(values = c(17, 15))

functional_plot


#abs(error_terms)<abs(one_step_end-0.5)
#bias comparison plot
bias_comparison <- ggplot() +
  geom_line(data = segment_data, aes(x = n, y = abs(one_step_estimate - 0.5)), linetype = "dotted") +
  geom_line(data = segment_data, aes(x = n, y = abs(plug_in_estimate - 0.5)), linetype = "dashed") +
  ylim(c(0, max(segment_data$plug_in_estimate) - 0.5)) +
  labs(
    title = "Absolute Bias Comparison",
    y = "|Bias|"
  ) +
  annotate(geom = "text", y = 0.5, x = 60, label = "plug-in estimator")+
  annotate(geom = "text", y = 0.2, x = 60, label = "one-step estimator")+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    text=element_text(size=16,  family="serif"),
    legend.position = c(0.25, 0.8),
    legend.box = "horizontal"
  )
                       
grid.arrange(functional_plot, bias_comparison, nrow = 2)

convergence_plot
segment_data$plug_in_estimate
segment_data$one_step_estimate

#Saving plots
wd <- getwd()
save_wd <- paste0(wd,'/Visually communicating IFs')



