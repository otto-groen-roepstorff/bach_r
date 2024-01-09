library(ggpattern)
library(tidyverse)

generate_vis_data <- function(n, ba_t = 0, bz_c = 0, bx_c = 0,
                              x_range = c(-1,1)){
  #variables
  a <- rbinom(n = n, size = 1, 1/2)  
  x <- runif(n, min = min(x_range), max = max(x_range))
  z <- runif(n, 0,1)
  #generate survival and censoring times
  surv_t <-  rgamma(n = n, shape = 7.5, scale = 1)/((exp(ba_t*a)))
  cens_t <- rgamma(n, shape = 6, scale = 1)/exp(bz_c*z + bx_c*x)
  
  #Observed values
  t_obs <- pmin(surv_t, cens_t)
  
  #Status
  delta <- surv_t < cens_t
  
  return(data.frame("True_survival_times" = surv_t, "Censoring_times" = cens_t, "Observed_survival_times" = t_obs, Uncensored = delta, A = a, X = x, Z = z))
}


data <- generate_vis_data(1000, ba_t = -log(2)) #Generate data
df1 <- pivot_longer(data = data, cols = c("True_survival_times", "Censoring_times", "Observed_survival_times"), names_to = "Source", ) #cahnge data format
df1$A <- as.factor(df1$A) #Make A into a factor 

ggplot(df1, aes(x = value))+
  geom_density_pattern(aes(pattern = A),
                       pattern_spacing = 0.03,
                       fill = "transparent", 
                       color = "black",
                       alpha = 0.4, 
                       pattern_density = 0.4) +
  facet_grid(.~Source) + theme_minimal()+
  theme(
    axis.text.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  )+labs(title = "Effect of censoring on observed survival times")

