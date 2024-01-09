##################
#Packages
##################
library(tidyverse)
library(ggplot2)
library(survival)
library(survminer)
#-----------------


########
#OPTIONS
#-------
########
n <- 100000                          #number of simulations
n_alpha <- 100                       #Coarseness of alphas
alph <- seq(1,2.5,length = n_alpha) #Defining the alphas we want to observe
gam <- c(0.2, 1.1, 2)               #Choose the gamma values of interest for the exponential distributions
n_gamma <- length(gam)
set.seed(100)

###############################
#RANDOMNESS
#the only randomness introduced
#------------------------------
###############################
t <- runif(n, 0, 1)

###########################
#DATA STORAGE AND FUNCTIONS
#--------------------------
###########################
#Splitting the data into two halves
t0 <- t[1:(n/2)]                    #data that follows cox-model
t1 <- t[(n/2+1):n]                  #data that does not follow cox-model

###################
#Quantile functions
###################
F_0 <- function(t){                 #quantile function for h0
  -log(1-t)
}
F_1 <- function(t,a){               #generalized quantile function for h1
  (-log(1-t))^(1/a)                 #with variable a representing alpha
}

#############
#Data storage
#############
F_0_sim <- F_0(t0)                  #storage of baseline values of t as a vector


#############################################################
#SIMULATION
#------------------------------------------------------------
#do not change, only run
#Generates stopping times based on the presets provided above
#############################################################
n_total <- length(gam) * length(gam) * length(alph)
values <- vector("numeric", length = n_total * 4)
index <- 1

for (gamma0 in gam){
  for (gamma1 in gam){
    C_0 <- rexp(n/2, gamma0)            #Generating random censor times for h0
    C_1 <- rexp(n/2, gamma1)            #Generating random censor times for h1
    C_full <- c(C_0, C_1)               #Collecting all random censor times in one vector
    
    for (i in alph){
      
      F_1_sim <- F_1(t1, i)             #Transforming uniform t's to F_1 based on alpha
      full_data <- c(F_0_sim,           #Collecting all t's into one vektor 
                     F_1_sim)           #alpha-dependent t's
      
      ##########################
      #SIMULATING CENSORED TIMES
      ##########################
      
      censored_times <- pmin(full_data, C_full) #censoring the original times

      ################################
      #Creating data frame from matrix
      ################################
      censored_df <- data.frame(Time = censored_times,          #list of (possibly) censored times
                                Censored = full_data < C_full,      #Indicator for censoring
                                Distribution = rep(0:1, each = n/2) #indicator for original distribution
                                )
        
      
      #######################
      #Running cox-regression
      #######################
      cox <- coxph(Surv(Time,Censored)~Distribution, data = censored_df)
      
      estimated_beta <- cox$coefficients #returning beta-values
      
      values[index:(index + 3)] <- c(estimated_beta, i, gamma0, gamma1)
      index <- index + 4
      
    }
  }
}


# Create a matrix from 'values' data
analysis_matrix <- matrix(
  data = values, 
  byrow = TRUE,
  nrow = n_total,
  ncol = 4
)

# Create a data frame for ggplot
analysis_matrix_df <- data.frame(analysis_matrix)

# Rename the columns for clarity
colnames(analysis_matrix_df) <- c('Beta', 'Alpha', 'Gamma0', 'Gamma1')

# Convert data types as needed
analysis_matrix_df <- analysis_matrix_df %>%
  mutate(
    Beta = as.numeric(Beta),
    Alpha = as.numeric(Alpha),
    Gamma0 = as.factor(Gamma0),
    Gamma1 = as.factor(Gamma1)
  )

# Create a new column for gamma values as a combination of gamma0 and gamma1
analysis_matrix_df$Gamma_Values <- paste0('(', analysis_matrix_df$Gamma0, ', ', analysis_matrix_df$Gamma1, ')')

# Create a ggplot with customized settings
ggplot(analysis_matrix_df, aes(x = Alpha, y = Beta, color = Gamma0, linetype = Gamma1)) +
  geom_line(linewidth = 0.6) +
  
  # Modify axis labels
  labs(
    x = expression(theta),
    y = expression(hat(beta))
  ) +
  

  # Remove grid lines and customize the axis line
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    text=element_text(size=16,  family="serif"),
    legend.position = c(0.2, 0.2),
    legend.box = "horizontal"
  ) +
  
  # Customize color and linetype scales
  #scale_color_manual(values = c("blue", "red","darkgreen")) +
  scale_color_grey() +
  scale_linetype_manual(values = c("dotted", "dashed", "solid"))
