library(tidyverse)
library(ggplot2)
library(survival)
library(survminer)
library(readr)
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
library(timereg)
library(xtable)
library(stargazer)
library(texreg)


################################
#   Reading and preparing data
################################
df <- read_csv("~/Desktop/Studie/3/Bachelor/Data/Operation_hos_boern/data_til_sara.csv")

attach(df)

#Constructing global time, indicating time until event or censoring
df$time=df$time_until_exam
df$time=df$time+runif(length(df$time),-1,1)/10000 #Adding random noise not to have jumps at the same time
df$time[!is.na(df$time_until_glau)]=df$time_until_glau[!is.na(df$time_until_glau)]
df$status <- df$glau_eye_single == 1


#Creating centered variables
df$age_at_surg_center <- scale(df$age_at_surg, center = TRUE, scale = FALSE)
df$axis_lenght_center <- scale(df$axis_lenght, center = TRUE, scale = FALSE)
df$young <- df$age_at_surg < 6*30


################################################################
#             Cumulative baseline hazard
################################################################

fit_all <- coxph(Surv(time,status) ~ ster_regi + iol_single + age_at_surg + axis_lenght + num_re_op, data=df)
fit_0 <- coxph(Surv(time, status) ~ 1, data = df)
cumbasehaz_coxph <- basehaz(fit_all, center = FALSE)[,1]

base_covs <- df %>% mutate(ster_regi = 0, iol_single = 0, age_at_surg = 0, axis_lenght = 0, num_re_op = 0)
fit_0_surv <- survfit(fit_all, newdata = base_covs)

times_coxph <- fit_0_surv$time
cumbasehaz_coxph <- fit_0_surv$cumhaz[,1]
std.err_coxph <- fit_0_surv$std.err[,1]

ci_lower_cumbasehaz_coxph <- cumbasehaz_coxph - std.err_coxph
ci_upper_cumbasehaz_coxph <- cumbasehaz_coxph + std.err_coxph



sample_list <- list()
for (i in 1:10000){
  sample_indices <- sample(nrow(df), size = 5000, replace = TRUE)
  sample_data <- df[sample_indices, ]
  
  fit_sample <- coxph(Surv(time,status) ~ ster_regi + iol_single + age_at_surg + axis_lenght + num_re_op, data=sample_data)
  
  cumbasehaz_coxph_fit <- basehaz(fit_sample, centered = FALSE)[,1]
  times <- basehaz(fit_sample, centered = FALSE)[,2]
  sample_list[[length(sample_list) + 1]] <- list(times, cumbasehaz_coxph_fit)
}

#plot(basehaz(fit_all)[,2], cumbasehaz_coxph, type = 'l')
#for (i in 1:100){
#  lines(sample_list[[i]][[1]], sample_list[[i]][[2]], lty = 2)
#}

means <- c()
for (i in 1:10000){
  means[i] <- mean(sample_list[[i]][[2]])
}

ci_lower_index <- which(means == max(means[which(means < quantile(means, 0.025))]))
ci_upper_index <- which(means == min(means[which(means > quantile(means, 0.975))]))

plot(basehaz(fit_all)[,2], cumbasehaz_coxph, type = 'l', ylim = c(0,1))

plot(times_coxph, cumbasehaz_coxph, type = 'l', ylim = c(0,1), 
     ylab = 'Cummulative baseline hazard', 
     xlab = 'Days', 
     main = 'Cummulative baseline hazard')
lines(times_coxph, ci_lower_cumbasehaz_coxph, lty = 2)
lines(times_coxph, ci_upper_cumbasehaz_coxph, lty = 2)
lines(sample_list[[ci_lower_index]][[1]], sample_list[[ci_lower_index]][[2]], lty = 4)
lines(sample_list[[ci_upper_index]][[1]], sample_list[[ci_upper_index]][[2]], lty = 5)


################################################################
#             Fitting different models
################################################################






plot((fit_cox_aalen$cum)[,1], (fit_cox_aalen$cum)[,2], col = 'blue')





summary(fit_all)

#In the above model only age_at_surg and num_re_op are significant

fit_selected1 <- coxph(Surv(time,status) ~ factor(ster_regi) + age_at_surg + num_re_op, data=df) 
summary(fit_selected1)

#In the above model iol_single is significant. ster_regi is still not significant

fit_selected2 <- coxph(Surv(time,status) ~ factor(iol_single) + age_at_surg + num_re_op, data=df)
summary(fit_selected2)

fit_ster_reg <- coxph(Surv(time,status) ~ factor(iol_single) + age_at_surg + num_re_op, data=df)




################################################################
#             Check proportional hazard assumption
################################################################

#Fitting models for plotting
full_surv <- survfit(Surv(time, status) ~ factor(ster_regi) + factor(iol_single) + age_at_surg + axis_lenght + num_re_op, data = df)
ster_regi_surv <- survfit(Surv(time, status) ~ factor(ster_regi), data = df)
iol_single_surv <- survfit(Surv(time, status) ~ factor(iol_single), data = df)
young_surv <- survfit(Surv(time, status) ~ factor(young), data = df)


################################################################
#     Kaplan Meier curves for specific covariates
################################################################



plot(ster_regi_surv$time, ster_regi_surv$strata)

#Kaplan Meier curve with ster_regi
ggsurvplot(ster_regi_surv,
                           pval = FALSE, conf.int = F,
                           palette = c("darkblue", "darkred"),
                           risk.table.col = "strata", # Change risk table color by groups
                           ggtheme = theme_bw(), # Change ggplot2 theme
                           xlim = c(0,3000),
                           ylim = c(0.8, 1))

#We see that the survival of the high dose group seems to higher in the first 1500 days but afterwards the opposite seems to be the case
#Note that the curves cross indicating a violation of the constant proportional hazards


#Kaplan Meier curve with iol_single
KM_iol_single <- ggsurvplot(iol_single_surv,
                           pval = FALSE, conf.int = F,
                           risk.table.col = "strata", # Change risk table color by groups
                           ggtheme = theme_bw(), # Change ggplot2 theme
                           xlim = c(0,3000),
                           ylim = c(0.0, 1))

#We notice that introducing an artifical lens seem to increase survival probability.
#The curves are not crossing


#Kaplan Meier curve with six_months_or_young
KM_young <- ggsurvplot(young_surv,
                           pval = FALSE, conf.int = F,
                           risk.table.col = "strata", # Change risk table color by groups
                           ggtheme = theme_bw(), # Change ggplot2 theme
                           xlim = c(0,3000),
                           ylim = c(0.8, 1))
#We see that it seems that older children (older than 6 months) tend to have a higher survival probability
#than younger children (younger than 6 months)
#Curves are not crossing

KM_plots <- list('KM_ster_regi' = KM_ster_regi, 
                 'KM_iol_single' = KM_iol_single, 
                 'KM_young' = KM_young)

arrange_ggsurvplots(KM_plots, ncol = 3, nrow = 1)

################################################################
#                       Cumulative hazard
################################################################

#estimating cumulative baseline hazard for fit_all:
cum_base_haz <- survfit(fit_all, newdata = data.frame(ster_regi = 0, iol_single = 0, age_at_surg = 0, axis_lenght = 0, num_re_op = 0))

fit_all_ca <- cox.aalen(Surv(time,status) ~ prop(ster_regi) + 
                       prop(iol_single) + prop(age_at_surg) + prop(axis_lenght) + prop(num_re_op), data=df)

fit_all_ca$cum[,2] - cum_base_haz
cum_base_haz

ggplot(cum_base_haz, aes(x = time, y = hazard)) + 
  geom_line() + 
  theme_bw() +
  xlab('Time in days') +
  ylab('Cumulative baseline hazard')



#ster_regi
CH_ster_regi <- ggsurvplot(ster_regi_surv,
                                        risk.table.col = "strata", # Change risk table color by groups
                                        ggtheme = theme_bw(), # Change ggplot2 theme
                                        fun = "cumhaz",
                                        xlim = c(0,3000),
                                        ylim = c(0,0.2))
#Once again the crossing of the cumulative hazard rates indicates a violation of constant proportional hazards

#iol_single
CH_iol_single <- ggsurvplot(iol_single_surv,
                               risk.table.col = "strata", # Change risk table color by groups
                               ggtheme = theme_bw(), # Change ggplot2 theme
                               fun = "cumhaz",
                               xlim = c(0,3000),
                               ylim = c(0,0.2))

#young
CH_young <- ggsurvplot(young_surv,
                                        risk.table.col = "strata", # Change risk table color by groups
                                        ggtheme = theme_bw(), # Change ggplot2 theme
                                        fun = "cumhaz",
                                        xlim = c(0,3000),
                                        ylim = c(0,0.2))


CH_plots <- list('CH_cox_ster_regi' = CH_ster_regi, 
                 'CH_cox_iol_single' = CH_iol_single, 
                 'CH_cox_young' = CH_young)
arrange_ggsurvplots(CH_plots, ncol = 2, nrow = 2, risk.table.height = 0.4)




################################################################
#       Cumulative hazard and difference in cumulative hazard ratio
################################################################

#ster_regi
fit_aalen_ster_regi = aalen(Surv(time,status) ~ factor(ster_regi) + factor(iol_single) + factor(young), max.time = 1800, data=df)
summary(fit_aalen_ster_regi)

+ factor(iol_single) +  + axis_lenght + num_re_op

par(mfrow=c(2,2))
plot(ster_regi_surv, fun = "cumhaz", col = c("steelblue", "red"), xlim = c(0,3000), xlab = 'Time', ylab = 'Cumulative hazard')
legend("bottomright", legend = c("ster_regi0", "ster_regi1"), col = c("steelblue", "red"), lty=1)
plot(fit_aalen_ster_regi)


#iol_single
fit_aalen_iol_single = aalen(Surv(time, glau_eye_single==1) ~ factor(iol_single), data=df)
par(mfrow=c(2,2))
plot(iol_single_surv, fun = "cumhaz", col = c("steelblue", "red"), xlim = c(0,3000), xlab = 'Time', ylab = 'Cumulative hazard')
legend("bottomright", legend = c("iol_single0", "iol_single1"), col = c("steelblue", "red"), lty=1)
plot(fit_aalen_iol_single, xlim = c(0,3000))





################################################################
#       Lin et al. and Wei plot
################################################################
#Fits cox model unweighted
cox_model <- cox.aalen(Surv(time, glau_eye_single==1) ~ prop(ster_regi) + prop(age_at_surg_center) + prop(iol_single) + prop(axis_lenght_center) + prop(num_re_op), 
                       data = df)
summary(cox_model)
par(mfrow=c(2,3))
plot(cox_model,score=1)


#Fits cox model weighted
cox_model_w <- cox.aalen(Surv(time, glau_eye_single==1) ~ prop(ster_regi) + prop(age_at_surg_center) + prop(iol_single) + prop(axis_lenght_center) + prop(num_re_op), 
                       weighted.test = 1, data = df)
summary(cox_model_w)
par(mfrow=c(2,3))
plot(cox_model_w,score=1)




################################################################
#                   End of current analysis
################################################################






#Fitting stratified focusing on ster_regi
cox_strat_model <- cox.aalen(Surv(time, status) ~ -1 + ster_regi + prop(age_at_surg) + prop(iol_single) + prop(axis_lenght) + prop(num_re_op),
                       data = df,
                       max.time = 7,
                       n.sim = 100)
summary(cox_model)
par(mfrow=c(1,2))
plot(cox_model,score=1) 



cox_model <- coxph(Surv(time, status) ~ -1 + ster_regi + prop(age_at_surg) + prop(iol_single) + prop(axis_lenght) + prop(num_re_op), data = df)
plot(cox.aalen(cox_model, data = df), var = 'ster_regi')




cox_strat_model_full <- cox.aalen(Surv(time, status) ~ -1 + ster_regi + prop(age_at_surg) + prop(iol_single) + prop(axis_lenght) + prop(num_re_op),
                                  data = df, covariance = 1)





















#Kaplan Meier without covariates

fit.all = survfit(Surv(time, glau_eye_single==1) ~ 1, data=df)
plot(fit.all, xlab="Time since operation", mark.time=F)


#Kaplan Meier with ster_regi
fit.sterregi = survfit(Surv(time, glau_eye_single==1) ~ ster_regi, data=df)

plot(fit.sterregi, xlab="Time since operation", 
     ylim=c(0,1), col=c("red", "green"), mark.time=F)
legend(3500, 0.9, legend=c("ster_regi=0", "ster_regi=1"), text.col=c("red", "green"))


#Kaplan Meier with iol_single
fit.iol = survfit(Surv(time, glau_eye_single==1) ~ iol_single, data=df)

plot(fit.iol, xlab="Time since operation", 
     ylim=c(0,1), col=c("red", "green"), mark.time=F)
legend(3500, 0.9, legend=c("iol=0", "iol=1"), text.col=c("red", "green"))





#LOG RANK TEST
survdiff(Surv(time, glau_eye_single==1) ~ ster_regi, data=df) #ikke signifikant
#H0 (ingen forskel mellem ster_regi grupperne), accepteres (p-værdi på 0,6)

survdiff(Surv(time, glau_eye_single==1) ~ iol_single, data=df) #signifikant

survdiff(Surv(time, glau_eye_single==1) ~ axis_lenght, data=df)




#COX
cox.sterregi = coxph(Surv(time, glau_eye_single==1) ~ factor(ster_regi), data=df)

summary(cox.sterregi) #effekten af ster_regi er ikke signifikant. exp(coef)=0,846, relative risk mellem ster_regi = 0 og 1

table(df$ster_regi)


#med iol_single (kunstig linse)
cox.kunstig_linse = coxph(Surv(time, glau_eye_single==1) ~ factor(iol_single), data=df)

summary(cox.kunstig_linse) #exp(coef)=0,12 - dem med ios=1 har meget lavere risiko for eventet end dem med ios=0. 
#netop 0,12 gange mindre hazard risk.



##FLERE COVARIATER
cox.more = coxph(Surv(time, glau_eye_single==1) ~ factor(ster_regi)+factor(iol_single), data=df)
summary(cox.more)
#iol_single er signifikant, det er ster_regi ikke. exp(coef) for ster_regi er 0.6, den er 0,11 for ios_single


cox.more = coxph(Surv(time, glau_eye_single==1) ~ factor(ster_regi)+factor(iol_single)+axis_lenght, data=df)
summary(cox.more)

cox.more = coxph(Surv(time, glau_eye_single==1) ~ factor(ster_regi)+factor(iol_single)+axis_lenght + age_at_surg, data=df)
summary(cox.more)

#med log af age
cox.log_age = coxph(Surv(time, glau_eye_single==1) ~ log(age_at_surg), data=df)
summary(cox.log_age)




#Validation
fit.cox = cox.aalen(Surv(time, glau_eye_single==1) ~ prop(factor(ster_regi))+prop(factor(iol_single)), data=df)
summary(fit.cox)

par(mfrow=c(2,2))
plot(fit.cox, score=T, xlab="Time")

#ser ikke så godt ud...

##time-varying effects

fit.aalen = aalen(Surv(time, glau_eye_single==1) ~ factor(ster_regi), data=df)
summary(fit.aalen)

par(mfrow=c(2,2))
plot(fit.aalen)



fit.aalen.iol = aalen(Surv(time, glau_eye_single==1) ~ factor(iol_single), data=df, max.time=3000)
plot(fit.aalen.iol)












