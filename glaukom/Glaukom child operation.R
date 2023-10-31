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


################################
#   Reading and preparing data
################################
df <- read_csv("~/Desktop/Studie/3/Bachelor/Data/Operation_hos_boern/data_til_sara.csv")

#Constructing global time, indicating time until event or censoring
df$time=df$time_until_exam
df$time=df$time+runif(length(df$time),-1,1)/10000 #Adding random noise not to have jumps at the same time
df$time[!is.na(df$time_until_glau)]=df$time_until_glau[!is.na(df$time_until_glau)]
df$status <- df$glau_eye_single == 1


#Creating centered variables
df$age_at_surg_center <- scale(df$age_at_surg, center = TRUE, scale = FALSE)
df$axis_lenght_center <- scale(df$axis_lenght, center = TRUE, scale = FALSE)
df$six_month_or_young <- df$age_at_surg < 6*30



################################################################
#             Fitting different models
################################################################

fit_all <- coxph(Surv(time,status) ~ factor(ster_regi) + factor(iol_single) + age_at_surg + axis_lenght + num_re_op, data=df)
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
ster_regi_surv <- survfit(Surv(time, status) ~ factor(ster_regi), data = df)
iol_single_surv <- survfit(Surv(time, status) ~ factor(iol_single), data = df)
six_month_or_young_surv <- survfit(Surv(time, status) ~ factor(six_month_or_young), data = df)


################################################################
#     Kaplan Meier curves for specific covariates
################################################################

#Kaplan Meier curve with ster_regi
KM_ster_regi <- ggsurvplot(ster_regi_surv,
                           pval = TRUE, conf.int = F,
                           risk.table.col = "strata", # Change risk table color by groups
                           ggtheme = theme_bw(), # Change ggplot2 theme
                           xlim = c(0,3000),
                           ylim = c(0.8, 1))

#We see that the survival of the high dose group seems to higher in the first 1500 days but afterwards the opposite seems to be the case
#Note that the curves cross indicating a violation of the constant proportional hazards


#Kaplan Meier curve with iol_single
KM_iol_single <- ggsurvplot(iol_single_surv,
                           pval = TRUE, conf.int = F,
                           risk.table.col = "strata", # Change risk table color by groups
                           ggtheme = theme_bw(), # Change ggplot2 theme
                           xlim = c(0,5000),
                           ylim = c(0.0, 1))

#We notice that introducing an artifical lens seem to increase survival probability.
#The curves are not crossing


#Kaplan Meier curve with six_months_or_young
KM_six_month_or_young <- ggsurvplot(six_month_or_young_surv,
                           pval = TRUE, conf.int = F,
                           risk.table.col = "strata", # Change risk table color by groups
                           ggtheme = theme_bw(), # Change ggplot2 theme
                           xlim = c(0,3000),
                           ylim = c(0.8, 1))
#We see that it seems that older children (older than 6 months) tend to have a higher survival probability
#than younger children (younger than 6 months)
#Curves are not crossing

KM_plots <- list('KM_ster_regi' = KM_ster_regi, 
                 'KM_iol_single' = KM_iol_single, 
                 'KM_six_month_or_young' = KM_six_month_or_young)

arrange_ggsurvplots(KM_plots, ncol = 2, nrow = 2, risk.table.height = 0.4)

################################################################
#                       Cumulative hazard
################################################################
#ster_regi
CH_ster_regi <- ggsurvplot(cox_ster_regi,
                                        risk.table.col = "strata", # Change risk table color by groups
                                        ggtheme = theme_bw(), # Change ggplot2 theme
                                        fun = "cumhaz",
                                        xlim = c(0,3000),
                                        ylim = c(0,0.2))
#Once again the crossing of the cumulative hazard rates indicates a violation of constant proportional hazards

#iol_single
CH_iol_single <- ggsurvplot(cox_iol_single,
                               risk.table.col = "strata", # Change risk table color by groups
                               ggtheme = theme_bw(), # Change ggplot2 theme
                               fun = "cumhaz",
                               xlim = c(0,3000),
                               ylim = c(0,0.2))

#six_month_or_young
CH_six_month_or_young <- ggsurvplot(cox_six_month_or_young,
                                        risk.table.col = "strata", # Change risk table color by groups
                                        ggtheme = theme_bw(), # Change ggplot2 theme
                                        fun = "cumhaz",
                                        xlim = c(0,3000),
                                        ylim = c(0,0.2))


CH_plots <- list('CH_cox_ster_regi' = CH_ster_regi, 
                 'CH_cox_iol_single' = CH_iol_single, 
                 'CH_cox_six_month_or_young' = CH_six_month_or_young)
arrange_ggsurvplots(CH_plots, ncol = 2, nrow = 2, risk.table.height = 0.4)




################################################################
#       Cumulative hazard and difference in cumulative hazard ratio
################################################################

#ster_regi
fit_aalen_ster_regi = aalen(Surv(time, glau_eye_single==1) ~ factor(ster_regi), data=df)
par(mfrow=c(2,2))
plot(cox_ster_regi, fun = "cumhaz", col = c("steelblue", "red"), xlim = c(0,3000), xlab = 'Time', ylab = 'Cumulative hazard')
legend("bottomright", legend = c("ster_regi0", "ster_regi1"), col = c("steelblue", "red"), lty=1)
plot(fit_aalen_ster_regi, xlim = c(0,3000))


#iol_single
fit_aalen_iol_single = aalen(Surv(time, glau_eye_single==1) ~ factor(iol_single), data=df)
par(mfrow=c(2,2))
plot(cox_iol_single, fun = "cumhaz", col = c("steelblue", "red"), xlim = c(0,3000), xlab = 'Time', ylab = 'Cumulative hazard')
legend("bottomright", legend = c("iol_single0", "iol_single1"), col = c("steelblue", "red"), lty=1)
plot(fit_aalen_iol_single, xlim = c(0,3000))





################################################################
#       Lin et al. and Wei plot
################################################################
#Fits cox model
cox_model <- cox.aalen(Surv(time, glau_eye_single==1) ~ prop(ster_regi) + prop(age_at_surg_center) + prop(iol_single) + prop(axis_lenght_center) + prop(num_re_op), 
                       data = df)
summary(cox_model)
par(mfrow=c(2,3))
plot(cox_model,score=1)


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












