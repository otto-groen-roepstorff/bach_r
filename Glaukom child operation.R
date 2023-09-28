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
df$time[!is.na(df$time_until_glau)]=df$time_until_glau[!is.na(df$time_until_glau)]
df$status <- df$time != df$time_until_exam


#Changing data types for certain covariates
df$op_place <- factor(df$op_place)
df$iol_single <- factor(df$iol_single)
df$glau_eye_single <- factor(df$glau_eye_single)
df$ster_regi <- factor(df$ster_regi)
df$op_teknik <- factor(df$op_teknik)



################################
#   Check proportional hazard assumption based on stratified Cox model
################################

#We focus on the covariates highlighted by the researchers.

#Fitting stratified focusing on ster_regi
cox_strat_model <- cox.aalen(Surv(time, status) ~ -1 + ster_regi + prop(age_at_surg) + prop(iol_single) + prop(axis_lenght) + prop(num_re_op),
                             data = df, covariance = 1)



#Fits cox model
cox_model <- cox.aalen(Surv(time, glau_eye_single==1) ~ prop(ster_regi) + prop(age_at_surg) + prop(iol_single) + prop(axis_length) + prop(num_re_op), 
                       data = df)
summary(cox_model)
par(mfrow=c(2,3))
plot(cox_model,score=1) 

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






#Creating Kaplan-Meier survival curves. Kaplan-Meier estimate the survival function S(t)
fit.all <- survfit(Surv(time, glau_eye_single == 1) ~ 1, data = df)
plot(fit.all)


survfit2(Surv(time, glau_eye_single == 1) ~ 1, data = df) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  )



View(head(df))

df$opt_tek_trt=interaction(df$op_teknik,df$ster_regi)
table(df$opt_tek_trt) # 0 in combination opt.tek=1 (århus) and trt=0
set.seed(99)

df$time=df$time+runif(length(df$time),-1,1)/10000 #Hvorfor tilføjer man noget random noise til tiderne?





#Preproduced analysis
fit=coxph(Surv(time,glau_eye_single==1)~factor(ster_regi)+factor(iol_single)+age_at_surg+axis_lenght,data=df) 
summary(fit)
fit=coxph(Surv(time,glau_eye_single==1)~factor(opt_tek_trt)+factor(iol_single)+age_at_surg+axis_lenght,data=df) 
summary(fit)
fit=coxph(Surv(time,glau_eye_single==1)~factor(opt_tek_trt)+factor(iol_single),data=df) 
summary(fit)


plot(survfit(Surv(time,glau_eye_single==1)~1,data=df))
plot(survfit(Surv(time,glau_eye_single==1)~ster_regi,data=df))
plot(survfit(Surv(time,glau_eye_single==1)~opt_tek_trt,data=df))
plot(survfit(Surv(time,glau_eye_single==1) ~ ster_regi, data=df))

#Crossing Kaplan-Meier curves is a clear indication that the constant hazard ratio assumption does not hold

















par(mfrow=c(1,1))


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

fit.aalen = aalen(Surv(time, glau_eye_single==1) ~ factor(ster_regi), data=df, max.time=3000)
summary(fit.aalen)

par(mfrow=c(2,2))
plot(fit.aalen)



fit.aalen.iol = aalen(Surv(time, glau_eye_single==1) ~ factor(iol_single), data=df, max.time=3000)
plot(fit.aalen.iol)












