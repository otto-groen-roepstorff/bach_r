library(tidyverse)
library(ggplot2)
library(survival)
library(survminer)
library(readr)
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)


############################
#   Reading and preparing data
############################
df <- read_csv("~/Desktop/Studie/3/Bachelor/Data/Operation_hos_boern/data_til_sara.csv")

#Constructing global time, indicating time until event or censoring
df$time=df$time_until_exam
df$time[!is.na(df$time_until_glau)]=df$time_until_glau[!is.na(df$time_until_glau)]

#Reformatting date variables
df <- df %>% 
  mutate(birth = ymd(birth),
         cat_surg_date = ymd(cat_surg_date),
         glau_date = ymd(glau_date),
         va_date_after = ymd(va_date_after),
         axis_date = ymd(axis_date)
)

#Creating survival objects and curves
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

#Crossing Kaplan-Meier curves is a clear indication that the constant hazard ratio assumption does not hold


