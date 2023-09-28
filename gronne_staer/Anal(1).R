library(tidyverse)
library(survival)
data<- read.csv("øjendata.csv", header= TRUE)
head(data)


data$time=data$time_until_exam
data$time[!is.na(data$time_until_glau)]=data$time_until_glau[!is.na(data$time_until_glau)]

head(data)
data$opt_tek_trt=interaction(data$op_teknik,data$ster_regi)
table(data$opt_tek_trt) # 0 in combination opt.tek=1 (århus) and trt=0
set.seed(99)
data$time=data$time+runif(length(data$time),-1,1)/10000

# now analyses:

table(data$ster_regi,data$op_place, dnn = c("ster_regi","op_place"))
?table

cbind(data$glau_eye_single,data$time,data$time_until_glau,data$time_until_exam)[1:40,]

fit=coxph(Surv(time,glau_eye_single==1)~factor(ster_regi)+factor(iol_single)+age_at_surg+axis_lenght,data=data) 
summary(fit)
fit=coxph(Surv(time,glau_eye_single==1)~factor(opt_tek_trt)+factor(iol_single)+age_at_surg+axis_lenght,data=data) 
summary(fit)
fit=coxph(Surv(time,glau_eye_single==1)~factor(opt_tek_trt)+factor(iol_single),data=data) 
summary(fit)


plot(survfit(Surv(time,glau_eye_single==1)~1,data=data))
plot(survfit(Surv(time,glau_eye_single==1)~ster_regi,data=data))
plot(survfit(Surv(time,glau_eye_single==1)~opt_tek_trt,data=data))



head(data)
summary(data)
data$X


a <- letters[1:3]
table(a, sample(a))                    # dnn is c("a", "")
table(a, sample(a), deparse.level = 0) # dnn is c("", "")
table(a, sample(a), deparse.level = 2) # dnn is c("a", "sample(a)")
