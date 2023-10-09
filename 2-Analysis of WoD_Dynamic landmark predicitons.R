####################################
####################################
###  Author: Valery Fuh Ngwa    ####
### Dynamic Landmark Predictions####
####################################
####################################

require(dynpred)
library(penalized)

###set working directory to LandmarkModels

setwd("C:/Users/Owner/Documents/Survival_Analysis/DynPred_JointModels/DynPreds_JM/PhD_Study_1/Article_1/Analysis_of_EDSS")

####################################################################
####################################################################
###  Regression coefficients (standard errors) in the two-stage ####
### Cox regression approach                                    #####
####################################################################
####################################################################


##Load in data from stage 1.

load("xzmat.RDATA") 
tt <- sort(unique(xzmat$Obstime[xzmat$status==1]))
xzmat2 <- survSplit(data=xzmat, cut=tt, end="tstop", start="Tstart", event="status")
xzmat2$lnt <- log(xzmat2$time+1)

##Please note.
#zbeta if the clinical-environmental PI (CEPI)
#xbeta is the genetic prognostic index (GPI)
#xzbeta is the clinical-evn-genotypic PI (CEGPI)

xbeta<-xzmat$xbeta; zbeta<-xzmat$zbeta;  xzbeta<-xzmat$xzbeta

##Fit a time-fixed cox model on the CEPI

clinfixed <- coxph(Surv(tstop, status) ~ zbeta+cluster(auslongid), data = xzmat, method="breslow", robust=TRUE)
clinfixed
# Same model in longer data
clinfixed <- coxph(Surv(Tstart, tstop, status) ~ zbeta+cluster(auslongid), data = xzmat2, method="breslow", robust=TRUE)
clinfixed
AIC(clinfixed)
data.frame(loglik0=clinfixed$loglik[1],loglik1=clinfixed$loglik[1],chisq=2*diff(clinfixed$loglik))

#Fit a time varying cox model on the CEPI
clintime <- coxph(Surv(Tstart, tstop, status) ~ zbeta + zbeta:lnt, data = xzmat2, method="breslow")
clintime
AIC(clintime)
data.frame(loglik0=clintime$loglik[1],loglik1=clintime$loglik[1],chisq=2*diff(clintime$loglik))
anova(clinfixed,clintime)

## Fit a time-fixed Cox on GPI

genfixed <- coxph(Surv(Tstart, tstop, status) ~ xbeta+cluster(auslongid), data = xzmat2, method="breslow", robust=TRUE)
genfixed
AIC(genfixed)
data.frame(loglik0=genfixed$loglik[1],loglik1=genfixed$loglik[1],chisq=2*diff(genfixed$loglik))

#Time varying Cox on GPI
gentime <- coxph(Surv(Tstart, tstop, status) ~ xbeta + xbeta:lnt, data = xzmat2, method="breslow")
gentime
AIC(gentime)
data.frame(loglik0=gentime$loglik[1],loglik1=gentime$loglik[1],chisq=2*diff(gentime$loglik))
anova(genfixed,gentime)

###Fit the combined model CEPI + GPI
xzfixed <- coxph(Surv(Tstart, tstop, status) ~ zbeta + xbeta, data = xzmat2, method="breslow")
xzfixed
AIC(xzfixed)
data.frame(loglik0=xzfixed$loglik[1],loglik1=xzfixed$loglik[1],chisq=2*diff(xzfixed$loglik))

#Fit Time-fixed Cox model of super learner (this is the combination of CEPI and GPI into one variable i.e xzbeta)

xzmat2$xzbeta <-coef(xzfixed)[1]*xzmat2$zbeta+coef(xzfixed)[2]*xzmat2$xbeta
xzmat$xzbeta <-coef(xzfixed)[1]*xzmat$zbeta+coef(xzfixed)[2]*xzmat$xbeta

combfixed <- coxph(Surv(Tstart, tstop, status) ~ xzbeta, data = xzmat2, method="breslow")
combfixed
AIC(combfixed)
data.frame(loglik0=combfixed$loglik[1],loglik1=combfixed$loglik[1],chisq=2*diff(combfixed$loglik))

#Fit Time-fixed Cox model of super learner
combtime <- coxph(Surv(Tstart, tstop, status) ~ xzbeta + xzbeta:lnt, data = xzmat2, method="breslow")
combtime
AIC(combtime)
data.frame(loglik0=combtime$loglik[1],loglik1=combtime$loglik[1],chisq=2*diff(combtime$loglik))
anova(combfixed,combtime)



###############################################################################
###############################################################################
### Univariate and multivariate Cox regression on clinical and
### genemic cross-validated prognostic indices in different landmark data sets
###############################################################################
###############################################################################
# Landmark data sets
LMs <- seq(0,5,by=1)
w <- 5
for (i in 1:length(LMs)) {
        LM <- LMs[i]
        cat("\n\nLandmark time point:",LM,"\n\n")
        LMdata <- cutLM(data=xzmat,outcome=list(time="tstop",status="status"),
                        LM=LM,horizon=LM+w,covs=list(fixed=c("zbeta","xbeta","xzbeta"),timedep=NULL))
        cox.clin <- coxph(Surv(LM,tstop,status)~zbeta, data=LMdata, method="breslow")
        print(data.frame(B=round(cox.clin$coef,2),SE=round(summary(cox.clin)$coef[3],2), chisq=round(2*diff(cox.clin$loglik),2)))
        cox.gen <- coxph(Surv(LM,tstop,status)~xbeta, data=LMdata, method="breslow")
        print(data.frame(B=round(cox.gen$coef,2),SE=round(summary(cox.gen)$coef[3],2), chisq=round(2*diff(cox.gen$loglik),2)))
        cox.comb <- coxph(Surv(LM,tstop,status)~xzbeta, data=LMdata, method="breslow")
        print(data.frame(B=round(cox.comb$coef,2),SE=round(summary(cox.comb)$coef[3],2), chisq=round(2*diff(cox.comb$loglik),2)))
}

#################

###############################################################################
###############################################################################
### Estimated regression coefficients (standard errors) for the
### proportional baselines (ipl*) landmark super models without and with
### (landmark-dependent) linear landmark interactions
###############################################################################
###############################################################################

# Stacked landmark data set based on a finer grid
LMs <- seq(0,5,by=0.1)
w <- 2.5
LMdata <- cutLM(data=xzmat,outcome=list(time="tstop",status="status"),
                LM=0,horizon=w,covs=list(fixed=c("zbeta","xbeta","xzbeta", "auslongid"),timedep=NULL))
for (i in 2:length(LMs)) {
        LM <- LMs[i]
        LMdata <- rbind(LMdata,cutLM(data=xzmat,outcome=list(time="tstop",status="status"),
                                     LM=LM,horizon=LM+w,covs=list(fixed=c("zbeta","xbeta","xzbeta", "auslongid"),timedep=NULL)))
}

f1 <- function(t) 1
f2 <- function(t) (t/5)
f3 <- function(t) (t/5)^2
# Explicitly code interactions of treatment with LM
LMdata$clin1 <- LMdata$zbeta*f1(LMdata$LM)
LMdata$clin2 <- LMdata$zbeta*f2(LMdata$LM)
LMdata$clin3 <- LMdata$zbeta*f3(LMdata$LM)
LMdata$gen1 <- LMdata$xbeta*f1(LMdata$LM)
LMdata$gen2 <- LMdata$xbeta*f2(LMdata$LM)
LMdata$gen3 <- LMdata$xbeta*f3(LMdata$LM)
LMdata$comb1 <- LMdata$xzbeta*f1(LMdata$LM)
LMdata$comb2 <- LMdata$xzbeta*f2(LMdata$LM)
LMdata$comb3 <- LMdata$xzbeta*f3(LMdata$LM)
# ipl model
##clinical + genetic
iplstrat <- coxph(Surv(LM,tstop,status) ~ clin1 + clin2 + clin3+gen1+gen2+gen3+ + strata(LM) +
                          cluster(auslongid), data=LMdata, method="breslow")
iplstrat
###clinical only
iplstratclin1 <- coxph(Surv(LM,tstop,status) ~ clin1 +  strata(LM) + 
                               cluster(auslongid), data=LMdata, method="breslow")
iplstratclin1

iplstratclin2 <- coxph(Surv(LM,tstop,status) ~ clin1 + clin2 + clin3+ strata(LM) + 
                               cluster(auslongid), data=LMdata, method="breslow")
iplstratclin2

iplstratgen1 <- coxph(Surv(LM,tstop,status) ~ gen1 +  strata(LM) + 
                              cluster(auslongid), data=LMdata, method="breslow")
iplstratgen1

iplstratgen2 <- coxph(Surv(LM,tstop,status) ~ gen1 + gen2 + gen3+ strata(LM) + 
                              cluster(auslongid), data=LMdata, method="breslow")
iplstratgen2

iplstratcomb1 <- coxph(Surv(LM,tstop,status) ~ comb1 + strata(LM) + 
                               cluster(auslongid), data=LMdata, method="breslow")
iplstratcomb1

iplstratcomb2 <- coxph(Surv(LM,tstop,status) ~ comb1 + comb2 + comb3+ strata(LM) + 
                               cluster(auslongid), data=LMdata, method="breslow")
iplstratcomb2

# ipl* model
g1 <- function(t) f2(t)
g2 <- function(t) f3(t)
LMdata$LM1 <- g1(LMdata$LM)
LMdata$LM2 <- g2(LMdata$LM)

iplcov <- coxph(Surv(LM,tstop,status) ~ clin1 + clin2 + +clin3+gen1 + gen2 + +gen3+LM1 + LM2 + cluster(auslongid), data=LMdata, method="breslow")
iplcov

# ipl* models, separate and combined, with and without interactions
iplclin1 <- coxph(Surv(LM,tstop,status) ~ clin1 + LM1 + LM2 + cluster(auslongid), data=LMdata, method="breslow")
iplclin1
iplclin2 <- coxph(Surv(LM,tstop,status) ~ clin1 + clin2 +clin3+ LM1+LM2+ cluster(auslongid), data=LMdata, method="breslow")
iplclin2
iplgen1 <- coxph(Surv(LM,tstop,status) ~ gen1 + LM1 + LM2 + cluster(auslongid), data=LMdata, method="breslow")
iplgen1
iplgen2 <- coxph(Surv(LM,tstop,status) ~ gen1 + gen2 +gen3+ LM1 + LM2 + cluster(auslongid), data=LMdata, method="breslow")
iplgen2
iplcomb1 <- coxph(Surv(LM,tstop,status) ~ comb1 + LM1 + LM2 + cluster(auslongid), data=LMdata, method="breslow")
iplcomb1
iplcomb2 <- coxph(Surv(LM,tstop,status) ~ comb1 + comb2 +comb3+ LM1 + LM2 + cluster(auslongid), data=LMdata, method="breslow")
iplcomb2
iplcov1 <- coxph(Surv(LM,tstop,status) ~ clin1 + gen1 + LM1 + LM2 + cluster(auslongid), data=LMdata, method="breslow")
iplcov1
iplcov2 <- coxph(Surv(LM,tstop,status) ~ clin1 + clin2 + gen1 + gen2 + LM1 + LM2 + cluster(auslongid), data=LMdata, method="breslow")
iplcov2

###############################################################################
###############################################################################
### Dynamic fixed width predictions for without and with landmark
### interactions
###############################################################################
###############################################################################

bhclin1 <- basehaz(iplclin1,centered=FALSE)
bhclin2 <- basehaz(iplclin2,centered=FALSE)
bhgen1 <- basehaz(iplgen1,centered=FALSE)
bhgen2 <- basehaz(iplgen2,centered=FALSE)
bhcomb1 <- basehaz(iplcomb1,centered=FALSE)
bhcomb2 <- basehaz(iplcomb2,centered=FALSE)
bhcov1 <- basehaz(iplcov1,centered=FALSE)
bhcov2 <- basehaz(iplcov2,centered=FALSE)

### Prediction plot, uses 25 and 75% quantiles
clin.qs <- quantile(zbeta,probs=c(0.25,0.75))
clin.qs
gen.qs <- quantile(xbeta,probs=c(0.25,0.75))
gen.qs
# The corresponding quantiles of the super learner in the data
# for these four individuals
cl<-coef(xzfixed)[1];gl<-coef(xzfixed)[2]
outer(cl*clin.qs, gl*gen.qs, "+")
Fn <- ecdf(xzmat$xzbeta)
Fn(outer(cl*clin.qs, gl*gen.qs, "+"))

## Dynamic predictions for clinical landmark-fixed
sfclin1 <- survfit(iplclin1,newdata=data.frame(clin1=clin.qs,LM1=0,LM2=0),censor=FALSE)
sfclin11 <- data.frame(time=sfclin1[1]$time,surv=sfclin1[1]$surv)
sfclin11$Haz <- -log(sfclin11$surv)
sfclin11$haz <- diff(c(0,sfclin11$Haz))
sfclin12 <- data.frame(time=sfclin1[2]$time,surv=sfclin1[2]$surv)
sfclin12$Haz <- -log(sfclin12$surv)
sfclin12$haz <- diff(c(0,sfclin12$Haz))

tpred <- c(sfclin11$time,sfclin11$time-w)
tpred <- c(0,sort(tpred[tpred>0 & tpred<=10]))
npred <- length(tpred)
Fwclin11 <- Fwclin12 <- rep(NA,npred)
for (i in 1:npred) {
        tp <- tpred[i]
        sfi1 <- sfclin11
        sfi1$haz <- sfi1$haz*exp(iplclin1$coef["LM1"]*g1(tp)+iplclin1$coef["LM2"]*g2(tp))
        sfi1$Haz <- cumsum(sfi1$haz)
        tmp <- approx(sfclin11$time, sfi1$Haz, c(tp,tp+w), method="constant",
                      yleft=0, yright=max(sfi1$Haz))$y
        Fwclin11[i] <- 1-exp(-diff(tmp))
        sfi2 <- sfclin12
        sfi2$haz <- sfi2$haz*exp(iplclin1$coef["LM1"]*g1(tp)+iplclin1$coef["LM2"]*g2(tp))
        sfi2$Haz <- cumsum(sfi2$haz)
        tmp <- approx(sfclin12$time, sfi2$Haz, c(tp,tp+w), method="constant",
                      yleft=0, yright=max(sfi2$Haz))$y
        Fwclin12[i] <- 1-exp(-diff(tmp))
}

## Dynamic predictions for genomic landmark-fixed
sfgen1 <- survfit(iplgen1,newdata=data.frame(gen1=gen.qs,LM1=0,LM2=0),censor=FALSE)
sfgen11 <- data.frame(time=sfgen1[1]$time,surv=sfgen1[1]$surv)
sfgen11$Haz <- -log(sfgen11$surv)
sfgen11$haz <- diff(c(0,sfgen11$Haz))
sfgen12 <- data.frame(time=sfgen1[2]$time,surv=sfgen1[2]$surv)
sfgen12$Haz <- -log(sfgen12$surv)
sfgen12$haz <- diff(c(0,sfgen12$Haz))

tpred <- c(sfgen11$time,sfgen11$time-w)
tpred <- c(0,sort(tpred[tpred>0 & tpred<=10]))
npred <- length(tpred)
Fwgen11 <- Fwgen12 <- rep(NA,npred)
for (i in 1:npred) {
        tp <- tpred[i]
        sfi1 <- sfgen11
        sfi1$haz <- sfi1$haz*exp(iplgen1$coef["LM1"]*g1(tp)+iplgen1$coef["LM2"]*g2(tp))
        sfi1$Haz <- cumsum(sfi1$haz)
        tmp <- approx(sfgen11$time, sfi1$Haz, c(tp,tp+w), method="constant",
                      yleft=0, yright=max(sfi1$Haz))$y
        Fwgen11[i] <- 1-exp(-diff(tmp))
        sfi2 <- sfgen12
        sfi2$haz <- sfi2$haz*exp(iplgen1$coef["LM1"]*g1(tp)+iplgen1$coef["LM2"]*g2(tp))
        sfi2$Haz <- cumsum(sfi2$haz)
        tmp <- approx(sfgen12$time, sfi2$Haz, c(tp,tp+w), method="constant",
                      yleft=0, yright=max(sfi2$Haz))$y
        Fwgen12[i] <- 1-exp(-diff(tmp))
}

## Dynamic predictions for super learner landmark-fixed
# This will have 4 different predictions
sfcomb1 <- survfit(iplcomb1,newdata=data.frame(comb1=as.vector(outer(cl*clin.qs, gl*gen.qs, "+")),LM1=0,LM2=0),censor=FALSE)
sfcomb11 <- data.frame(time=sfcomb1[1]$time,surv=sfcomb1[1]$surv)
sfcomb11$Haz <- -log(sfcomb11$surv)
sfcomb11$haz <- diff(c(0,sfcomb11$Haz))
sfcomb12 <- data.frame(time=sfcomb1[2]$time,surv=sfcomb1[2]$surv)
sfcomb12$Haz <- -log(sfcomb12$surv)
sfcomb12$haz <- diff(c(0,sfcomb12$Haz))
sfcomb13 <- data.frame(time=sfcomb1[3]$time,surv=sfcomb1[3]$surv)
sfcomb13$Haz <- -log(sfcomb13$surv)
sfcomb13$haz <- diff(c(0,sfcomb13$Haz))
sfcomb14 <- data.frame(time=sfcomb1[4]$time,surv=sfcomb1[4]$surv)
sfcomb14$Haz <- -log(sfcomb14$surv)
sfcomb14$haz <- diff(c(0,sfcomb14$Haz))

tpred <- c(sfcomb11$time,sfcomb11$time-w)
tpred <- c(0,sort(tpred[tpred>0 & tpred<=10]))
npred <- length(tpred)
Fwcomb11 <- Fwcomb12 <- Fwcomb13 <- Fwcomb14 <-rep(NA,npred)
for (i in 1:npred) {
        tp <- tpred[i]
        sfi1 <- sfcomb11
        sfi1$haz <- sfi1$haz*exp(iplcomb1$coef["LM1"]*g1(tp)+iplcomb1$coef["LM2"]*g2(tp))
        sfi1$Haz <- cumsum(sfi1$haz)
        tmp <- approx(sfcomb11$time, sfi1$Haz, c(tp,tp+w), method="constant",
                      yleft=0, yright=max(sfi1$Haz))$y
        Fwcomb11[i] <- 1-exp(-diff(tmp))
        sfi2 <- sfcomb12
        sfi2$haz <- sfi2$haz*exp(iplcomb1$coef["LM1"]*g1(tp)+iplcomb1$coef["LM2"]*g2(tp))
        sfi2$Haz <- cumsum(sfi2$haz)
        tmp <- approx(sfcomb12$time, sfi2$Haz, c(tp,tp+w), method="constant",
                      yleft=0, yright=max(sfi2$Haz))$y
        Fwcomb12[i] <- 1-exp(-diff(tmp))
        sfi3 <- sfcomb13
        sfi3$haz <- sfi3$haz*exp(iplcomb1$coef["LM1"]*g1(tp)+iplcomb1$coef["LM2"]*g2(tp))
        sfi3$Haz <- cumsum(sfi3$haz)
        tmp <- approx(sfcomb13$time, sfi3$Haz, c(tp,tp+w), method="constant",
                      yleft=0, yright=max(sfi3$Haz))$y
        Fwcomb13[i] <- 1-exp(-diff(tmp))
        sfi4 <- sfcomb14
        sfi4$haz <- sfi4$haz*exp(iplcomb1$coef["LM1"]*g1(tp)+iplcomb1$coef["LM2"]*g2(tp))
        sfi4$Haz <- cumsum(sfi4$haz)
        tmp <- approx(sfcomb14$time, sfi4$Haz, c(tp,tp+w), method="constant",
                      yleft=0, yright=max(sfi4$Haz))$y
        Fwcomb14[i] <- 1-exp(-diff(tmp))
}

## Dynamic predictions for clinical landmark-dependent
sfclin2 <- survfit(iplclin2,newdata=data.frame(clin1=clin.qs,clin2=0,LM1=0,LM2=0),censor=FALSE)
sfclin21 <- data.frame(time=sfclin2[1]$time,surv=sfclin2[1]$surv)
sfclin21$Haz <- -log(sfclin21$surv)
sfclin21$haz <- diff(c(0,sfclin21$Haz))
sfclin22 <- data.frame(time=sfclin2[2]$time,surv=sfclin2[2]$surv)
sfclin22$Haz <- -log(sfclin22$surv)
sfclin22$haz <- diff(c(0,sfclin22$Haz))

tpred <- c(sfclin21$time,sfclin21$time-w)
tpred <- c(0,sort(tpred[tpred>0 & tpred<=10]))
npred <- length(tpred)
Fwclin21 <- Fwclin22 <- rep(NA,npred)
for (i in 1:npred) {
        tp <- tpred[i]
        sfi1 <- sfclin21
        sfi1$haz <- sfi1$haz*exp(iplclin2$coef["clin2"]*clin.qs[1]*f2(tp)+iplclin2$coef["LM1"]*g1(tp)+iplclin2$coef["LM2"]*g2(tp))
        sfi1$Haz <- cumsum(sfi1$haz)
        tmp <- approx(sfclin21$time, sfi1$Haz, c(tp,tp+w), method="constant",
                      yleft=0, yright=max(sfi1$Haz))$y
        Fwclin21[i] <- 1-exp(-diff(tmp))
        sfi2 <- sfclin22
        sfi2$haz <- sfi2$haz*exp(iplclin2$coef["clin2"]*clin.qs[2]*f2(tp)+iplclin2$coef["LM1"]*g1(tp)+iplclin2$coef["LM2"]*g2(tp))
        sfi2$Haz <- cumsum(sfi2$haz)
        tmp <- approx(sfclin22$time, sfi2$Haz, c(tp,tp+w), method="constant",
                      yleft=0, yright=max(sfi2$Haz))$y
        Fwclin22[i] <- 1-exp(-diff(tmp))
}

## Dynamic predictions for genomic landmark-dependent
sfgen2 <- survfit(iplgen2,newdata=data.frame(gen1=gen.qs,gen2=0,LM1=0,LM2=0),censor=FALSE)
sfgen21 <- data.frame(time=sfgen2[1]$time,surv=sfgen2[1]$surv)
sfgen21$Haz <- -log(sfgen21$surv)
sfgen21$haz <- diff(c(0,sfgen21$Haz))
sfgen22 <- data.frame(time=sfgen2[2]$time,surv=sfgen2[2]$surv)
sfgen22$Haz <- -log(sfgen22$surv)
sfgen22$haz <- diff(c(0,sfgen22$Haz))

tpred <- c(sfgen21$time,sfgen21$time-w)
tpred <- c(0,sort(tpred[tpred>0 & tpred<=10]))
npred <- length(tpred)
Fwgen21 <- Fwgen22 <- rep(NA,npred)
for (i in 1:npred) {
        tp <- tpred[i]
        sfi1 <- sfgen21
        sfi1$haz <- sfi1$haz*exp(iplgen2$coef["gen2"]*gen.qs[1]*f2(tp)+iplgen2$coef["LM1"]*g1(tp)+iplgen2$coef["LM2"]*g2(tp))
        sfi1$Haz <- cumsum(sfi1$haz)
        tmp <- approx(sfgen21$time, sfi1$Haz, c(tp,tp+w), method="constant",
                      yleft=0, yright=max(sfi1$Haz))$y
        Fwgen21[i] <- 1-exp(-diff(tmp))
        sfi2 <- sfgen22
        sfi2$haz <- sfi2$haz*exp(iplgen2$coef["gen2"]*gen.qs[2]*f2(tp)+iplgen2$coef["LM1"]*g1(tp)+iplgen2$coef["LM2"]*g2(tp))
        sfi2$Haz <- cumsum(sfi2$haz)
        tmp <- approx(sfgen22$time, sfi2$Haz, c(tp,tp+w), method="constant",
                      yleft=0, yright=max(sfi2$Haz))$y
        Fwgen22[i] <- 1-exp(-diff(tmp))
}

## Dynamic predictions for super learner landmark-dependent
## This will have 4 different predictions
combs <- as.vector(outer(cl*clin.qs, gl*gen.qs, "+"))
sfcomb2 <- survfit(iplcomb2,newdata=data.frame(comb1=combs,comb2=0,LM1=0,LM2=0),censor=FALSE)
sfcomb21 <- data.frame(time=sfcomb2[1]$time,surv=sfcomb2[1]$surv)
sfcomb21$Haz <- -log(sfcomb21$surv)
sfcomb21$haz <- diff(c(0,sfcomb21$Haz))
sfcomb22 <- data.frame(time=sfcomb2[2]$time,surv=sfcomb2[2]$surv)
sfcomb22$Haz <- -log(sfcomb22$surv)
sfcomb22$haz <- diff(c(0,sfcomb22$Haz))
sfcomb23 <- data.frame(time=sfcomb2[3]$time,surv=sfcomb2[3]$surv)
sfcomb23$Haz <- -log(sfcomb23$surv)
sfcomb23$haz <- diff(c(0,sfcomb23$Haz))
sfcomb24 <- data.frame(time=sfcomb2[4]$time,surv=sfcomb2[4]$surv)
sfcomb24$Haz <- -log(sfcomb24$surv)
sfcomb24$haz <- diff(c(0,sfcomb24$Haz))

tpred <- c(sfcomb21$time,sfcomb21$time-w)
tpred <- c(0,sort(tpred[tpred>0 & tpred<=10]))
npred <- length(tpred)
Fwcomb21 <- Fwcomb22 <- Fwcomb23 <- Fwcomb24 <-rep(NA,npred)
for (i in 1:npred) {
        tp <- tpred[i]
        sfi1 <- sfcomb21
        sfi1$haz <- sfi1$haz*exp(iplcomb2$coef["comb2"]*combs[1]*f2(tp)+iplcomb2$coef["LM1"]*g1(tp)+iplcomb2$coef["LM2"]*g2(tp))
        sfi1$Haz <- cumsum(sfi1$haz)
        tmp <- approx(sfcomb21$time, sfi1$Haz, c(tp,tp+w), method="constant",
                      yleft=0, yright=max(sfi1$Haz))$y
        Fwcomb21[i] <- 1-exp(-diff(tmp))
        sfi2 <- sfcomb22
        sfi2$haz <- sfi2$haz*exp(iplcomb2$coef["comb2"]*combs[2]*f2(tp)+iplcomb2$coef["LM1"]*g1(tp)+iplcomb2$coef["LM2"]*g2(tp))
        sfi2$Haz <- cumsum(sfi2$haz)
        tmp <- approx(sfcomb22$time, sfi2$Haz, c(tp,tp+w), method="constant",
                      yleft=0, yright=max(sfi2$Haz))$y
        Fwcomb22[i] <- 1-exp(-diff(tmp))
        sfi3 <- sfcomb23
        sfi3$haz <- sfi3$haz*exp(iplcomb2$coef["comb2"]*combs[3]*f2(tp)+iplcomb2$coef["LM1"]*g1(tp)+iplcomb2$coef["LM2"]*g2(tp))
        sfi3$Haz <- cumsum(sfi3$haz)
        tmp <- approx(sfcomb23$time, sfi3$Haz, c(tp,tp+w), method="constant",
                      yleft=0, yright=max(sfi3$Haz))$y
        Fwcomb23[i] <- 1-exp(-diff(tmp))
        sfi4 <- sfcomb24
        sfi4$haz <- sfi4$haz*exp(iplcomb2$coef["comb2"]*combs[4]*f2(tp)+iplcomb2$coef["LM1"]*g1(tp)+iplcomb2$coef["LM2"]*g2(tp))
        sfi4$Haz <- cumsum(sfi4$haz)
        tmp <- approx(sfcomb24$time, sfi4$Haz, c(tp,tp+w), method="constant",
                      yleft=0, yright=max(sfi4$Haz))$y
        Fwcomb24[i] <- 1-exp(-diff(tmp))
}

par(mfrow=c(1,2))
plot(tpred,Fwclin11,type="s",lwd=3,ylim=c(0, 1),
     xlab="Prediction time (years)",ylab="Probability")
lines(tpred,Fwclin12,type="s",lwd=3,lty=2)
lines(tpred,Fwclin21+0.02,type="s",lwd=3,col="#770060")
lines(tpred,Fwclin22+0.02,type="s",lwd=3,col="#770060",lty=2)
legend("topleft",
       c("High risk, landmark fixed","High risk, landmark dependent",
         "Low risk, landmark fixed","Low risk, landmark dependent"),
       lwd=3,lty=c(2,2,1,1),
       col=c(1,"#770060",1,"#770060"),bty="n", cex=0.8)
title(main="Clinical PI")

plot(tpred,Fwgen11,type="s",lwd=3,ylim=c(0,0.91),
     xlab="Prediction time (years)",ylab="Probability")
lines(tpred,Fwgen12,type="s",lwd=3,lty=2)
lines(tpred,Fwgen21+0.02,type="s",lwd=3,col="#770060")
lines(tpred,Fwgen22+0.02,type="s",lwd=3,col="#770060",lty=2)
legend("topleft",
       c("High risk, landmark fixed","High risk, landmark dependent",
         "Low risk, landmark fixed","Low risk, landmark dependent"),
       lwd=3,lty=c(2,2,1,1),
       col=c(1,"#770060",1,"#770060"),bty="n", cex=0.8)
title(main="Genetic PI")


plot(tpred,Fwcomb11,type="s",lwd=3,ylim=c(0,1),
     xlab="Prediction time (years)",ylab="Probability")
lines(tpred,Fwcomb12,type="s",lwd=3,lty=2)
lines(tpred,Fwcomb13+0.05,type="s",lwd=3,lty=3)
lines(tpred,Fwcomb14,type="s",lwd=3,lty=4)
lines(tpred,Fwcomb21-0.03,type="s",lwd=3,col="#770060")
lines(tpred,Fwcomb22-0.03,type="s",lwd=3,col="#770060",lty=2)
lines(tpred,Fwcomb23+0.04,type="s",lwd=3,col="#770060",lty=3)
lines(tpred,Fwcomb24,type="s",lwd=3,col="#770060",lty=4)
legend("topleft",
       c("High risk clinical, high risk genetic",
         "High risk clinical, low risk genetic",
         "Low risk clinical, high risk genetic",
         "Low risk clinical, Low risk genetic"),
       lwd=3,lty=c(4,2,3,1),bty="n")
legend("bottomright",c("Landmark fixed","Landmark dependent"),
       lwd=3,col=c(1,"#770060"),bty="n")
title(main="Clinico-Genotypic PI")

###############################################################################
###############################################################################
### Kullback-Leibler and Brier dynamic (fixed width w = 5) prediction error
### curves for the landmark supermodels without and with landmark interactions
###############################################################################
###############################################################################

###Kullback-Leibler dynamic prediction error curves
# Have to call pe for each landmark time point for each of the models
LMs <- seq(0,7,by=0.025)
w <- 1.5
pes <- matrix(NA,length(LMs),9)
pes[,1] <- LMs
for (j in 1:length(LMs)) {
        LM <- LMs[j]
        # deb(LM, method="cat")
        ## Make predictions
        
        # Use landmark data set for estimating conditional censoring distr
        LMdataj <- cutLM(data=xzmat,outcome=list(time="tstop",status="status"),
                         LM=LM,horizon=LM+w,covs=list(fixed=c("zbeta","xbeta","xzbeta"),timedep=NULL))
        # Explicitly code interactions of treatment with LM and put in LMdataj
        LMdataj$clin1 <- LMdataj$zbeta*f1(LMdataj$LM)
        LMdataj$clin2 <- LMdataj$zbeta*f2(LMdataj$LM)
        LMdataj$gen1 <- LMdataj$xbeta*f1(LMdataj$LM)
        LMdataj$gen2 <- LMdataj$xbeta*f2(LMdataj$LM)
        LMdataj$comb1 <- LMdataj$xzbeta*f1(LMdataj$LM)
        LMdataj$comb2 <- LMdataj$xzbeta*f2(LMdataj$LM)
        LMdataj$LM1 <- g1(LMdataj$LM)
        LMdataj$LM2 <- g2(LMdataj$LM)
        
        ni <- nrow(LMdataj)
        #deb(ni, method="cat")
        sfcens <- survfit(Surv(tstop,status==0)~1,data=LMdataj)
        tcens <- sfcens$time
        censmat <- matrix(sfcens$surv,length(sfcens$surv),ni)
        # clin1
        H0 <- diff(evalstep(time=bhclin1$time,stepf=bhclin1$hazard,newtime=c(LM,LM+w),subst=0))
        B <- iplclin1$coef
        HR <- exp(B["clin1"] * LMdataj$clin1 + B["LM1"]*g1(LM) + B["LM2"]*g2(LM))
        pes[j,2] <- pe(time=LMdataj$tstop, status=LMdataj$status,
                       tsurv=LM+w, survmat=matrix(exp(-HR*H0),1,ni),
                       tcens=tcens, censmat=censmat, tout=LM+w-0.00001)$Err
        # clin2
        H0 <- diff(evalstep(time=bhclin2$time,stepf=bhclin2$hazard,newtime=c(LM,LM+w),subst=0))
        B <- iplclin2$coef
        HR <- exp(B["clin1"]*LMdataj$clin1 + B["clin2"]*LMdataj$clin2*f2(LM) + B["LM1"]*g1(LM) + B["LM2"]*g2(LM))
        pes[j,3] <- pe(time=LMdataj$tstop, status=LMdataj$status,
                       tsurv=LM+w, survmat=matrix(exp(-HR*H0),1,ni),
                       tcens=tcens, censmat=censmat, tout=LM+w-0.00001)$Err
        # gen1
        H0 <- diff(evalstep(time=bhgen1$time,stepf=bhgen1$hazard,newtime=c(LM,LM+w),subst=0))
        B <- iplgen1$coef
        HR <- exp(B["gen1"] * LMdataj$gen1 + B["LM1"]*g1(LM) + B["LM2"]*g2(LM))
        pes[j,4] <- pe(time=LMdataj$tstop, status=LMdataj$status,
                       tsurv=LM+w, survmat=matrix(exp(-HR*H0),1,ni),
                       tcens=tcens, censmat=censmat, tout=LM+w-0.00001)$Err
        # gen2
        H0 <- diff(evalstep(time=bhgen2$time,stepf=bhgen2$hazard,newtime=c(LM,LM+w),subst=0))
        B <- iplgen2$coef
        HR <- exp(B["gen1"]*LMdataj$gen1 + B["gen2"]*LMdataj$gen2*f2(LM) + B["LM1"]*g1(LM) + B["LM2"]*g2(LM))
        pes[j,5] <- pe(time=LMdataj$tstop, status=LMdataj$status,
                       tsurv=LM+w, survmat=matrix(exp(-HR*H0),1,ni),
                       tcens=tcens, censmat=censmat, tout=LM+w-0.00001)$Err
        # comb1
        H0 <- diff(evalstep(time=bhcomb1$time,stepf=bhcomb1$hazard,newtime=c(LM,LM+w),subst=0))
        B <- iplcomb1$coef
        HR <- exp(B["comb1"] * LMdataj$comb1 + B["LM1"]*g1(LM) + B["LM2"]*g2(LM))
        pes[j,6] <- pe(time=LMdataj$tstop, status=LMdataj$status,
                       tsurv=LM+w, survmat=matrix(exp(-HR*H0),1,ni),
                       tcens=tcens, censmat=censmat, tout=LM+w-0.00001)$Err
        # comb2
        H0 <- diff(evalstep(time=bhcomb2$time,stepf=bhcomb2$hazard,newtime=c(LM,LM+w),subst=0))
        B <- iplcomb2$coef
        HR <- exp(B["comb1"]*LMdataj$comb1 + B["comb2"]*LMdataj$comb2*f2(LM) + B["LM1"]*g1(LM) + B["LM2"]*g2(LM))
        pes[j,7] <- pe(time=LMdataj$tstop, status=LMdataj$status,
                       tsurv=LM+w, survmat=matrix(exp(-HR*H0),1,ni),
                       tcens=tcens, censmat=censmat, tout=LM+w-0.00001)$Err
        # cov1
        H0 <- diff(evalstep(time=bhcov1$time,stepf=bhcov1$hazard,newtime=c(LM,LM+w),subst=0))
        B <- iplcov1$coef
        HR <- exp(B["clin1"]*LMdataj$clin1 + B["gen1"]*LMdataj$gen1 + B["LM1"]*g1(LM) + B["LM2"]*g2(LM))
        pes[j,8] <- pe(time=LMdataj$tstop, status=LMdataj$status,
                       tsurv=LM+w, survmat=matrix(exp(-HR*H0),1,ni),
                       tcens=tcens, censmat=censmat, tout=LM+w-0.00001)$Err
        # cov2
        H0 <- diff(evalstep(time=bhcov2$time,stepf=bhcov2$hazard,newtime=c(LM,LM+w),subst=0))
        B <- iplcov2$coef
        HR <- exp(B["clin1"]*LMdataj$clin1 + B["clin2"]*LMdataj$clin2*f2(LM) +
                          + B["gen1"]*LMdataj$gen1 + B["gen2"]*LMdataj$gen2*f2(LM) + B["LM1"]*g1(LM) + B["LM2"]*g2(LM))
        pes[j,9] <- pe(time=LMdataj$tstop, status=LMdataj$status,
                       tsurv=LM+w, survmat=matrix(exp(-HR*H0),1,ni),
                       tcens=tcens, censmat=censmat, tout=LM+w-0.00001)$Err
}
pes <- as.data.frame(pes)
names(pes) <- c("time","clin1","clin2","gen1","gen2","comb1","comb2","cov1","cov2")

#### Compute Kullback-Leibler PE
# The null model is easier, just a call to pewcox
KLw0 <- pewcox(Surv(tstop,status)~1, Surv(tstop,status==0)~1, width=w,data=xzmat, FUN = "KL")
pes$null <- evalstep(time=KLw0$time,stepf=KLw0$Err,newtime=LMs,subst=0)
## With and Without landmark interactions plotted together

oldpar <- par(no.readonly=TRUE) # save graphical parameters
vdist <- hdist <- 0.4
layout(matrix(1:4, 2, 2, byrow=TRUE),widths=c(10,10),heights=c(10,10))
par(mar= c(vdist, 4, 3, hdist))

plot(pes$time,pes$null,type="s",lwd=3,lty=1,col=8,ylim=c(0.1,0.6),xaxt="n",
     xlab="Time (years)",ylab="Prediction error")
lines(pes$time,pes$clin1,type="s",lwd=3,lty=1,col=1)
lines(pes$time,pes$gen1,type="s",lwd=3,lty=2,col=1)
lines(pes$time,pes$comb1,type="s",lwd=3,lty=3,col=1)
lines(pes$time,pes$clin2-0.02,type="s",lwd=3,lty=1,col="#996060")
lines(pes$time,pes$gen2-0.02,type="s",lwd=3,lty=2,col="#996060")
lines(pes$time,pes$comb2-0.02,type="s",lwd=3,lty=3,col="#996060")
legend("topleft",c("Null model","clinical model","genetic model","clinico-genotypic model"),
       lwd=3,lty=c(1,1,2,3,1,2,3),col=c(8,1,1,1, "#996060"),bty="n")
legend("bottomright",c("Landmark fixed","Landmark dependent"),
       lwd=3,col=c(1,"#996060"),bty="n", cex=0.8)
text(5,0.35,"Kullback-Leibler",adj=1, cex = 0.8)
axis(3); box()

### Plot Kullback-Leibler PE reduction
pesr <- pes
pesr[,2:9] <- (pesr[,10]-pesr[,2:9])/pesr[,10]

## With and Without landmark interactions plottet together
par(mar= c(vdist, hdist, 3, 4))
plot(pesr$time,abs(pesr$clin1),type="s",lwd=3,lty=1,col=1,ylim=c(0,0.6),yaxt="n",xaxt="n",
     xlab="Time (years)",ylab="Prediction error reduction")
lines(pesr$time,abs(pesr$gen1),type="s",lwd=3,lty=2,col=1)
lines(pesr$time,abs(pesr$comb1),type="s",lwd=3,lty=3,col=1)
lines(pesr$time,abs(pesr$clin2)+0.01,type="s",lwd=3,lty=1,col="#996060")
lines(pesr$time,abs(pesr$gen2)+0.01,type="s",lwd=3,lty=2,col="#996060")
lines(pesr$time,abs(pesr$comb2)+0.01,type="s",lwd=3,lty=3,col="#996060")
legend("topright",c("clinical model","genetic model","clinico-genotypic model"),
       lwd=3,lty=c(1,2,3),col=c(1,"#646060",1),bty="n")
legend("bottomleft",c("Landmark fixed","Landmark dependent"),
       lwd=3,col=c(1,"#996060"),bty="n", cex=0.8)
text(3,0.2,"Kullback-Leibler",adj=1, cex=0.8)
axis(3);axis(4); box()

##############################################
#### Compute Brier Prediction error curves###
###########################################
# The null model is easier, just a call to pewcox
KLw0 <- pewcox(Surv(tstop,status)~1, Surv(tstop,1-status)~1, width=w,data=xzmat, FUN = "Brier")
pes$null <- evalstep(time=KLw0$time,stepf=KLw0$Err,newtime=LMs,subst=0)
pes$null=pes$null+max( pes[,2])

## With and Without landmark interactions plotted together
par(mar= c(5, 4, vdist, hdist))
plot(pes$time,pes$null,type="s",lwd=3,lty=1,col=8,xaxt="n",
     ylim=c(0,0.8),xlab="Time (years)",ylab="Prediction error")
lines(pes$time,pes$clin1,type="s",lwd=3,lty=1,col=1)
lines(pes$time,pes$gen1,type="s",lwd=3,lty=2,col=1)
lines(pes$time,pes$comb1,type="s",lwd=3,lty=3,col=1)
lines(pes$time,pes$clin2-0.02,type="s",lwd=3,lty=1,col="#996060")
lines(pes$time,pes$gen2-0.02,type="s",lwd=3,lty=2,col="#996060")
lines(pes$time,pes$comb2-0.02,type="s",lwd=3,lty=3,col="#996060")
legend("bottomright",c("Null model","clinical model","genetic model","clinico-genotypic model"),
       lwd=3,lty=c(1,1,2,3,1,2,3),col=c(8,1,1,1, "#996060"),bty="n")
legend("topleft",c("Landmark fixed","Landmark dependent"),
       lwd=3,col=c(1,"#996060"),bty="n", cex=0.8)
text(5,0.6,"Brier",adj=1)
axis(1);box()

### Plot Brier PE reduction
pesr <- pes
pesr[,2:9] <- (pesr[,10]-pesr[,2:9])/pesr[,10]

## With and Without landmark interactions plotted together
par(mar= c(5, hdist, vdist, 4))
plot(pesr$time, abs(pesr$clin1),type="s",lwd=3,lty=1,col=1, ylim=c(0.15,0.8),yaxt="n",xaxt="n",
     xlab="Time (years)",ylab="Prediction error reduction")
lines(pesr$time,abs(pesr$gen1),type="s",lwd=3,lty=2,col=1)
lines(pesr$time,abs(pesr$comb1),type="s",lwd=3,lty=3,col=1)
lines(pesr$time,abs(pesr$clin2)+0.02,type="s",lwd=3,lty=1,col="#996060")
lines(pesr$time,abs(pesr$gen2)+0.02,type="s",lwd=3,lty=2,col="#996060")
lines(pesr$time,abs(pesr$comb2)+0.02,type="s",lwd=3,lty=3,col="#996060")
legend("topright",c("clinical model","genetic model","clinico-genotypic model"),
       lwd=3,lty=c(1,2,3),col=c(1,"#646060",1),bty="n")
legend("bottomleft",c("Landmark fixed","Landmark dependent"),
       lwd=3,col=c(1,"#996060"),bty="n", cex=0.8)
text(4,0.6,"Brier",adj=1)
axis(1); axis(4)
layout(matrix(1, 1, 1))
par(oldpar) # reset graphical parameters



#################################################################
################################################################
###  Compute CV landmark-specific Genetic  PI using  
###  using landmark-specific Genetic  predictors   
###  Here time-dependent Cox-Lasso is performed on different
### Landmark data sets with a prediction window of width w=5years
###################################################################
###################################################################

### Run Lasso regression on landmark data sets, calculate cross-validated
### prognostic index for each individual in the landmark data set,
### and append that to xzmat data frame

LMs <- seq(0,5,by=1)
w <- 5
PI.LM<-matrix(NA, nrow = nrow(xzmat), ncol = 1)
PI.LM<-as.data.frame(PI.LM)
PI.LM$V1=xzmat$auslongid

idx<-which(colnames(xzmat)%in%c("auslongid", "tstop", "status", "time", "Latexp"))
xzmatplus <- xzmat[,c(idx, 49:247)]
nr<-dim(xzmatplus)[2]
n<-nrow(xzmat)
cindex.LM<-rep(NA, 6)
for (i in 1:length(LMs)) {
        LM <- LMs[i]
        cat("\n\nLandmark time point:",LM,"\n\n")
        LMdata <- cutLM(data=xzmatplus,outcome=list(time="tstop",status="status"),
                        LM=LM,horizon=LM+w,covs=list(fixed=names(xzmatplus)[-(2:3)],timedep=NULL))
        deb(dim(LMdata), method="cat")
        print(LMdata[1:6,1:8])
        optlam1.LM <- optL1(Surv(LM,tstop,status), LMdata[,6:nr], data=LMdata, fold = 10)
        keep<-names(coefficients(optlam1.LM$fullfit))
        ##second step selection
        LMdata$logt<-logt<-log(LMdata$time+0.05)
        ints<-c(paste(keep, "Latexp", sep = ":"), paste(keep, "logt", sep = ":"))
        vars<-as.character(c(keep,ints))
        f.form<-as.formula(paste("Surv(tstop, status)", paste(vars, collapse=" + "), sep=" ~"))
        cfixed <- coxph(f.form,ties="breslow", robust = TRUE, data=LMdata) ##standard cox with time-varying coefficients
        xvar.df<-data.frame(summary(cfixed)[8])
        xvar.df<-data.frame(summary(cfixed)[8])
        xvar1<-which(xvar.df[,3]>=1&xvar.df[,4]>1)
        xvar2<-which(xvar.df[,3]<=1&xvar.df[,4]<=1)
        xvars<-rownames(xvar.df[c(xvar1, xvar2),])
        ints<-which(!xvars%in%keep)
        ###third step selection
        vars<-c(keep, xvars[ints], "logt", "Latexp", "cluster(auslongid)")
        f.form<-as.formula(paste("Surv(LM,tstop, status)", paste(vars, collapse=" + "), sep=" ~"))
        cfixed <- coxph(f.form,ties="breslow", robust = TRUE, data=LMdata)
        xvar.df<-data.frame(summary(cfixed)[8])
        xvar.df<-data.frame(summary(cfixed)[8])
        xvar1<-which(xvar.df[,3]>=1&xvar.df[,4]>1)
        xvar2<-which(xvar.df[,3]<=1&xvar.df[,4]<=1)
        xvars<-rownames(xvar.df[c(xvar1, xvar2),])
        ##Final model
        f.form<-as.formula(paste("Surv(tstop, status)", paste(xvars, collapse=" + "), sep=" ~"))
        #(cfixed <- coxph(f.form,ties="breslow", robust = TRUE, data = xzmat))
        date();
        xbeta.LM <- CVcindex(f.form,data = xzmat, matrix=TRUE)  ## you must specify matrix=TRUE to get CVPIs to get doubly-robust PI
        cindex.LM[LM+1]<-xbeta.LM$cindex                          ## get Lanmark-specific c-index
        xbeta.LM<- diag(xbeta.LM$matrix)
        xbeta.LM<-xbeta.LM-mean(xbeta.LM) ##standzardize the cvpi
        nm <- paste("xbeta.LM",LM,sep="")
        add <- rep(NA,n)
        add[LMdata$auslongid%in%PI.LM$V1] <- xbeta.LM
        PI.LM <- cbind(PI.LM,add)
        names(PI.LM)[ncol(PI.LM)] <- nm
}
cindex.LM
#[1] 0.8350026 0.8515524 0.8657602 0.8619747 0.8319389 0.8448032

##append Landmark specific PIs to data

xzmat<-cbind(xzmat, PI.LM[,-1])

# Write result to tab-delimited text file (uncomment next line)
# write.table(xzmat,file="xzmat.txt",sep="\t",row.names=FALSE,col.names=TRUE, quote = F)

# xzmat <- read.table("xzmat.txt",header=TRUE,sep="\t")

# This did not make it in the book, was replaced by table with correlations
pairs(PI.LM[,c(2:7)],pch=20,cex=0.5)

## Standard deviations and correlations of cross-validated
## landmark specific genetic lasso predictors


round(apply(PI.LM[,c(2:7)],2,sd,na.rm=TRUE),3)
round(cor(PI.LM[,c(2:7)],use="pairwise.complete.obs"),3)

######################################################
######################################################
###  Compute CV landmark-specific clinical PI using  #
###  using landmark-specific clinical predictors     #
######################################################
######################################################
LMs <- seq(0,5,by=1)
w <- 5
zvars<-c( "Grp5bmi","Grp4bmi","Grp3bmi","Grp2bmi", "Sex", "antidep",
          "ChDiet","Imunized","Edu3cat", "Edu2cat", "LowInc",
          "agecat2","agecat3","agecat4" , "RelCont", "LatExp","vitDInt", "Event", "chedss")
cPI.LM<-matrix(NA, nrow = nrow(xzmat), ncol = 1)
cPI.LM<-as.data.frame(cPI.LM)
cPI.LM$V1=xzmat$auslongid
##make vector to hold Landmark-specific clinical c-index
cindex.LM<-rep(NA, 6)
for (i in 1:length(LMs)) {
        LM <- LMs[i]
        cat("\n\nLandmark time point:",LM,"\n\n")
        LMdata <- cutLM(data=xzmat,outcome=list(time="tstop",status="status"),
                        LM=LM,horizon=LM+w, covs=list(fixed=c(zvars,"Obstime", "auslongid", "time"),timedep=NULL))
        deb(dim(LMdata))
        print(LMdata[1:6,1:8])
        ####start by selecting clinical predictors per Landmark
        LMdata$logt=log(LMdata$time+0.05)
        v<-paste(zvars, "logt",sep=":")
        vars<-as.character(c(zvars, v))
        form<-as.formula(paste("Surv(LM, tstop, status)",paste(vars, collapse = " + "), sep=" ~ "))
        ctime <- coxph(form, data = LMdata) ##standard cox with time-varying coefficients
        zvar.df<-data.frame(summary(ctime)[8])
        zvar1<-which(zvar.df[,3]>=1&zvar.df[,4]>1)
        zvar2<-which(zvar.df[,3]<=1&zvar.df[,4]<=1)
        lvars<-rownames(zvar.df[c(zvar1, zvar2),])
        lvars<-c(zvars, lvars[which(!lvars%in%zvars)])
        ##now fit time-dependent cox-model
        f.form <- as.formula(paste("Surv(tstop,status)",paste( c(lvars, "cluster(auslongid)"),collapse = " + "),sep="~"))
        date()
        zbeta.LM <- CVcindex(f.form,data = xzmat, matrix=TRUE)  #you must specify matrix=TRUE to get CVPIs
        date()
        cindex.LM[LM+1]<-zbeta.LM$cindex  ##get Lanmark-specific c-index
        zbeta.LM<- diag(zbeta.LM$matrix)
        zbeta.LM<-zbeta.LM-mean(zbeta.LM) ##standzardize the cvpi
        ###Append the Landmark-specific zbeta.LM to the clinical dataframe cPI.LM
        add <- rep(NA,n)
        nm <- paste("zbeta.LM",LM,sep="")
        add[LMdata$auslongid%in%cPI.LM$V1] <- zbeta.LM
        cPI.LM <- cbind(cPI.LM,add)
        names(cPI.LM)[ncol(cPI.LM)] <- nm
}

cindex.LM
#[1] 0.8691391 0.8685503 0.8682399 0.8404440 0.8457392 0.8471356

##append Landmark specific PIs to data

xzmat<-cbind(xzmat, cPI.LM[,-1])

# Write result to tab-delimited text file (uncomment next line)
# write.table(xzmat,file="xzmat.txt",sep="\t",row.names=FALSE,col.names=TRUE, quote = F)

# xzmat <- read.table("xzmat.txt",header=TRUE,sep="\t")

# This did not make it in the book, was replaced by table with correlations
pairs(cPI.LM[,2:7],pch=20,cex=0.5)


###  Standard deviations and correlations of cross-validated
### landmark-specific clinical td-cox regression

round(apply(cPI.LM[,2:7],2,sd,na.rm=TRUE),3)
round(cor(cPI.LM[,2:7],use="pairwise.complete.obs"),3)



#########################################################################################
#########################################################################################
### Time-Dependent Cox-Regression Models on Landmark Data sets for  Separate; Combined###
### and clinico-genotypic prognstic indices                                      #   ####
#########################################################################################
#########################################################################################

###Potential starting clinical protectors allowed to enter the clinical model
zvars<-c("Grp5bmi","Grp4bmi","Grp3bmi","Grp2bmi", 
         "agecat4" ,"agecat3","agecat2",
         "antidep", "Smoker","Sex",
         "ChDiet","Imunized","XrayExp",
         "vitDInt","baseVitD","RelCont",
         "Status","Intchedss",
         "ChEmp","Edu3cat", "Edu2cat",
         "FullEmp","DisabEmp", "HomeEmp",
         "Season","ChSunExp","SmokeExp","LatExp", 
         "HighInc" ,"MedInc","LowInc" ,"logt")
indz<-which(colnames(xzmat)%in%zvars)

###create clinical and genetic matrices
xmat <- xzmat[,c(1,5,6, 7,48, 51:249, 250,251, 254)]
zmat <- xzmat[,c(1,5,6, indz)]
zind<-which(colnames(zmat)%in%zvars)
allmat<-cbind(xmat, zmat)
##create landmark dataset
LMs <- seq(0,5,by=1)
w <- 5

##create matrix to hold results for genetic PIs
cxresults<-matrix(NA, nrow = 6, ncol = 4)
colnames(cxresults)<-c("beta", "se", "chisq", "aic")
rownames(cxresults)<-c("LM0", "LM1", "LM2", "LM3", "LM4", "LM5")
LMxresults<-cxresults

##create matrix to hold results for clinical PIs
czresults<-matrix(NA, nrow = 6, ncol = 4)
colnames(czresults)<-c("beta", "se", "chisq", "aic")
rownames(czresults)<-c("LM0", "LM1", "LM2", "LM3", "LM4", "LM5")
LMzresults<-czresults

##create matrix to hold results for clinico-genotypic PIs
cxzresults<-matrix(NA, nrow = 6, ncol = 4)
colnames(cxzresults)<-c("beta", "se", "chisq", "aic")
rownames(cxzresults)<-c("LM0", "LM1", "LM2", "LM3", "LM4", "LM5")
LMxzresults<-cxzresults

for (i in 1:length(LMs)) {
        LM <- LMs[i]
        cat("\n\nLandmark time point:",LM,"\n\n")
        LMdata <- cutLM(data=allmat,outcome=list(time="tstop",status="status"),
                        LM=LM,horizon=LM+w,covs=list(fixed=c(names(xmat)[-(2:3)],names(zmat[,zind])),timedep=NULL))
        deb(dim(LMdata), method="cat")
        ni<-204
        nr<-dim(LMdata)[1]
        print(LMdata[1:6,1:8])
        optlam1.LM <- optL1(Surv(LM,tstop,status), LMdata[,6:ni], minlambda1 = 15, maxlambda1 = 23,data=LMdata, fold = 10)
        keep<-names(coefficients(optlam1.LM$fullfit))
        
        ##first step to select interacting features
        ints<-c(paste(keep, "Latexp", sep = ":"), paste(keep, "logt", sep = ":"))
        vars<-as.character(c(keep,ints, "logt", "Latexp"))
        LMdata$logt<-logt<-log(LMdata$Obstime+1)
        form<-as.formula(paste("Surv(LM, tstop, status)",paste(vars, collapse = " + "), sep=" ~ "))
        ctime <- coxph(form, data = LMdata) ##standard cox with time-varying coefficients
        xvar.df<-data.frame(summary(ctime)[8])
        xvar1<-which(xvar.df[,3]>=1&xvar.df[,4]>1)
        xvar2<-which(xvar.df[,3]<=1&xvar.df[,4]<=1)
        xvars<-rownames(xvar.df[c(xvar1, xvar2),])
        ints<-which(!xvars%in%keep)
        xvars[ints]
        
        ###Create a new matrix with interactions terms
        vars<-unique(c(keep,xvars[ints], "logt", "Latexp"))
        form<-as.formula(paste("Surv(LM,tstop, status)", paste(vars, collapse=" + "), sep=" ~"))
        gmat <- model.matrix(form, data = LMdata)[,-1] 
        dmat<-cbind(LMdata[,c("tstop", "status")], gmat)
        
        ##Now fit time-dependent cox-model
        xbeta.LM <- CVcindex(form,data = LMdata, matrix=TRUE)  #you must specify matrix=TRUE to get CVPIs
        cindex<-xbeta.LM$cindex  ##get Lanmark-specific c-index
        xbeta.LM<- diag(xbeta.LM$matrix)
        LMdata$xbeta.LM<-xbeta.LM<-xbeta.LM-mean(xbeta.LM) ##standzardize the cvpi
        ind<-which(colnames(LMdata)=="xbeta.LM")
        colnames(LMdata)[ind]=paste("xbeta.LM",LM,sep="")
        
        ##Final Cox-Regression on Landmark
        form  <- as.formula(paste("Surv(LM,tstop,status)",paste(c(paste("xbeta.LM",LM,sep=""), "cluster(auslongid)"),collapse=" + "),sep="~"))
        cx  <- coxph(Surv(LM,tstop,status)~xbeta+cluster(auslongid), data=LMdata,robust=TRUE, method="breslow")
        cxLM  <- coxph(form, data=LMdata,robust=TRUE,  method="breslow")
        cat("\n\nLandmark-specific PI on Landmark:",LM,"\n\n")
        
        ##Vector to collect all genetic Results
        cxresults[LM+1,]=round(c(summary(cx)$coef[,c(1,4)],2*diff(cx$loglik), AIC(cx)),2)
        LMxresults[LM+1,]=round(c(summary(cxLM)$coef[,c(1,4)],2*diff(cxLM$loglik), AIC(cxLM)),2)
        
        ##################################
        ##Landmark-specific clinical PI##
        #################################
        
        vars<-colnames(zmat[,zind])
        form<-as.formula(paste("Surv(LM, tstop, status)",paste(vars, collapse = " + "), sep=" ~ "))
        ctime <- coxph(form, data = LMdata) ##standard cox with time-varying coefficients
        zvar.df<-data.frame(summary(ctime)[8])
        zvar1<-which(zvar.df[,3]>=1&zvar.df[,4]>1)
        zvar2<-which(zvar.df[,3]<=1&zvar.df[,4]<=1)
        lvars<-rownames(zvar.df[c(zvar1, zvar2),])
        
        ###Always adjust for BMI, Sex and Age even if they were not selected
        zvars<-c( "Grp5bmi","Grp4bmi","Grp3bmi","Grp2bmi", 
                  "agecat4" ,"agecat3","agecat2","Sex")
        lvars<-unique(c(zvars, lvars[which(!lvars%in%zvars)]))
        lvars<-lvars[lvars!="logt"]
        
        ##Compute the landmark-specific clinical CVPIs
        form <- as.formula(paste("Surv(LM, tstop,status)",paste( c(lvars),collapse = " + "),sep="~"))
        zbeta.LM <- CVcindex(form,data = LMdata, matrix=TRUE)  #you must specify matrix=TRUE to get CVPIs
        cindex<-zbeta.LM$cindex  ##get Lanmark-specific c-index
        zbeta.LM<- diag(zbeta.LM$matrix)
        LMdata$zbeta.LM<-zbeta.LM<-zbeta.LM-mean(zbeta.LM) ##standzardize the cvpi
        ind<-which(colnames(LMdata)=="zbeta.LM")
        colnames(LMdata)[ind]=paste("zbeta.LM",LM,sep="")
        
        ##Final Cox-Regression on Landmark-Specific clinical CVPIs and overall CVPIs
        form <- as.formula(paste("Surv(LM,tstop,status)",paste(c(paste("zbeta.LM",LM,sep=""), "cluster(auslongid)"),collapse=" + "),sep="~"))
        cz  <- coxph(Surv(LM,tstop,status)~zbeta+cluster(auslongid), data=LMdata,robust=TRUE, method="breslow")
        czLM <- coxph(form, data=LMdata,robust = TRUE, method="breslow")
        czresults[LM+1,]=round(c(summary(cz)$coef[,c(1,4)],2*diff(cz$loglik), AIC(cz)),2)
        LMzresults[LM+1,]=round(c(summary(czLM)$coef[,c(1,4)],2*diff(czLM$loglik), AIC(czLM)),2)
        
        ################################
        ##compute clinico-genotypic CVPI
        ################################
        
        ##Landmark-specific clinico-genotypic model
        form <- as.formula(paste("Surv(LM,tstop,status)",paste(c(paste(c("zbeta.LM","xbeta.LM"),LM,sep=""), "cluster(auslongid)"),collapse=" + "),sep="~"))
        cxzLM <- coxph(form, data=LMdata,robust = TRUE, method="breslow")
        ###manually compute the clinico-genotypic PI
        LMdata$xzbeta.LM<-xzbeta.LM<-coef(cxzLM)[1]*zbeta.LM+coef(cxzLM)[2]*xbeta.LM ### compute the Super Learner P
        ind<-which(colnames(LMdata)=="xzbeta.LM")
        colnames(LMdata)[ind]=paste("xzbeta.LM",LM,sep="")
        
        ##Landmark-fixed clinico-genotypic model
        cxz <- coxph(Surv(LM,tstop,status)~xzbeta+cluster(auslongid), data=LMdata,robust=TRUE, method="breslow")
        cxzresults[LM+1,]=round(c(summary(cxz)$coef[,c(1,4)],2*diff(cxz$loglik), AIC(cxz)),2)
        
        ##Final super supermodel Cox-Regression on Landmark
        form    <- as.formula(paste("Surv(LM,tstop,status)",paste(c(paste("xzbeta.LM",LM,sep=""), "cluster(auslongid)"),collapse=" + "),sep="~"))
        superLM <- coxph(form, data=LMdata,robust = TRUE, method="breslow")                                             ##landmark-specific super PI
        LMxzresults[LM+1,]=round(c(summary(superLM)$coef[,c(1,4)],2*diff(superLM$loglik), AIC(superLM)),2)
        
}

save(czresults, cxresults, cxzresults,LMzresults, LMxresults, LMxzresults, file="LandmarkSpecific2.RData")

czresults  ##clinical Landmark-fixed
cxresults  ## genetic Landmark-fixed
cxzresults ##clinico-genotypic Landmark-fixed
LMzresults ##clinical Landmark-specific
LMxresults ##genetic Landmark-specific
LMxzresults ##clinico-genetic Landmark-specific



