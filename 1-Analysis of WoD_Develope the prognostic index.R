###############################################################################
###############################################################################
###
### This file contains R code for the analyses presented in the article 
### "Developing a clinical-environmental-genotypic prognostic index for 
###  relapsing-onset multiple sclerosis and clinically isolated syndrome" 
### This code was copied from Dynamic Prediction in Clinical Survival Analysis 
### (CRC Chapman & Hall) by Hans C. van Houwelingen and Hein Putter
### R code written by Hein Putter and modified by Valery Fuh Ngwa
### Correspondence to be sent to     valeryfuh.ngwa@utas.edu.au 
################################################################################


##############################################
### Analysis of EDSS: MS Disability progression
###############################################

####Load required packages and Prepare Data
require(dynpred)
library(survival)
library(lattice)
library(penalized)
library(coxme)
library(mfp)
library(xtable)
##set working directory
setwd("C:/Users/Owner/Documents/Survival_Analysis/DynPred_JointModels/DynPreds_JM/PhD_Study_1/Article_1/Analysis_of_EDSS")

###Read in clinical data
zmat<- read.delim("MasterUltimateLongSurv.txt") 

##Read genotype data
xmat<-read.delim("~/Survival_Analysis/trainG_N272.txt")  #read in 199 SNPs

##Merge clinical and genotype datra

xzmat<-merge(zmat, xmat, "auslongid")

##Define clinical and environmental variables
crude<-c("I(basebmi)", "fp(logt)","strata(edsscat)","strata(cdms)","baseVitD",
         "antidep", "Smoker","Sex","dage","tt(vitd)",
         "Imunized","RelCont",
         "Status","ChEmp","Edu3cat", "Edu2cat",
         "FullEmp","DisabEmp", "HomeEmp",
         "ChSunExp","LatExp","HighInc" ,"MedInc","LowInc", "MarjExp")

#######################################
## Table 1: Clinical Prognostic Factors
#######################################

form<-as.formula(paste("Surv(tstop, status)", paste(c(crude, "cluster(auslongid)"), collapse=" + "), sep=" ~"))
(cfixed <- coxph(form,ties="breslow",tt = function(vitd, t, ...) ifelse(t>vitd, 0,1), robust = TRUE, data = xzmat)) 
xtable(round(summary(cfixed)$coef[,c(1,4,6)],2))  ##Results of table S1

##step-wise variable selection with multi-factorial polynomials
tt = function(vitd, t, ...) ifelse(t<vitd, 0,1)  ## time-dependent effects of vitamin D
xzmat$vidtInt<-with(xzmat, tt(vitd, tstop))
crude<-crude[-10]
form<-as.formula(paste("Surv(tstop, status)", paste(crude, collapse=" + "), sep=" ~"))

##fit multi-factoirial model using mfp functio and set alpha =0.05

mpf=mfp(form, family = cox, data = xzmat, select=0.05, verbose=TRUE)

##Extract final variables selected
zvars<-subset(mpf$trafo,  !is.na(formula)==TRUE)
zvars<-zvars[,1]
adj.vars<-c("dage","Sex","basebmi", "LatExp", "tt(vitd)", "strata(edsscat)","strata(cdms)", "cluster(auslongid)")
core.vars<-unique(c( adj.vars,zvars))
f.form<-as.formula(paste("Surv( tstop, status)", paste(core.vars, collapse=" + "), sep=" ~"))
(cfixed <- coxph(f.form,ties="breslow",tt = function(vitd, t, ...) ifelse(t>vitd, 0,1), robust = TRUE, data = xzmat)) 
xtable(round(summary(cfixed)$coef[,c(1,4,6)],2)) ##generate the table using xtable function


#Refit the Cox models using the selected features 
adj.vars<-c("dage","Sex","basebmi", "LatExp", "vitd", "cluster(auslongid)")
core.vars<-unique(c( adj.vars,zvars))
f.form<-as.formula(paste("Surv(tstop, status)", paste(core.vars, collapse=" + "), sep=" ~"))
(cfixed <- coxph(f.form,ties="breslow", robust = TRUE, data = xzmat)) 

##Compute the Cross-validated clinical-environmental prognostic (CEPI) index using LOOCV
set.seed(25467)
zbeta <- CVcindex(f.form, data = xzmat, matrix=TRUE)  #You have to specify matrix=TRUE to get CVPIs
(c.index<-zbeta$cindex)                               ##get the c-index
zbeta<- diag(zbeta$matrix)
xzmat$zbeta<-zbeta-mean(zbeta)                        ##Standardize the clinical cvpi
mean(xzmat$zbeta)
sqrt(var(xzmat$zbeta))


#####################################################
## Compute the Cross-validated Genetic PI data GPI
####################################################

##Lasso penalised cox regression
optlam  <- optL1(Surv(tstop,status), penalized = xmat,minlambda1 = 0, maxlambda1 = 25, data=xzmat, fold = 20)
lamlass <- optlam$lambda # optimal lambda
keep<-names(coefficients(optlam$fullfit)!=0)  ##
fit <- profL1(Surv(tstop,status), penalized =  xmat[,keep], minlambda1 = 15, maxlambda1 = 35.80009, fold = optlam$fold, data=xzmat, steps=10)

#plot Lasso lamda path
par(mfrow=c(1,2))
plot(fit$lambda, fit$cvl, type="l", xlab="Cross-validated Partial log-likelihood") ## cross-validated Partial log-likelihood
plotpath(fit$fullfit)                ### plot selected coefficients

###Estimate unbiased regression estimates for the Lasso-selected variants
##Use backward selection with multifactorial polynomials

xzmat$logt<-logt<-log(xzmat$time+0.05)                ##for including time-interaction
ints<-c(paste(keep, "Latexp", sep = ":"), paste(keep, "logt", sep = ":"))
vars<-as.character(c(keep,ints, "logt", "age","Sex", "Latexp", "strata(edsscat)","strata(cdms)")) #adjust for age, sex, and stratify EDSS categories and concversion status: CDMS
form<-as.formula(paste("Surv(tstop, status)", paste(vars, collapse=" + "), sep=" ~"))
mpf=mfp(form, family = cox, data = xzmat, select=0.05, verbose=TRUE)

#######################################
## Table S2: Genetic Prognostic Variants#
#######################################

xvars <- subset(mpf$trafo,  !is.na(formula)==TRUE)
xvars <- xvars[,1]
## Create a single stratification variable using interaction()
xzmat$Z  <- with(xzmat, interaction(edsscat, cdms))
xzmat$Z  <- relevel(xzmat$Z, ref = "0.0") ## reference level is 0(cdms) x 0(EDSS) i.e 0.0.
adj.vars<-c("age","Sex", "Latexp","logt","strata(edsscat)","strata(cdms)", "cluster(auslongid)")
core.vars<-unique(c( adj.vars,xvars))
f.form<-as.formula(paste("Surv(tstop, status)", paste(core.vars, collapse=" + "), sep=" ~"))
(cfixed <- coxph(f.form, method="breslow", robust = TRUE, data = xzmat)) 
xtable(round(summary(cfixed)$coef[,c(1,4,6)],2))  ##Results of Table S2

##Compute Cross-Validated Genetic Prognostic Index using the selected variants in Table S2

core.vars <- core.vars[-c(5:6)]
f.form<-as.formula(paste("Surv(tstop, status)", paste(core.vars, collapse=" + "), sep=" ~"))
xbeta <- CVcindex(f.form,data = xzmat, matrix=TRUE)  ## you must specify matrix=TRUE to get CVPIs
(c.index<-xbeta$cindex)                              ## Get C-index
xbeta<- diag(xbeta$matrix)
xzmat$xbeta<-xbeta<-xbeta-mean(xbeta)                ## standzardize the GPI
mean(xzmat$xbeta)
sqrt(var(xzmat$xbeta))

###############################################
## Discrimination and Internal Calibration ###
###############################################

##Cox Regression on CPI and GPI
###############################

summary(cox<-coxph(Surv(tstop, status)~zbeta+xbeta+cluster(auslongid), robust=TRUE, ties= "breslow",data=xzmat))

### Manually compute the CGPI
#############################

xzmat$xzbeta<-xzbeta<-coef(cox)[1]*xzmat$zbeta+coef(cox)[2]*xzmat$xbeta  ### CGPI

## Make 4 groups of equal sizes with CPI and GPI
################################################

xzmat$group.zbeta <- as.numeric(cut(xzmat$zbeta,c(-Inf,quantile(xzmat$zbeta)[2:4],Inf)))
xzmat$group.xbeta <- as.numeric(cut(xzmat$xbeta,c(-Inf,quantile(xzmat$xbeta)[2:4],Inf)))
###write data
write.table(xzmat, file = "xzmat.txt", col.names = T, quote = F, row.names = F, sep = "\t")

###Get discrimination and Calibration Plots
#############################################

oldpar <- par(no.readonly=TRUE) # save graphical parameters
vdist <- hdist <- 0.4
layout(matrix(1:4, 2, 2, byrow=TRUE),widths=c(10,10),heights=c(10,10))
par(mar= c(vdist, 4, 3, hdist))
cx <- coxph(Surv(tstop, status) ~ group.zbeta+cluster(auslongid), data=xzmat,robust=TRUE, method="breslow")
nd <- data.frame(group.zbeta=1)
sf <- survfit(cx,newdata=nd,  censor = FALSE)
plot(sf,lwd=3,conf.int=TRUE,mark.time=FALSE,lty=1,xlab="Time (years)",xaxt="n",ylab="Survival function")
nd <- data.frame(group.zbeta=2)
sf <- survfit(cx,newdata=nd,  censor = FALSE)
lines(sf,lwd=3,mark.time=FALSE,conf.int=TRUE,lty=1,col=2)
nd <- data.frame(group.zbeta=3)
sf <- survfit(cx,newdata=nd,  censor = FALSE)
lines(sf,lwd=3,mark.time=FALSE,conf.int=TRUE,lty=1,col=3)
nd <- data.frame(group.zbeta=4)
sf <- survfit(cx,newdata=nd, censor = FALSE)
lines(sf,lwd=3,mark.time=FALSE,conf.int=TRUE,lty=1,col=4)
cxkm <- coxph(Surv(tstop, status) ~logt, data=xzmat,robust=TRUE, method="breslow")
nd <- data.frame(logt=mean(xzmat$logt))
sf <- survfit(cxkm,newdata=nd,censor = FALSE)
lines(sf,lwd=3,mark.time=FALSE,conf.int=TRUE,lty=1,col=8)
legend("bottomleft",
       c("Kaplan-Meiers","Percentile 0-25","Percentile 25-50","Percentile 50-75","Percentile 75-100"),
       lwd=3,col=c(8,1:4),lty=c(1,1,2,3,4),bty="n", cex=0.8)
axis(3); box()

# zbeta
par(mar= c(vdist, hdist, 3, 4))
cx <- coxph(Surv(tstop, status) ~ group.xbeta+cluster(auslongid), data=xzmat,robust=TRUE, method="breslow")
nd <- data.frame(group.xbeta=1)
sf <- survfit(cx,newdata=nd, times=5, censor = FALSE)
plot(sf,lwd=3,conf.int=TRUE,mark.time=FALSE,lty=1,xlab="Time (years)",yaxt="n", xaxt="n",ylab="Survival function")
nd <- data.frame(group.xbeta=2)
sf <- survfit(cx,newdata=nd,  censor = FALSE)
lines(sf,lwd=3,mark.time=FALSE,conf.int=TRUE,lty=1,col=2)
nd <- data.frame(group.xbeta=3)
sf <- survfit(cx,newdata=nd, censor = FALSE)
lines(sf,lwd=3,mark.time=FALSE,conf.int=FALSE,lty=1,col=3)
nd <- data.frame(group.xbeta=4)
sf <- survfit(cx,newdata=nd,  censor = FALSE)
lines(sf,lwd=3,mark.time=FALSE,conf.int=TRUE,lty=1,col=4)
cxkm <- coxph(Surv(tstop, status) ~logt, data=xzmat,robust=TRUE, method="breslow")
nd <- data.frame(logt=mean(xzmat$logt))
sf <- survfit(cxkm,newdata=nd, censor = FALSE)
lines(sf,lwd=3,mark.time=FALSE,conf.int=TRUE,lty=1,col=8)
legend("bottomleft",
       c("Kaplan-Meiers","Percentile 0-25","Percentile 25-50","Percentile 50-75","Percentile 75-100"),
       lwd=3,col=c(8,1:4),lty=c(1,1,2,3,4),bty="n", cex=0.8)
axis(3);axis(4); box()

### td-cox regression on clinical cvpi
par(mar= c(5, 4, vdist, hdist))
czbeta <- coxph(Surv(tstop, status) ~ zbeta+cluster(auslongid), data=xzmat,robust=TRUE, method="breslow")
nd <- data.frame(zbeta=-2*sd(zbeta))
sf <- survfit(czbeta,newdata=nd)
plot(sf,lwd=3,conf.int=TRUE,mark.time=FALSE,lty=1,xlab="Time (years)",ylab="Survival function", xaxt="n")
nd <- data.frame(zbeta=-1*sd(zbeta))
sf <- survfit(czbeta,newdata=nd)
lines(sf,lwd=3,conf.int=TRUE,lty=1,mark.time=FALSE, col=2)
nd <- data.frame(zbeta=0)
sf <- survfit(czbeta,newdata=nd)
lines(sf,lwd=3,conf.int=TRUE,lty=1,mark.time=FALSE, col=3)
lines(survfit(Surv(tstop, status)~1,data=xzmat),mark.time=FALSE,lwd=2,lty=1,col=8) ##Kaplan meiers with subgroups
nd <- data.frame(zbeta=sd(zbeta))
sf <- survfit(czbeta,newdata=nd)
lines(sf,lwd=3,conf.int=TRUE,lty=1,mark.time=FALSE, col=4)
nd <- data.frame(zbeta=2*sd(zbeta))
sf <- survfit(czbeta,newdata=nd)
lines(sf,lwd=3,conf.int=TRUE,lty=1,mark.time=FALSE, col=5)
legend("bottomleft",c("Survival 2 SD below mean","Survival 1 SD below mean",
                      "Survival Average  PI", "Kaplan-Meiers","Survival 1 SD above mean","Survival 2 SD above mean"),
       lwd=3,lty=c(2,3,1,1,4,5),col=c(1:3,8,4:5), bty="n", cex=0.8)
axis(1);box()
### td-cox regression on genetic cvpi
par(mar= c(5, hdist, vdist, 4))
cxbeta <- coxph(Surv(tstop, status) ~ xbeta+cluster(auslongid), data=xzmat,robust=TRUE, method="breslow")
nd <- data.frame(xbeta=-2*sd(xbeta))
sf <- survfit(cxbeta,newdata=nd)
plot(sf,lwd=3,conf.int=TRUE,mark.time=FALSE,lty=1,xlab="Time (years)",ylab="Survival function",yaxt="n",xaxt="n")
nd <- data.frame(xbeta=-1*sd(xbeta))
sf <- survfit(cxbeta,newdata=nd)
lines(sf,lwd=3,conf.int=TRUE,lty=1,mark.time=FALSE, col=2)
nd <- data.frame(xbeta=0)
sf <- survfit(cxbeta,newdata=nd)
lines(sf,lwd=3,conf.int=TRUE,lty=1,mark.time=FALSE, col=3)
cxkm <- coxph(Surv(tstop, status) ~ logt, data=xzmat,robust=TRUE, method="breslow")
nd <- data.frame(logt=mean(logt))
sf <- survfit(cxkm,newdata=nd)
lines(sf,lwd=3,mark.time=FALSE,conf.int=TRUE,lty=1,col=8)
nd <- data.frame(xbeta=sd(xbeta))
sf <- survfit(cxbeta,newdata=nd)
lines(sf,lwd=3,conf.int=TRUE,lty=1,mark.time=FALSE, col=4)
nd <- data.frame(xbeta=2*sd(xbeta))
sf <- survfit(cxbeta,newdata=nd)
lines(sf,lwd=3,conf.int=TRUE,lty=1,mark.time=FALSE, col=5)
legend("bottomleft",c("Survival 2 SD below mean","Survival 1 SD below mean",
                      "Survival Average  PI", "Kaplan-Meiers","Survival 1 SD above mean","Survival 2 SD above mean"),
       lwd=3,lty=c(2,3,1,1,4,5),col=c(1:3,8,4:5), bty="n", cex = 0.8)

axis(1); axis(4)
layout(matrix(1, 1, 1))
par(oldpar) # reset graphical parameters
####################################################


####################################################
### Produce results for Risk Deciles of the CVPIs  #
####################################################

###make decile groups for genetic, clinical and clinico-genotypic PIs
xzmat$group.xbeta<- as.numeric(cut(xbeta,c(-Inf,quantile(xbeta, seq(0,100, 10)/100)[2:11],Inf)))
xzmat$group.zbeta<- as.numeric(cut(zbeta,c(-Inf,quantile(zbeta, seq(0,100, 10)/100)[2:11],Inf)))
xzmat$group.xzbeta<- as.numeric(cut(xzbeta,c(-Inf,quantile(xzbeta, seq(0,100, 10)/100)[2:11],Inf)))
###convert the groups to factors variables
xzmat$group.xbeta<-as.factor(xzmat$group.xbeta)
xzmat$group.zbeta<-as.factor(xzmat$group.zbeta)
xzmat$group.xzbeta<-as.factor(xzmat$group.xzbeta)

### td-cox regress on deciles of PIs
DF <- within(xzmat, group.zbeta <- relevel(group.zbeta, ref = 5))
DF <- within(xzmat, group.xbeta <- relevel(group.xbeta, ref = 5))
DF <- within(xzmat, group.xzbeta <- relevel(group.xzbeta, ref = 5))

##Reference group is the 10 percentile
AIC(coxa<-coxph(Surv(tstop, status)~group.zbeta+cluster(auslongid), robust=TRUE,
                ties= "breslow",data=xzmat))
AIC(coxb<-coxph(Surv(tstop, status)~group.xbeta+cluster(auslongid), robust=TRUE,
                ties= "breslow",data=xzmat))
AIC(coxc<-coxph(Surv(tstop, status)~group.xzbeta+cluster(auslongid), robust=TRUE,
                ties= "breslow",data=xzmat))
coefz=as.numeric(summary(coxa)$coef[,1])
coefx=as.numeric(summary(coxb)$coef[,1])
coefc=as.numeric(summary(coxc)$coef[,1])
##manually create risk deciles, with 10% being the reference group
x<-c( 20, 30, 40, 50, 60, 70, 80,90, 100)
px=data.frame(x=x,yx=c(coefx))
pz=data.frame(x=x,yz=c(coefz))
pxz=data.frame(x=x,yxz=c(coefc))

###Plot the log-odds ratios for each group
par(mfrow=c(1,2))
rbeta<-c(px[,2], pz[,2], pxz[,2])

plot(pxz, pch=25,ylim=c(range(rbeta)),
     col="#667837",lwd=3, lty=3,cex=1.5, ylab="log-odds ratios",main="Ref=10%", xlab="Deciles of Prognostic Indices (%)" )
points(px,pch=23, cex = 1.5, lwd=3, col="#998200")
points(pz,pch=24, cex = 1.5, lwd=3, col="#770060")
lines(px, lty=3, lwd=3,col="#998200" )
lines(pxz, lty=4, lwd=3, col="#667837")
lines(pz,  lty=5, lwd=3, col="#770060")
legend("topleft", legend = c("Clinico-Genotypic PI", "Genetic PI",
                             "clinical PI"), lty=c(4,3,5), bty="n", lwd=3, pch=c(25,23,24),
       col = c("#667837","#998200","#770060"))

##Reference group is the 50 percentile
AIC(coxa<-coxph(Surv(tstop, status)~group.zbeta+cluster(auslongid), robust=TRUE,
                ties= "breslow",data=DF))
AIC(coxb<-coxph(Surv(tstop, status)~group.xbeta+cluster(auslongid), robust=TRUE,
                ties= "breslow",data=DF))
AIC(coxc<-coxph(Surv(tstop, status)~group.xzbeta+cluster(auslongid), robust=TRUE,
                ties= "breslow",data=DF))
coefz=as.numeric(summary(coxa)$coef[,1])
coefx=as.numeric(summary(coxb)$coef[,1])
coefc=as.numeric(summary(coxc)$coef[,1])

########################################
### manually create risk deciles,      #
### with 10% being the reference group #
########################################

x<-c(10, 20, 30, 40, 60, 70, 80,90, 100)
px=data.frame(x=x,yx=c(coefx))
pz=data.frame(x=x,yz=c(coefz))
pxz=data.frame(x=x,yxz=c(coefc))

rbeta<-c(px[,2], pz[,2], pxz[,2])

plot(pxz, pch=25,ylim=c(range(rbeta)), 
     col="#667837",lwd=3, lty=3,cex=1.5, ylab="log-odds ratios",main="Ref=50%", xlab="Deciles of Prognostic Indices (%)" )
points(px,pch=23, cex = 1.5, lwd=3, col="#998200")
points(pz,pch=24, cex = 1.5, lwd=3, col="#770060")
lines(px, lty=3, lwd=3,col="#998200" )
lines(pxz, lty=4, lwd=3, col="#667837")
lines(pz,  lty=5, lwd=3, col="#770060")
legend("topleft", legend = c("Clinico-Genotypic PI", "Genetic PI",
                             "clinical PI"), lty=c(4,3,5), bty="n", lwd=3, pch=c(25,23,24),
       col = c("#667837","#998200","#770060"))

#################################
#"#667837", "#998200", "#770060"#
#################################

#####################################################
### Make Scatter Plots of Genetic and Clinical PIs  #
#####################################################

### Ridge and Lasso PI
#######################

oldpar <- par(no.readonly=TRUE) # save graphical parameters
layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), widths=c(2,1), heights=c(1,2))
# Scatterplot
par(mar=c(5,5,1,1))
plot(zbeta,xbeta,xlim=c(min(zbeta),max(zbeta)),ylim=c(-2.0,3.8),xlab="",ylab="")
mtext("Clinical PI",side=1,line=3,cex=1.5,font=2)
mtext("Genetic PI",side=2,line=3,cex=1.5,font=2)
# Histogram of ridge
par(mar=c(0,5,2,1))
yhist <- hist(zbeta, breaks=seq(min(zbeta)-1,max(zbeta)+1,0.09), plot=FALSE)
top <- max(yhist$counts)
barplot(yhist$counts, axes=TRUE, ylim=c(0, top), space=0, horiz=FALSE, col = "blue")
leg <- legend("topright","",bty="n")$text
text(leg$x,leg$y,
     paste("Mean =",format(round(mean(zbeta),2),nsmall=2),"\nSD =",round(sd(zbeta),2)),
     adj=c(1,1))
par(mar=c(5,0,1,2))
yhist <- hist(xbeta, breaks=seq(min(xbeta)-1,max(xbeta)+1,0.09), plot=FALSE)
top <- max(yhist$counts)
barplot(yhist$counts, axes=TRUE, xlim=c(0, top), space=0, horiz=TRUE, col = "red")
leg <- legend("topright","",bty="n")$text
text(leg$x,leg$y,
     paste("Mean =",format(round(mean(xbeta),2),nsmall=2),"\nSD =",round(sd(xbeta),2)),
     adj=c(1,1))
par(oldpar) # reset graphical parameters

########################################
## Plot Individual trajectories for PIs#
########################################
library(lattice)

xyplot(biomrk1~ Obstime|Sex, data = xLong, type = "l", groups = auslongid, lty = 2,lwd=2,
       col = "#762783", xlab = "Follow-up Times (Obstime)", lines = TRUE, points = FALSE,
       ylab = expression("Genetic Prognostic Index"))

xyplot(zbeta~ Obstime| Sex, data = xzmat, type = "l", groups = auslongid, lty = 2,lwd=2,
       col = "#762783", xlab = "Follow-up Times (Obstime)", lines = TRUE, points = FALSE,
       ylab = expression("Clinical Prognostic Index"))

xyplot(xzbeta~ Obstime| Sex, data = xzmat, type = "l", groups = auslongid, lty = 2,lwd=2,
       col = "#762783", xlab = "Follow-up Times (Obstime)", lines = TRUE, points = FALSE,
       ylab = expression("Clinico-Genotypic Prognostic Index"))

##Make Survival Average profiles of first four subjects
library(ggplot2)
library(reshape2)
df= subset(xzmat, auslongid%in%c("2400100", "2400300", "2400500", "2400600"))
df=data.frame(Genetic=df$xbeta, Clinical=df$zbeta, id=df$auslongid, Obstime=df$Obstime, "Clinico-Genotypic"=df$xzbeta)
df2 <- melt(data = df, id.vars = c("id", "Obstime"))
ggplot(data=df2, aes(x = Obstime, y = value, colour=variable)) +
  geom_line(size = 1.5, lty=2) +
  facet_wrap(~id)+
  ylab(expression("Prognostic Indices"))+
  xlab("Follow-up Time (years)")

##########################################################
# Compute AIC; Summary stats,  and Prediction R2 for CVPIs#
##########################################################

##clinical PI
AIC(coxz <- coxph(Surv(tstop, status)~zbeta+cluster(auslongid),robust = TRUE, ties = "breslow", data = xzmat)) ##final cox-clinical model
summary(coxz) ##final clinical results 
logtest <- -2 * (coxz$loglik[1] - coxz$loglik[2])
(rsq = 1 - exp(-logtest/coxz$nevent))


## Genetic PI
AIC(coxx<- coxph(Surv(tstop,status) ~ xbeta+ cluster(auslongid),
                 robust = T, data=xzmat, ties = "breslow",))
summary(coxx) ##final genetic results
logtest <- -2 * (coxx$loglik[1] - coxx$loglik[2])
(rsq = 1 - exp(-logtest/coxx$nevent))


## clinico-genetic PI

AIC(coxcomb<- coxph(Surv( tstop,status) ~xbeta+zbeta+ cluster(auslongid),
                    robust = T, data=xzmat, method="breslow"))
summary(coxcomb)        ##final clinico-genetic results
logtest <- -2 * (coxcomb$loglik[1] - coxcomb$loglik[2])
(rsq = 1 - exp(-logtest/coxcomb$nevent))


###############################################################################
###############################################################################
### Kullback-Leibler and Brier prediction error curves (left) and 
### prediction error reduction curves (right) for the null model (Kaplan-Meier)
### and for the three models
###############################################################################
###############################################################################

#################################
## Start with Kullback-Leibler ##
#################################

KL0 <- pecox(Surv( tstop, status)~1,Surv( tstop, 1-status)~1,data=xzmat,FUN = "KL")
KL0$Err[is.na(KL0$Err)] <- 0
KL1 <- pecox(Surv( tstop, status)~zbeta+cluster(auslongid),Surv( tstop, 1-status)~1,data=xzmat,FUN = "KL")
KL1$Err[is.na(KL1$Err)] <- 0
KL2 <- pecox(Surv( tstop, status)~xbeta+cluster(auslongid),Surv( tstop,  1-status)~1,data=xzmat,FUN = "KL")
KL2$Err[is.na(KL2$Err)] <- 0
KL3 <- pecox(Surv( tstop, status)~xzbeta+cluster(auslongid),Surv( tstop, 1-status)~1,data=xzmat,FUN = "KL")
KL3$Err[is.na(KL3$Err)] <- 0

##compute pred error reduction for brier
KL <- data.frame(time=KL0$time,Err0=KL0$Err,Err1=KL1$Err,Err2=KL2$Err,Err3=KL3$Err)
KL$ErrRed1 <- (KL$Err0-KL$Err1)/KL$Err0
KL$ErrRed2 <- (KL$Err0-KL$Err2)/KL$Err0
KL$ErrRed3 <- (KL$Err0-KL$Err3)/KL$Err0

###Compute PE for Brier
#########################
BR0 <- pecox(Surv(tstop, status)~1,Surv(tstop, 1-status)~1,data=xzmat,FUN = "Brier")
BR0$Err[is.na(BR0$Err)] <- 0
BR1 <- pecox(Surv(tstop, status)~zbeta,Surv(tstop, 1-status)~1,data=xzmat,FUN = "Brier")
BR1$Err[is.na(BR1$Err)] <- 0
BR2 <- pecox(Surv(tstop, status)~xbeta,Surv(tstop, 1-status)~1,data=xzmat,FUN = "Brier")
BR2$Err[is.na(BR2$Err)] <- 0
BR3 <- pecox(Surv(tstop, status)~xzbeta,Surv(tstop, 1-status)~1,data=xzmat,FUN = "Brier")
BR3$Err[is.na(BR3$Err)] <- 0

##compute pred error reduction for brier
BR <- data.frame(time=BR0$time,Err0=BR0$Err,Err1=BR1$Err,Err2=BR2$Err,Err3=BR3$Err)
BR$ErrRed1 <- (BR$Err0-BR$Err1)/BR$Err0
BR$ErrRed2 <- (BR$Err0-BR$Err2)/BR$Err0
BR$ErrRed3 <- (BR$Err0-BR$Err3)/BR$Err0
BR$ErrRed1[is.na(BR$ErrRed1)==TRUE]=0
### Plot prediction error curves for kullback-leibler

oldpar <- par(no.readonly=TRUE) # save graphical parameters
vdist <- hdist <- 0.4
layout(matrix(1:4, 2, 2, byrow=TRUE),widths=c(10,10),heights=c(10,10))
par(mar= c(vdist, 4, 3, hdist))

plot(KL0$time,KL0$Err,type="s",lwd=3,lty=3,xlim=c(0,15),col=8,ylim=c(0, 0.7),xaxt="n",
     xlab="Time (years)",ylab="Prediction error")
lines(KL1$time,KL1$Err,type="s",lwd=3,lty=3, col="#667837")
lines(KL2$time,KL2$Err,type="s",lwd=3,lty=3, col="#998200")
lines(KL3$time,KL3$Err,type="s",lwd=3,lty=3, col="#770060")
### Plot lines for Brier##
lines(BR1$time,BR1$Err,type="s",lwd=3,lty=1, col="#667837")
lines(BR2$time,BR2$Err,type="s",lwd=3,lty=1, col="#998200")
lines(BR3$time,BR3$Err,type="s",lwd=3,lty=1, col="#770060")
legend("topleft",c("Null model","clinical model","genetic model","clinico-genotypic model"),
       lwd=3,lty=c(1,1,1,1),col=c(8,"#667837","#998200","#770060"),bty="n", cex=0.8)
text(10,0.35,"Kullback-Leibler",adj=1)
text(6,0.05,"Brier",adj=1)
axis(3); box()

###Plot Prediction error reduction curves for kullback-leibler
par(mar= c(vdist, hdist, 3, 4))
plot(KL$time,abs(KL$ErrRed1),type="s", lty=3, lwd=3,xlim=c(0,15),ylim=c(0,0.7),col="#667837",xaxt="n",yaxt="n",
     xlab="Time (years)",ylab="Prediction error reduction")
lines(KL$time,KL$ErrRed2,type="s",lty=3, lwd=3,col="#998200")
lines(KL$time,KL$ErrRed3,type="s",lty=3, lwd=3,col="#770060")
### Plot lines curves for Brier error reduction curves
lines(BR$time,BR$ErrRed1,type="s",lty=1, lwd=3,col="#667837")
lines(BR$time,BR$ErrRed2,type="s",lty=1, lwd=3,col="#998200")
lines(BR$time,BR$ErrRed3,type="s",lty=1, lwd=3,col="#770060")
legend("topleft",
       c("clinical model","genetic model","clinico-genotypic model"),
       lwd=3,lty=c(1,1,1),col=c("#667837","#998200","#770060"),bty="n", cex=0.8)
text(8,0.1,"Kullback-Leibler",adj=1)
text(8,0.5,"Brier",adj=1)
axis(3);axis(4); box()

###############################################################################
###  Dynamic prediction error curves (w = 2.5) for Kullback-Leibler
### and Briers
###############################################################################

dynpe1.zbeta <- pewcox(Surv(tstop, status) ~ zbeta,
                       Surv(tstop, 1-status) ~ 1, data = xzmat, width = 2.5, FUN = "Brier")
dynpe0.Brier <- pewcox(Surv(tstop, status) ~ 1,
                       Surv(tstop, 1-status) ~ 1, data = xzmat, width = 2.5, FUN = "Brier")

dynpe1.xbeta <- pewcox(Surv(tstop, status) ~ xbeta,
                       Surv(tstop, 1-status) ~ 1, data = xzmat, width = 2.5, FUN = "Brier")
dynpe1.xzbeta <- pewcox(Surv(tstop, status) ~ xzbeta,
                        Surv(tstop, 1-status) ~ 1, data = xzmat, width = 2.5, FUN = "Brier")

###Kullback Leibler
dynpe1.zbeta.kl <- pewcox(Surv(tstop, status) ~ zbeta,
                          Surv(tstop, 1-status) ~ 1, data = xzmat, width = 2.5, FUN = "KL")
dynpe0.kl <- pewcox(Surv(tstop, status) ~ 1,
                    Surv(tstop, 1-status) ~ 1, data = xzmat, width = 2.5, FUN = "KL")

dynpe1.xbeta.kl <- pewcox(Surv(tstop, status) ~ xbeta,
                          Surv(tstop, 1-status) ~ 1, data = xzmat, width = 2.5, FUN = "KL")
dynpe1.xzbeta.kl <- pewcox(Surv(tstop, status) ~ xzbeta,
                           Surv(tstop, 1-status) ~ 1, data = xzmat, width = 2.5, FUN = "KL")

###create data frames of prediction errors for brier
dynpe.zbeta <- data.frame(time=dynpe1.zbeta$time, Err0=dynpe0.Brier$Err, Err1=dynpe1.zbeta$Err)
dynpe.zbeta$ErrRed <- abs((dynpe.zbeta$Err0-dynpe.zbeta$Err1)/dynpe.zbeta$Err0)

dynpe.xbeta <- data.frame(time=dynpe1.xbeta$time, Err0=dynpe0.Brier$Err, Err1=dynpe1.xbeta$Err)
dynpe.xbeta$ErrRed <- abs((dynpe.xbeta$Err0-dynpe.xbeta$Err1)/dynpe.xbeta$Err0)

dynpe.xzbeta <- data.frame(time=dynpe1.xzbeta$time, Err0=dynpe0.Brier$Err, Err1=dynpe1.xzbeta$Err)
dynpe.xzbeta$ErrRed <- abs((dynpe.xzbeta$Err0-dynpe.xzbeta$Err1)/dynpe.xzbeta$Err0)

###create data frames of prediction errors for Kullback leibler

dynpe.zbeta.kl <- data.frame(time=dynpe1.zbeta.kl$time, Err0=dynpe0.kl$Err, Err1=dynpe1.zbeta.kl$Err)
dynpe.zbeta.kl$ErrRed <- abs((dynpe.zbeta.kl$Err0-dynpe.zbeta.kl$Err1)/dynpe.zbeta.kl$Err0)

dynpe.xbeta.kl <- data.frame(time=dynpe1.xbeta.kl$time, Err0=dynpe0.kl$Err, Err1=dynpe1.xbeta.kl$Err)
dynpe.xbeta.kl$ErrRed <- abs((dynpe.xbeta.kl$Err0-dynpe.xbeta.kl$Err1)/dynpe.xbeta.kl$Err0)

dynpe.xzbeta.kl <- data.frame(time=dynpe1.xzbeta.kl$time, Err0=dynpe0.kl$Err, Err1=dynpe1.xzbeta.kl$Err)
dynpe.xzbeta.kl$ErrRed <- abs((dynpe.xzbeta.kl$Err0-dynpe.xzbeta.kl$Err1)/dynpe.xzbeta.kl$Err0)

# Cut off last part
nt <- length(dynpe.zbeta$time)
dynpe.zbeta <- subset(dynpe.zbeta,time<=15)
dynpe.xbeta <- subset(dynpe.xbeta,time<=15)
dynpe.xzbeta <- subset(dynpe.xzbeta,time<=15)

dynpe.zbeta.kl <- subset(dynpe.zbeta.kl,time<=15)
dynpe.xbeta.kl <- subset(dynpe.xbeta.kl,time<=15)
dynpe.xzbeta.kl <- subset(dynpe.xzbeta.kl,time<=15)

dynpe.zbeta.kl$Err0[is.na(dynpe.zbeta.kl$Err0)] <- 0
dynpe.zbeta.kl$Err1[is.na(dynpe.zbeta.kl$Err1)] <- 0
dynpe.zbeta.kl$ErrRed[is.na(dynpe.zbeta.kl$ErrRed)] <- 0

dynpe.xbeta.kl$Err0[is.na(dynpe.xbeta.kl$Err0)] <- 0
dynpe.xbeta.kl$Err1[is.na(dynpe.xbeta.kl$Err1)] <- 0
dynpe.xbeta.kl$ErrRed[is.na(dynpe.xbeta.kl$ErrRed)] <- 0

dynpe.xzbeta.kl$Err0[is.na(dynpe.xzbeta.kl$Err0)] <- 0
dynpe.xzbeta.kl$Err1[is.na(dynpe.xzbeta.kl$Err1)] <- 0
dynpe.xzbeta.kl$ErrRed[is.na(dynpe.xzbeta.kl$ErrRed)] <- 0


##Plot kulback leibler Dynamic error plots
par(mar= c(5, 4, vdist, hdist))
plot(dynpe.zbeta.kl$time,dynpe.zbeta.kl$Err0,type="s",xlim=range(dynpe.zbeta$time),lty=8, col=8,xaxt="n",
     ylim=c(0,max(cbind(dynpe.zbeta[,c("Err0","Err1")],dynpe.xzbeta.kl[,c("Err0","Err1")]))),
     lwd=3,xlab="Time (years)",ylab="Dynamic prediction error")
lines(dynpe.zbeta.kl$time,dynpe.zbeta.kl$Err1,type="s",lwd=3,lty=3, col="#667837")
lines(dynpe.xbeta.kl$time,dynpe.xbeta.kl$Err1,type="s",lwd=3, lty=3,col="#998200")
lines(dynpe.xzbeta.kl$time,dynpe.xzbeta.kl$Err1,type="s",lwd=3,lty=3, col="#770060")
text(7,0.60,"Kullback-Leibler",adj=1)
###Plot lines for dynamic Brier PE curves
lines(dynpe.zbeta$time,dynpe.zbeta$Err0,type="s",lty=1, col=8,lwd=3)
lines(dynpe.zbeta$time,dynpe.zbeta$Err1,type="s",lwd=3,lty=1, col="#667837")
lines(dynpe.xbeta$time,dynpe.xbeta$Err1,type="s",lwd=3, lty=1,col="#998200")
lines(dynpe.xzbeta$time,dynpe.xzbeta$Err1,type="s",lwd=3,lty=1, col="#770060")
legend("topleft",c("Null model (KM)","clinical model", "genetic model",
                   "clinico-genotypic model"),lwd=3,lty=c(8,1,1,1),
       col=c(8,"#667837","#998200","#770060"),bty="n", cex=0.8)
text(7,0.05,"Brier",adj=1)
axis(1);box()


### Dynamic Error reduction curves
### td-cox regression on genetic cvpi
par(mar= c(5, hdist, vdist, 4))
plot(dynpe.zbeta.kl$time,dynpe.zbeta.kl$ErrRed,type="s",xlim=range(dynpe.zbeta.kl$time),yaxt="n",xaxt="n",
     ylim=c(0,0.45),lwd=3,lty=3, col="#667837",
     xlab="Time (years)",ylab="Dynamic prediction error reduction")
lines(dynpe.xbeta.kl$time,dynpe.xbeta.kl$ErrRed,type="s",lwd=3,lty=3, col="#998200")
lines(dynpe.xzbeta.kl$time,dynpe.xzbeta.kl$ErrRed,type="s",lwd=3, lty=3,col="#770060")
text(6, 0.05,"Kullback-Leibler",adj=1)

lines(dynpe.zbeta$time,dynpe.zbeta$ErrRed,type="s",lwd=3,lty=2, col="#667837")
lines(dynpe.xbeta$time,dynpe.xbeta$ErrRed,type="s",lwd=2, lty=1,col="#998200")
lines(dynpe.xzbeta$time,dynpe.xzbeta$ErrRed,type="s",lwd=3,lty=1, col="#770060")
legend("topleft",c("clinical model", "genetic model",
                   "clinico-genotypic model"),lwd=3,lty=c(1,1,1),
       col=c("#667837","#998200","#770060"),bty="n", cex=0.8)
text(6, 0.6,"Brier",adj=1)
axis(1); axis(4)
layout(matrix(1, 1, 1))
par(oldpar) # reset graphical parameters


##################################################################
### Get Tables of Interval-specific and total prediction errors ##
##################################################################

### Done only for the clinical prognostic index 
###Note: to avoid reducdancy of the codes, change zbeta to xbeta or xzbeta accordingly.
dynpe1.KL1 <- pewcox(Surv(tstop, status) ~ xbeta,Surv(tstop, 1-status) ~ 1, data = xzmat, width = 1, FUN = "KL")
dynpe0.KL1 <- pewcox(Surv(tstop, status) ~ 1,Surv(tstop, 1-status) ~ 1, data = xzmat, width = 1, FUN = "KL")
dynpe1.Brier1 <- pewcox(Surv(tstop, status) ~ xbeta,Surv(tstop, 1-status) ~ 1, data = xzmat, width = 1, FUN = "Brier")
dynpe0.Brier1 <- pewcox(Surv(tstop, status) ~ 1,Surv(tstop, 1-status) ~ 1, data = xzmat, width = 1, FUN = "Brier")

dynpe.Brier1 <- data.frame(time=dynpe1.Brier1$time, Err0=dynpe0.Brier1$Err, Err1=dynpe1.Brier1$Err)
dynpe.Brier1$ErrRed <- abs((dynpe.Brier1$Err0-dynpe.Brier1$Err1)/dynpe.Brier1$Err0)
dynpe.KL1 <- data.frame(time=dynpe1.KL1$time, Err0=dynpe0.KL1$Err, Err1=dynpe1.KL1$Err)
dynpe.KL1$ErrRed <- abs((dynpe.KL1$Err0-dynpe.KL1$Err1)/dynpe.KL1$Err0)

xzmat.km <- survfit(formula = Surv(tstop, status) ~ 1, data = xzmat)
KMstart <- evalstep(xzmat.km$time,xzmat.km$surv,0:12,subst=1)

res <- data.frame(Interval=c("0-1","1-2","2-3","3-4","4-5","5-6","6-7", "7-8", "8-9","9-10", "10-11", "11-12" ),
                  KMstart=KMstart[1:12],
                  Hazard=NA,
                  BrierKM=c(dynpe.Brier1$Err0[1],evalstep(dynpe.Brier1$time,dynpe.Brier1$Err0,1:5)),
                  BrierModel=c(dynpe.Brier1$Err1[1],evalstep(dynpe.Brier1$time,dynpe.Brier1$Err1,1:5)),
                  BrierRed=c(dynpe.Brier1$ErrRed[1],evalstep(dynpe.Brier1$time,dynpe.Brier1$ErrRed,1:5)),
                  KLKM=c(dynpe.KL1$Err0[1],evalstep(dynpe.KL1$time,dynpe.KL1$Err0,1:5)),
                  KLModel=c(dynpe.KL1$Err1[1],evalstep(dynpe.KL1$time,dynpe.KL1$Err1,1:5)),
                  KLRed=c(dynpe.KL1$ErrRed[1],evalstep(dynpe.KL1$time,dynpe.KL1$ErrRed,1:5)))
res$Hazard <- 1 - KMstart[-1]/KMstart[-13]

# Append totals
res <- rbind(res,data.frame(Interval="Total",KMstart=NA,Hazard=NA,
                            BrierKM=sum(res$KMstart*res$BrierKM),
                            BrierModel=sum(res$KMstart*res$BrierModel),
                            BrierRed=1-sum(res$KMstart*res$BrierModel)/sum(res$KMstart*res$BrierKM),
                            KLKM=sum(res$KMstart*res$KLKM),
                            KLModel=sum(res$KMstart*res$KLModel),
                            KLRed=1-sum(res$KMstart*res$KLModel)/sum(res$KMstart*res$KLKM)))
res
# round to three decimals
res[,2:9] <- round(res[,2:9],3)
res

#################################################
### DYanmic AUC(t) for the td-Cox models on PIs##
#################################################
par(mfrow=c(1,3))
AUC(Surv(tstop, status) ~ zbeta  +cluster(auslongid), data  = xzmat)
AUC(Surv(tstop, status) ~ xbeta  +cluster(auslongid), data  = xzmat)
AUC(Surv(tstop, status) ~ xzbeta +cluster(auslongid), data  = xzmat)


############################################################
### Dynamic Harrel's C-index  with w = 2.3 for the Cox model 
###########################################################
##compute and plot Dynamic Harrel's C-index
AUCw.coxdata <- AUCw(Surv(tstop, status) ~ zbeta  +cluster(auslongid), data  = xzmat, width = 2.3)
plot(AUCw.coxdata$time,AUCw.coxdata$AUCw,type="s",lwd=3,lty=2,ylim=c(0.5,0.9),col="#667837",
     xlab="Time (years)",ylab="Dynamic C-index")
AUCw.coxdata <- AUCw(Surv(tstop, status) ~ xbeta  +cluster(auslongid), data  = xzmat, width = 2.3)
lines(AUCw.coxdata$time,AUCw.coxdata$AUCw,type="s", lwd=3,lty=6, col="#998200")
AUCw.coxdata <- AUCw(Surv(tstop, status) ~ xzbeta  +cluster(auslongid), data  = xzmat, width = 2.3)
lines(AUCw.coxdata$time,AUCw.coxdata$AUCw,type="s", lwd=3, lty=3, col="#770060")
legend("topleft",c("clinical model", "genetic model",
                   "clinico-genotypic model"),lwd=3,lty=c(2,6,3),
       col=c("#667837","#998200","#770060"),bty="n")
text(3, 0.54,"C-index(t), w=2.3",adj=1)

############################
### Cross-validated C-index#
############################
## Done only for clinical PI
## Please change accordingly for xbeta, xzbeta

date()
CVcindex(Surv(tstop, status) ~ zbeta +cluster(auslongid),  data  = xzmat)
date()
CVcindex(Surv(tstop, status) ~ xbeta +cluster(auslongid),  data  = xzmat)
date()
CVcindex(Surv(tstop, status) ~ xzbeta+cluster(auslongid),  data  = xzmat)

date()
CVcindex(Surv(tstop, status) ~ zbeta +cluster(auslongid),  data  = xzmat, type="pair")
date();date()
CVcindex(Surv(tstop, status) ~ xbeta +cluster(auslongid),  data  = xzmat, type="pair")
date();date()
CVcindex(Surv(tstop, status) ~ xzbeta+cluster(auslongid),  data  = xzmat, type="pair")
date()

#############################################################
###Partial log-likelihoods with and without cross-validation#
#############################################################

## Done only for clinical PI
## Please change accordingly for xbeta, xzbeta

cfull <- coxph(Surv(tstop, status)   ~ zbeta+cluster(auslongid),method = "breslow", robust = TRUE, data = xzmat)
CVPLKM <- CVPL(Surv(tstop, status)    ~ zbeta+cluster(auslongid), data = xzmat, shrink = 1)
CVPLCox <- CVPL(Surv(tstop, status)   ~ zbeta+cluster(auslongid), data = xzmat, overall=TRUE)
CVPLCVCox <- CVPL(Surv(tstop, status) ~ zbeta+cluster(auslongid), data = xzmat)

cfull$loglik;
diff(cfull$loglik)
CVPLKM
CVPLCox
CVPLCVCox
CVPLCox - CVPLKM
CVPLCVCox - CVPLKM

### Save Data Set for used in next analysis of relapses 
write.table(xzmat, file = "xzmat.txt", col.names = T, quote = F, row.names = F, sep = "\t")

########################################################################

##### Continue to Landmark Predictions #######

