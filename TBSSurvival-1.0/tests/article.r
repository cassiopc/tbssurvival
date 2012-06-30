# TBSSurvival package for R (http://www.R-project.org)
# Copyright (C) 2012 Adriano Polpo, Cassio de Campos, Debajyoti Sinha
#                    Jianchang Lin and Stuart Lipsitz.
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

#######
####### Code to generate the results in the paper:
####### Transform Both Sides Model: A parametric approach
#######

library("TBSSurvival")
## IF YOU DONT HAVE THE LIBRARY, PLEASE DOWNLOAD IT FROM
## http://code.google.com/p/tbssurvival/
## AND INSTALL IT USING:
## R CMD install TBSSurvival_version.tar.gz
## OR INSIDE R:
## install.packages(""TBSSurvival_version.tar.gz",repos=NULL,type="source")
##

## To perform convergence analysis of BE set flag.convergence "TRUE".
## This takes some extra time to run however...
flag.convergence <- TRUE

initial.time <- proc.time()
# Seting initial seed to have the same results that we obtained.
set.seed(1234)

##########################################
## Section 2. Transform-both-side model ##
##########################################

i.fig <- 1 # counter figure number
## Figure 1(a)
## hazard functions with normal N(0,1) error
## plotted for lambda = 0.1, 1 and 2, 
name <- paste("tbs_fig_",i.fig,"a",".eps",sep="")
eps(name,width=5,height=5,paper="special",colormodel="gray")
t <- seq(0.01,10,0.01)
plot(t,htbs(t,lambda=0.1,xi=sqrt(1),beta=1,dist="norm"),type="l",lwd=2,lty=1,col="gray60",
     axes="FALSE",ylim=c(0,1),xlab="",ylab="")
lines(t,htbs(t,lambda=1.0,xi=sqrt(1),beta=1,dist="norm"),type="l",lwd=2,lty=2,col="gray40")
lines(t,htbs(t,lambda=2.0,xi=sqrt(1),beta=1,dist="norm"),type="l",lwd=2,lty=3,col="gray20")
title(xlab="time",ylab="h(t)",main=expression(epsilon ~ textstyle(paste("~",sep="")) ~ Normal(0,sigma^2==1)),cex.lab=1.2)
axis(1,at=seq(0,10,1),lwd=2,lwd.ticks=2,pos=0)
axis(2,at=seq(0,1,0.2),lwd=2,lwd.ticks=2,pos=0)
legend(6,0.95,c(expression(lambda == 0.1),
                expression(lambda == 1.0),
                expression(lambda == 2.0)),
              col=c("gray60","gray40","gray20"),lty=c(1,2,3),cex=1.1,lwd=2,bg="white")
dev.off()

## Figure 1(b)
## hazard functions with double-exp error (b=1)
## and lambda = 0.5, 1, and 2
name <- paste("tbs_fig_",i.fig,"b",".eps",sep="")
eps(name,width=5,height=5,paper="special",colormodel="gray")
t <- seq(0.01,10,0.01)
plot(t,htbs(t,lambda=0.5,xi=sqrt(1),beta=1,dist="doubexp"),type="l",lwd=2,lty=1,col="gray60",
     axes="FALSE",ylim=c(0,1),xlab="",ylab="")
lines(t,htbs(t,lambda=1.0,xi=sqrt(1),beta=1,dist="doubexp"),type="l",lwd=2,lty=2,col="gray40")
lines(t,htbs(t,lambda=2.0,xi=sqrt(1),beta=1,dist="doubexp"),type="l",lwd=2,lty=3,col="gray20")
title(xlab="time",ylab="h(t)",main=expression(epsilon ~ textstyle(paste("~",sep="")) ~ doubexp(b==1)),cex.lab=1.2)
axis(1,at=seq(0,10,1),lwd=2,lwd.ticks=2,pos=0)
axis(2,at=seq(0,1,0.2),lwd=2,lwd.ticks=2,pos=0)
legend(6,0.95,c(expression(lambda == 0.5),
                expression(lambda == 1.0),
                expression(lambda == 2.0)),
              col=c("gray60","gray40","gray20"),lty=c(1,2,3),cex=1.1,lwd=2,bg="white")
dev.off()


################################
##  Section 4 - Experiments   ##
################################


################################
# Section 4.1 - Real data sets #
################################

############################################################################################
# Data                                                                                     #
# Number of Cycles (in Thousands) of Fatigue Life for 67 of 72 Alloy T7987 Specimens that  #
# Failed Before 300 Thousand Cycles.                                                       #
# Meeker & Escobar, pp. 130-131 (1998)                                                     #
############################################################################################
data(alloyT7987)

## Summary statistic of the failure time in thousand cycles for the Alloy T7987 data set.
cat("Time AlloyT7987: Summary statistics\n")
print(summary(alloyT7987$time))


###########
## Part MLE
cat("\n"); cat("Maximum Likelihood Estimation:\n")

## MLE Estimation for each possible error distribution
tbs.mle.norm     <- tbs.survreg.mle(Surv(alloyT7987$time,alloyT7987$delta) ~ 1,dist="norm")
tbs.mle.doubexp  <- tbs.survreg.mle(Surv(alloyT7987$time,alloyT7987$delta) ~ 1,dist="doubexp")
tbs.mle.t        <- tbs.survreg.mle(Surv(alloyT7987$time,alloyT7987$delta) ~ 1,dist="t")
tbs.mle.cauchy   <- tbs.survreg.mle(Surv(alloyT7987$time,alloyT7987$delta) ~ 1,dist="cauchy")
tbs.mle.logistic <- tbs.survreg.mle(Surv(alloyT7987$time,alloyT7987$delta) ~ 1,dist="logistic")

## Building the Table of Example (MLE)
## This is only for ease of inclusion in the paper, as the results are all available
## in the corresponding R variables
text <- paste("\\begin{center}\n",
              "\\begin{table}\n",
              "\\caption{Some quantities of the TBS Model for the Alloy T7987 data set ",
                        "using MLE: AIC, BIC, parameter estimates and their standard ",
                        "deviations in parenthesis.}\n",
              "\\label{table_alloy-mle}\n",
              "\\centering\n",
              "\\begin{tabular}{c|cc|ccc}\n",
              "\\hline\n",
              "Error Distribution & AIC & BIC & $\\widehat{\\lambda}$ & $\\widehat{\\beta_0}$ & $\\widehat{\\xi}$ \\\\ \n",
              "\\hline \n",
              "\\hline \n",
              "Normal        & ",round(tbs.mle.norm$AIC,2),       " & ",round(tbs.mle.norm$BIC,2),        " & ",
                                 round(tbs.mle.norm$par[1],4),    " & ",round(tbs.mle.norm$par[3],4),     " & ",round(tbs.mle.norm$par[2],4),     " \\\\ \n",
              "& & & {\\footnotesize(",round(tbs.mle.norm$std.error[1],4),")} & {\\footnotesize(",round(tbs.mle.norm$std.error[3],4),")} & {\\footnotesize(",round(tbs.mle.norm$std.error[2],5),")} \\\\ \n",
              "DoubExp       & ",round(tbs.mle.doubexp$AIC,2),    " & ",round(tbs.mle.doubexp$BIC,2),     " & ",
                                 round(tbs.mle.doubexp$par[1],4), " & ",round(tbs.mle.doubexp$par[3],4),  " & ",round(tbs.mle.doubexp$par[2],4),  " \\\\ \n",
              "& & & {\\footnotesize(",round(tbs.mle.doubexp$std.error[1],4),")} & {\\footnotesize(",round(tbs.mle.doubexp$std.error[3],4),")} & {\\footnotesize(",round(tbs.mle.doubexp$std.error[2],5),")} \\\\ \n",
              "t-Student     & ",round(tbs.mle.t$AIC,2),          " & ",round(tbs.mle.t$BIC,2),           " & ",
                                 round(tbs.mle.t$par[1],4),       " & ",round(tbs.mle.t$par[3],4),        " & ",round(tbs.mle.t$par[2],4),        " \\\\ \n",
              "& & & {\\footnotesize(",round(tbs.mle.t$std.error[1],4),")} & {\\footnotesize(",round(tbs.mle.t$std.error[3],4),")} & {\\footnotesize(",round(tbs.mle.t$std.error[2],5),")} \\\\ \n",
              "Cauchy        & ",round(tbs.mle.cauchy$AIC,2),     " & ",round(tbs.mle.cauchy$BIC,2),      " & ",
                                 round(tbs.mle.cauchy$par[1],4),  " & ",round(tbs.mle.cauchy$par[3],4),   " & ",round(tbs.mle.cauchy$par[2],4),   " \\\\ \n",
              "& & & {\\footnotesize(",round(tbs.mle.cauchy$std.error[1],4),")} & {\\footnotesize(",round(tbs.mle.cauchy$std.error[3],4),")} & {\\footnotesize(",round(tbs.mle.cauchy$std.error[2],5),")} \\\\ \n",
              "Logistic      & ",round(tbs.mle.logistic$AIC,2),   " & ",round(tbs.mle.logistic$BIC,2),    " & ",
                                 round(tbs.mle.logistic$par[1],4)," & ",round(tbs.mle.logistic$par[3],4), " & ",round(tbs.mle.logistic$par[2],4), " \\\\ \n",
              "& & & {\\footnotesize(",round(tbs.mle.logistic$std.error[1],4),")} & {\\footnotesize(",round(tbs.mle.logistic$std.error[3],4),")} & {\\footnotesize(",round(tbs.mle.logistic$std.error[2],5),")} \\\\ \n",
              "\\hline\n",
              "\\end{tabular}\n",
              "\\end{table}\n",
              "\\end{center}\n",
              sep="")
cat(text, file="table_alloy-mle.tex")
rm(text)

## Estimating the Kaplan-Meier
## This is used to plot the graph with the comparison of the estimations
km <- survfit(formula = Surv(alloyT7987$time, alloyT7987$delta == 1) ~ 1)

i.fig <- i.fig+1
## Figure 2(a) - Reliability functions (MLE)
## Kaplan-Meier is also plotted
name <- paste("tbs_fig_",i.fig,"a",".eps",sep="")
eps(name,width=5,height=5,paper="special",colormodel="gray")
plot(km,ylab="",xlab="",xlim=c(min(alloyT7987$time),max(alloyT7987$time)),conf.int=FALSE,axes=FALSE,lty=1,lwd=1)
t <- seq(min(alloyT7987$time),max(alloyT7987$time),(max(alloyT7987$time)-min(alloyT7987$time)-0.01)/1000)
title(ylab="R(t)",xlab="t: number of cycles (in thousands)",main="Reliability function (MLE)",cex.lab=1.2)
axis(1,at=c(93,100,150,200,250,300),lwd=2,lwd.ticks=2,pos=0)
axis(2,lwd=2,lwd.ticks=2,pos=min(alloyT7987$time))
legend(170,0.95,c("Kaplan-Meier",
                 expression(textstyle(paste("TBS / ",sep="")) ~ epsilon ~ textstyle(paste("~",sep="")) ~ Normal)),
                 col=c(1,"gray20"),lty=c(1,1),cex=1.1,lwd=c(1,2),bg="white")
lines(t,1-ptbs(t,lambda=tbs.mle.norm$par[1],xi=tbs.mle.norm$par[2],
               beta=tbs.mle.norm$par[3],dist=tbs.mle.norm$error.dist),type="l",lwd=2,col="gray20",lty=1)
dev.off()

## Figure 2(b) - Hazard function (MLE)
name <- paste("tbs_fig_",i.fig,"b",".eps",sep="")
eps(name,width=5,height=5,paper="special",colormodel="gray")
t <- seq(0.01,350,(350-0.01)/1000)
 plot(t,htbs(t,l=tbs.mle.norm$par[1],xi=tbs.mle.norm$par[2],
                   beta=tbs.mle.norm$par[3],dist=tbs.mle.norm$error.dist),
      type="l",col="gray20",lwd=2,axes="FALSE",ylim=c(0,0.05),xlab="",ylab="")
title(ylab="h(t)",xlab="t: number of cycles (in thousands)",main="Hazard function (MLE)",cex.lab=1.2)
legend(30,0.048,c(expression(textstyle(paste("TBS / ",sep="")) ~ epsilon ~ textstyle(paste("~",sep="")) ~ Normal)),
                   col=c("gray20","gray10"),lty=c(1,2),cex=1.1,lwd=2,bg="white")
axis(1,lwd=2,lwd.ticks=2,pos=0)
axis(2,lwd=2,lwd.ticks=2,pos=0)
lines(t,htbs(t,lambda=tbs.mle.norm$par[1],xi=tbs.mle.norm$par[2],
             beta=tbs.mle.norm$par[3],dist=tbs.mle.norm$error.dist),type="l",col="gray20",lwd=2)
dev.off()

## Residual summary statistics of TBS model with normal error distribution (MLE).
cat("\n"); cat("Error: Summary statistics\n")
print(summary(tbs.mle.norm$error))


## Quantile Estimation
cat("\n"); cat("Quatile estimates:\n")
cat("Median time: ",round(exp(tbs.mle.norm$par[3]),2),"\n")
cat("95% c.i.: (",round(exp(tbs.mle.norm$par[3]+qnorm(0.025,0,1)*tbs.mle.norm$std.error[3]),2),",",
                  round(exp(tbs.mle.norm$par[3]+qnorm(0.975,0,1)*tbs.mle.norm$std.error[3]),2),")\n")


###########
## Part of the Bayesian estimation
cat("\n"); cat("Bayesian Estimation:\n")

## Bayesian Estimation
## We run for each of the five error functions
tbs.bayes.norm     <- tbs.survreg.be(Surv(alloyT7987$time,alloyT7987$delta) ~ 1,dist="norm",
                                     guess.beta=tbs.mle.norm$par[3],
                                     guess.lambda=tbs.mle.norm$par[1],
                                     guess.xi=tbs.mle.norm$par[2],
                                     burn=500000,jump=2000,size=1000,scale=0.07)
tbs.bayes.t        <- tbs.survreg.be(Surv(alloyT7987$time,alloyT7987$delta) ~ 1,dist="t",
                                     guess.beta=tbs.mle.t$par[3],
                                     guess.lambda=tbs.mle.t$par[1],
                                     guess.xi=tbs.mle.t$par[2],
                                     burn=500000,jump=2000,size=1000,scale=0.07)
tbs.bayes.cauchy   <- tbs.survreg.be(Surv(alloyT7987$time,alloyT7987$delta) ~ 1,dist="cauchy",
                                     guess.beta=tbs.mle.cauchy$par[3],
                                     guess.lambda=tbs.mle.cauchy$par[1],
                                     guess.xi=tbs.mle.cauchy$par[2],
                                     burn=500000,jump=2000,size=1000,scale=0.07)
tbs.bayes.doubexp  <- tbs.survreg.be(Surv(alloyT7987$time,alloyT7987$delta) ~ 1,dist="doubexp",
                                     guess.beta=tbs.mle.doubexp$par[3],
                                     guess.lambda=tbs.mle.doubexp$par[1],
                                     guess.xi=tbs.mle.doubexp$par[2],
                                     burn=500000,jump=2000,size=1000,scale=0.07)
tbs.bayes.logistic <- tbs.survreg.be(Surv(alloyT7987$time,alloyT7987$delta) ~ 1,dist="logistic",
                                     guess.beta=tbs.mle.logistic$par[3],
                                     guess.lambda=tbs.mle.logistic$par[1],
                                     guess.xi=tbs.mle.logistic$par[2],
                                     burn=500000,jump=2000,size=1000,scale=0.07)

## Here it is the code for Building the Table that goes into the paper
## As before, all the info is already available in the R variables. The
## table is just for presentation purposes.
text <- paste("\\begin{center}\n",
              "\\begin{table}\n",
              "\\caption{Some quantities of the TBS Model for the Alloy T7987 data set using",
                        "BE: DIC, parameter estimates and their standard deviations in parenthesis.}\n",
              "\\label{table_alloy-bayes}\n",
              "\\centering\n",
              "\\begin{tabular}{c|c|ccc}\n",
              "\\hline\n",
              "Error Distribution & DIC & $\\widehat{\\lambda}$ & $\\widehat{\\beta_0}$ & $\\widehat{\\xi}$ \\\\ \n",
              "\\hline \n",
              "\\hline \n",
              "Normal        & ",round(tbs.bayes.norm$DIC,2),       " & ",
                                 round(tbs.bayes.norm$par[1],4),    " & ",round(tbs.bayes.norm$par[3],4),     " & ",round(tbs.bayes.norm$par[2],4),     " \\\\ \n",
              "& & {\\footnotesize(",round(tbs.bayes.norm$par.sd[1],4),")} & {\\footnotesize(",round(tbs.bayes.norm$par.sd[3],4),")} & {\\footnotesize(",round(tbs.bayes.norm$par.sd[2],5),")} \\\\ \n",
              "DoubExp       & ",round(tbs.bayes.doubexp$DIC,2),    " & ",
                                 round(tbs.bayes.doubexp$par[1],4), " & ",round(tbs.bayes.doubexp$par[3],4),  " & ",round(tbs.bayes.doubexp$par[2],4),  " \\\\ \n",
              "& & {\\footnotesize(",round(tbs.bayes.doubexp$par.sd[1],4),")} & {\\footnotesize(",round(tbs.bayes.doubexp$par.sd[3],4),")} & {\\footnotesize(",round(tbs.bayes.doubexp$par.sd[2],5),")} \\\\ \n",
              "t-Student     & ",round(tbs.bayes.t$DIC,2),          " & ",
                                 round(tbs.bayes.t$par[1],4),       " & ",round(tbs.bayes.t$par[3],4),        " & ",round(tbs.bayes.t$par[2],4),        " \\\\ \n",
              "& & {\\footnotesize(",round(tbs.bayes.t$par.sd[1],4),")} & {\\footnotesize(",round(tbs.bayes.t$par.sd[3],4),")} & {\\footnotesize(",round(tbs.bayes.t$par.sd[2],5),")} \\\\ \n",
              "Cauchy        & ",round(tbs.bayes.cauchy$DIC,2),     " & ",
                                 round(tbs.bayes.cauchy$par[1],4),  " & ",round(tbs.bayes.cauchy$par[3],4),   " & ",round(tbs.bayes.cauchy$par[2],4),   " \\\\ \n",
              "& & {\\footnotesize(",round(tbs.bayes.cauchy$par.sd[1],4),")} & {\\footnotesize(",round(tbs.bayes.cauchy$par.sd[3],4),")} & {\\footnotesize(",round(tbs.bayes.cauchy$par.sd[2],5),")} \\\\ \n",
              "Logistic      & ",round(tbs.bayes.logistic$DIC,2),   " & ",
                                 round(tbs.bayes.logistic$par[1],4)," & ",round(tbs.bayes.logistic$par[3],4), " & ",round(tbs.bayes.logistic$par[2],4), " \\\\ \n",
              "& & {\\footnotesize(",round(tbs.bayes.logistic$par.sd[1],4),")} & {\\footnotesize(",round(tbs.bayes.logistic$par.sd[3],4),")} & {\\footnotesize(",round(tbs.bayes.logistic$par.sd[2],5),")} \\\\ \n",
              "\\hline\n",
              "\\end{tabular}\n",
              "\\end{table}\n",
              "\\end{center}\n",
              sep="")
cat(text, file="table_alloy-bayes.tex")
rm(text)

## Evaluating the Survival and Hazard functions
## This code uses the output of the Bayesian estimation with the logistic error
## to build the survival and hazard functions. The logistic is chosen because it
## is the one with the best results for the particular data set, so we use it
## to illustrate the curves in the next plots.
aux.hazard   <- matrix(0,length(alloyT7987$time),length(tbs.bayes.logistic$post[,1]))
aux.survival <- matrix(0,length(alloyT7987$time),length(tbs.bayes.logistic$post[,1]))
for (j in 1:length(tbs.bayes.logistic$post[,1])) {
  aux.hazard[,j]   <- c(htbs(alloyT7987$time, lambda=tbs.bayes.logistic$post[j,1], xi=tbs.bayes.logistic$post[j,2],
                             beta=tbs.bayes.logistic$post[j,3:length(tbs.bayes.logistic$post[1,])],
                             dist="logistic"))
  aux.survival[,j] <- c(1-ptbs(alloyT7987$time, lambda=tbs.bayes.logistic$post[j,1], xi=tbs.bayes.logistic$post[j,2],
                               beta=tbs.bayes.logistic$post[j,3:length(tbs.bayes.logistic$post[1,])],
                               dist="logistic"))
}
be.hazard   <- matrix(0,length(alloyT7987$time),3)
be.survival <- matrix(0,length(alloyT7987$time),3)
for (i in 1:length(alloyT7987$time)) {
  be.hazard[i,]   <- c(  mean(aux.hazard[i,]),HPDinterval(as.mcmc(  aux.hazard[i,]),0.95))
  be.survival[i,] <- c(mean(aux.survival[i,]),HPDinterval(as.mcmc(aux.survival[i,]),0.95))
}
rm(aux.hazard,aux.survival)


i.fig <- i.fig+1
## Figure 3(a) - Reliability functions / Boundary (Bayes)
## This code generates the plot of the Kaplan-Meier and TBS BE estimation with
## the logistic error distribution. It also plots the 95% HPD credible interval
name <- paste("tbs_fig_",i.fig,"a",".eps",sep="")
eps(name,width=5,height=5,paper="special",colormodel="gray")
plot(km,ylab="",xlab="",xlim=c(min(alloyT7987$time),max(alloyT7987$time)),conf.int=FALSE,axes=FALSE,lty=1,lwd=1)
t <- seq(min(alloyT7987$time),max(alloyT7987$time),(max(alloyT7987$time)-min(alloyT7987$time)-0.01)/1000)
title(ylab="R(t)",xlab="t: number of cycles (in thousands)",main="Reliability function (BE)",cex.lab=1.2)
axis(1,at=c(93,100,150,200,250,300),lwd=2,lwd.ticks=2,pos=0)
axis(2,lwd=2,lwd.ticks=2,pos=min(alloyT7987$time))
legend(170,0.95,c("Kaplan-Meier",
                 expression(textstyle(paste("TBS / ",sep="")) ~ epsilon ~ textstyle(paste("~",sep="")) ~ Logistic),
                 "95% HPD Interval"),
                 col=c(1,"gray20","gray20"),lty=c(1,1,2),cex=1.1,lwd=c(1,2,2),bg="white")
lines(alloyT7987$time,be.survival[,1],type="l",lwd=2,col="gray20",lty=1)
lines(alloyT7987$time,be.survival[,2],type="l",lwd=2,col="gray20",lty=2)
lines(alloyT7987$time,be.survival[,3],type="l",lwd=2,col="gray20",lty=2)
dev.off()

## Figure 3(b) - Hazard functions (Bayes)
## Again, the hazard function plus the 95% HPD interval for the TBS BE with logistic error
name <- paste("tbs_fig_",i.fig,"b",".eps",sep="")
eps(name,width=5,height=5,paper="special",colormodel="gray")
 plot(alloyT7987$time,be.hazard[,1],type="l",col="gray20",lwd=2,axes="FALSE",ylim=c(0,0.05),xlab="",ylab="")
title(ylab="h(t)",xlab="t: number of cycles (in thousands)",main="Hazard function (BE)",cex.lab=1.2)
legend(100,0.048,c(expression(textstyle(paste("TBS / ",sep="")) ~ epsilon ~ textstyle(paste("~",sep="")) ~ Logistic),
                   "95% HPD Interval"),
                   col=c("gray20","gray20"),lty=c(1,2),cex=1.1,lwd=2,bg="white")
axis(1,at=c(80,100,150,200,250,300),lwd=2,lwd.ticks=2,pos=0)
axis(2,lwd=2,lwd.ticks=2,pos=min(alloyT7987$time))
lines(alloyT7987$time,be.hazard[,2],type="l",lwd=2,col="gray20",lty=2)
lines(alloyT7987$time,be.hazard[,3],type="l",lwd=2,col="gray20",lty=2)
dev.off()


#########################
## Quantile Estimation ##

## Quantile Estimation
cat("\n"); cat("Quatile estimates:\n")
cat("Median time: ",round(median(exp(tbs.bayes.logistic$post[,3])),2),"\n")
hpd <- HPDinterval(as.mcmc(exp(tbs.bayes.logistic$post[,3])),0.95)
cat("95% HPD c.i.: (",round(hpd[1],2),",",round(hpd[2],2),")\n")

## Convergence analysis of the chain
if (flag.convergence) {
  # To perform convergence analysis set "TRUE" in the if above.
  # See the first lines of this code.
  
  i.conv <- 1
  ## Figure - ACF and Time series graphics
  ## Although this generates graphs regarding the convergence analysis,
  ## they do not appear in the body of the paper.
  name <- paste("tbs_fig-conv_",i.conv,".eps",sep="")
  eps(name,width=5,height=5,paper="special",colormodel="gray")
  par(mfrow=c(2,3))
  plot(ts(tbs.bayes.logistic$post[,1]),xlab="iteration",ylab=expression(lambda),main="",type="l")
  plot(ts(tbs.bayes.logistic$post[,2]),xlab="iteration",ylab=expression(beta[0]),main="",type="l")
  plot(ts(tbs.bayes.logistic$post[,3]),xlab="iteration",ylab=expression(xi),main="",type="l")
  acf(tbs.bayes.logistic$post[,1],main=expression(lambda),ci.col="gray40")
  acf(tbs.bayes.logistic$post[,2],main=expression(beta[0]),ci.col="gray40")
  acf(tbs.bayes.logistic$post[,3],main=expression(xi),ci.col="gray40")
  par(mfrow=c(1,1))
  dev.off()
  
  # Analysis with 4 chains
  chain.logistic1 <- tbs.survreg.be(Surv(alloyT7987$time,alloyT7987$delta) ~ 1,dist="logistic",
                                    guess.beta=tbs.mle.logistic$par[3],
                                    guess.lambda=tbs.mle.logistic$par[1],
                                    guess.xi=tbs.mle.logistic$par[2],
                                    burn=1,jump=1,size=250000,scale=0.06)
  chain.logistic2 <- tbs.survreg.be(Surv(alloyT7987$time,alloyT7987$delta) ~ 1,dist="logistic",
                                    guess.beta=6,
                                    guess.lambda=1,
                                    guess.xi=1,
                                    burn=1,jump=1,size=250000,scale=0.06)
  chain.logistic3 <- tbs.survreg.be(Surv(alloyT7987$time,alloyT7987$delta) ~ 1,dist="logistic",
                                    guess.beta=4,
                                    guess.lambda=2,
                                    guess.xi=2,
                                    burn=1,jump=1,size=250000,scale=0.06)
  chain.logistic4 <- tbs.survreg.be(Surv(alloyT7987$time,alloyT7987$delta) ~ 1,dist="logistic",
                                    guess.beta=5.5,
                                    guess.lambda=0.5,
                                    guess.xi=0.5,
                                    burn=1,jump=1,size=250000,scale=0.06)
  
  n <- length(chain.logistic1$post[,1])
  aux <- seq(1,n,500) 
  ## Figure: Ergodic Means - lambda
  ## Again, the graphs about convergence are generated, but do not appear in the paper.
  i.conv <- i.conv+1
  name <- paste("tbs_fig-conv_",i.conv,".eps",sep="")
  eps(name,width=5,height=5,paper="special",colormodel="gray")
  aux.1 <- cumsum(chain.logistic1$post[,1])/seq(1,n,1)
  aux.2 <- cumsum(chain.logistic2$post[,1])/seq(1,n,1)
  aux.3 <- cumsum(chain.logistic3$post[,1])/seq(1,n,1)
  aux.4 <- cumsum(chain.logistic4$post[,1])/seq(1,n,1)
  plot(aux,aux.1[aux],type="l",xlab="iteration",ylab=expression(lambda),main="",
       ylim=c(min(aux.1,aux.2,aux.3,aux.4),max(aux.1,aux.2,aux.3,aux.4)))
  lines(aux,aux.2[aux],type="l",col="gray30")
  lines(aux,aux.3[aux],type="l",col="gray20")
  lines(aux,aux.4[aux],type="l",col="gray10")
  dev.off()

  #Figure: Ergodic Means - xi
  i.conv <- i.conv+1
  name <- paste("tbs_fig-conv_",i.conv,".eps",sep="")
  eps(name,width=5,height=5,paper="special",colormodel="gray")
  aux.1 <- cumsum(chain.logistic1$post[,2])/seq(1,n,1)
  aux.2 <- cumsum(chain.logistic2$post[,2])/seq(1,n,1)
  aux.3 <- cumsum(chain.logistic3$post[,2])/seq(1,n,1)
  aux.4 <- cumsum(chain.logistic4$post[,2])/seq(1,n,1)
   plot(aux,aux.1[aux],type="l",xlab="iteration",ylab=expression(xi),main="",
        ylim=c(min(aux.1,aux.2,aux.3,aux.4),max(aux.1,aux.2,aux.3,aux.4)))
  lines(aux,aux.2[aux],type="l",col="gray30")
  lines(aux,aux.3[aux],type="l",col="gray20")
  lines(aux,aux.4[aux],type="l",col="gray10")
  dev.off()

  #Figure: Ergodic Means - beta
  i.conv <- i.conv+1
  name <- paste("tbs_fig-conv_",i.conv,".eps",sep="")
  eps(name,width=5,height=5,paper="special",colormodel="gray")
  aux.1 <- cumsum(chain.logistic1$post[,3])/seq(1,n,1)
  aux.2 <- cumsum(chain.logistic2$post[,3])/seq(1,n,1)
  aux.3 <- cumsum(chain.logistic3$post[,3])/seq(1,n,1)
  aux.4 <- cumsum(chain.logistic4$post[,3])/seq(1,n,1)
   plot(aux,aux.1[aux],type="l",xlab="iteration",ylab=expression(beta[0]),main="",
        ylim=c(min(aux.1,aux.2,aux.3,aux.4),max(aux.1,aux.2,aux.3,aux.4)))
  lines(aux,aux.2[aux],type="l",col="gray30")
  lines(aux,aux.3[aux],type="l",col="gray20")
  lines(aux,aux.4[aux],type="l",col="gray10")
  dev.off()
  rm(aux,aux.1,aux.2,aux.3,aux.4)
  
  ## Gelman and Rubin Statistics
  ## It is presented in the paper to mention the convergence
  ## of the parameters estimated with BE
  R <- rep(0,3)
  for (j in 1:3) {
    W <- (var(chain.logistic1$post[,j])+var(chain.logistic2$post[,j])+var(chain.logistic3$post[,j])+var(chain.logistic4$post[,j]))/4
    mean.theta1 <- mean(chain.logistic1$post[,j])
    mean.theta2 <- mean(chain.logistic2$post[,j])
    mean.theta3 <- mean(chain.logistic3$post[,j])
    mean.theta4 <- mean(chain.logistic4$post[,j])
    mean.theta <- mean(c(chain.logistic1$post[,j],chain.logistic2$post[,j],chain.logistic3$post[,j],chain.logistic4$post[,j]))
    B <- (n/3)*((mean.theta1-mean.theta)^2+(mean.theta2-mean.theta)^2+(mean.theta3-mean.theta)^2+(mean.theta4-mean.theta)^2)
    R[j] <- (((1-1/n)*W+B/n)/W)
  }

  cat("\n"); cat("Gelman-Rubin Statistics:\n")
  cat("lambda: ",round(R[1],3),"\n")
  cat("    xi: ",round(R[2],3),"\n")
  cat("  beta: ",round(R[3],3),"\n")
}

############################################################################################
# Data                                                                                     #
# Colon Cancer - library(survival)                                                         #
############################################################################################
cat("\n"); cat("\n"); cat("Data - Colon\n")
data(colon)

###########
# Part MLE
cat("\n"); cat("Maximum Likelihood Estimation:\n")

# Estimating TBS model with normal error distribution
fit.mle <- tbs.survreg.mle(Surv(colon$time,colon$status) ~ colon$node4,dist="norm")

# Kaplan-Meier estimates
km <- survfit(formula = Surv(colon$time,colon$status) ~ colon$node4)

# evaluating the estimated survival functions for both groups
axis.x <- seq(0.01,3000,1)
x0 <- 1-ptbs(axis.x,lambda=fit.mle$par[1],xi=fit.mle$par[2],beta=sum(fit.mle$par[3]),dist="norm")
x1 <- 1-ptbs(axis.x,lambda=fit.mle$par[1],xi=fit.mle$par[2],beta=sum(fit.mle$par[3:4]),dist="norm")

## Figure 4(a): estimated survival functions (MLE)
## This plots the Kaplan-Meier and the TBS curves for the colon data set,
## with the patients split into 2 groups based on the variable node4
i.fig <- i.fig+1
name <- paste("tbs_fig_",i.fig,".eps",sep="")
eps(name,width=5,height=5,paper="special",colormodel="gray")
plot(km,xlim=c(0,3000),conf.int=FALSE,axes=FALSE,lty=1,lwd=1,mark.time=FALSE,xlab="",ylab="")
title(ylab="S(t)",xlab="time",main="Survival function (MLE)",cex.lab=1.2)
lines(axis.x,x0,type="l",lwd=3,col="gray50",lty=2)
lines(axis.x,x1,type="l",lwd=3,col="gray50",lty=1)
axis(1,lwd=2,lwd.ticks=2,pos=0)
axis(2,lwd=2,lwd.ticks=2,pos=0)
legend(100,0.22,c("TBS | X=0","TBS | X=1","Kaplan-Meier"),lty=c(2,1,1),lwd=c(3,3,1),
       col=c("gray30","gray30",1),cex=0.7)
dev.off()

## Quantile Estimation
## This results of the quantile estimation with MLE are not displayed in the paper,
## because they are very similar to the BE results, which we present.
cat("\n"); cat("Quatile estimates:\n")
cat("Median time | X=0: ",round(exp(fit.mle$par[3]),2),"\n")
cat("95% c.i.: (",round(exp(fit.mle$par[3]+qnorm(0.025,0,1)*fit.mle$std.error[3]),2),",",
                  round(exp(fit.mle$par[3]+qnorm(0.975,0,1)*fit.mle$std.error[3]),2),")\n")
cat("\n")
cat("Median time | X=1: ",round(exp(sum(fit.mle$par[3:4])),2),"\n")
cat("95% c.i.: (",round(exp(sum(fit.mle$par[3:4])+qnorm(0.025,0,1)*sqrt(sum(fit.mle$std.error[3:4]^2))),2),",",
                  round(exp(sum(fit.mle$par[3:4])+qnorm(0.975,0,1)*sqrt(sum(fit.mle$std.error[3:4]^2))),2),")\n")
cat("\n")
cat("Odds Median: ",round(exp(fit.mle$par[4]),2),"\n")
cat("95% c.i.: (",round(exp(fit.mle$par[4]+qnorm(0.025,0,1)*fit.mle$std.error[4]),2),",",
                  round(exp(fit.mle$par[4]+qnorm(0.975,0,1)*fit.mle$std.error[4]),2),")\n")


###########
# Part BE
cat("\n"); cat("Bayesian Estimation:\n")

## Estimating TBS model with normal error distribution
## The same as before for the colon data set with node4 as discriminative covariate,
## but now using the Bayesian estimation
fit.be <- tbs.survreg.be(Surv(colon$time,colon$status) ~ colon$node4,dist="norm",
                         guess.lambda=fit.mle$par[1],
                         guess.xi=fit.mle$par[2],
                         guess.beta=fit.mle$par[3:4],
                         burn=50000,
                         jump=500,
                         size=1000,
                         scale=0.05)

## Estimating the survival function
## Computing the survival function for presenting in figure 4(b)
axis.x <- seq(0.01,3000,1)
aux.x0 <- matrix(NA,length(axis.x),length(fit.be$post[,1]))
aux.x1 <- matrix(NA,length(axis.x),length(fit.be$post[,1]))
for (j in 1:length(fit.be$post[,1])) {
  aux.x0[,j] <- 1-ptbs(axis.x,lambda=fit.be$post[j,1],xi=fit.be$post[j,2],beta=fit.be$post[j,3],dist="norm")
  aux.x1[,j] <- 1-ptbs(axis.x,lambda=fit.be$post[j,1],xi=fit.be$post[j,1],beta=sum(fit.be$post[j,3:4]),dist="norm")
}
survival.x0 <- matrix(NA,length(axis.x),3)
survival.x1 <- matrix(NA,length(axis.x),3)
for (i in 1:length(axis.x)) {
  survival.x0[i,] <- c(mean(aux.x0[i,]),HPDinterval(as.mcmc(aux.x0[i,]),0.95))
  survival.x1[i,] <- c(mean(aux.x1[i,]),HPDinterval(as.mcmc(aux.x1[i,]),0.95))
}
rm(aux.x0,aux.x1,i,j)

## Figure 4(b): estimated survival functions (BE)
## It also shows the Kaplan-Meier and the 95% HPD credible intervals
i.fig <- i.fig+1
name <- paste("tbs_fig_",i.fig,".eps",sep="")
eps(name,width=5,height=5,paper="special",colormodel="gray")
plot(km,xlim=c(0,3000),conf.int=FALSE,axes=FALSE,lty=1,lwd=1,mark.time=FALSE,xlab="",ylab="")
title(ylab="S(t)",xlab="time",main="Survival function (BE)",cex.lab=1.2)
lines(axis.x,survival.x0[,1],type="l",lwd=3,col="gray30",lty=2)
lines(axis.x,survival.x0[,2],type="l",lwd=3,col="gray70",lty=3)
lines(axis.x,survival.x0[,3],type="l",lwd=3,col="gray70",lty=3)
lines(axis.x,survival.x1[,1],type="l",lwd=3,col="gray30",lty=1)
lines(axis.x,survival.x1[,2],type="l",lwd=3,col="gray70",lty=3)
lines(axis.x,survival.x1[,3],type="l",lwd=3,col="gray70",lty=3)
axis(1,lwd=2,lwd.ticks=2,pos=0)
axis(2,lwd=2,lwd.ticks=2,pos=0)
legend(100,0.25,c("TBS | X=0","TBS | X=1","HPD CI 95%","Kaplan-Meier"),lty=c(2,1,3,1),lwd=c(3,3,3,1),
       col=c("gray30","gray30","gray70",1),cex=0.7)
dev.off()

## Quantile Estimation
## This results are discussed in the final part of section 4.1 of the paper
cat("\n"); cat("Quatile estimates:\n")
median.0 <- mean(exp(fit.be$post[,3]))
hpd <- HPDinterval(as.mcmc(exp(fit.be$post[,3])),0.95)
cat("Median time | X=0: ",round(median.0,2),"\n")
cat("95% HPD c.i.: (",round(hpd[1],2),",",round(hpd[2],2),")\n")
median.1 <- mean(exp(apply(fit.be$post[,3:4],1,sum)))
hpd <- HPDinterval(as.mcmc(exp(apply(fit.be$post[,3:4],1,sum))),0.95)
cat("Median time | X=1: ",round(median.1,2),"\n")
cat("95% HPD c.i.: (",round(hpd[1],2),",",round(hpd[2],2),")\n")
O <- mean(exp(fit.be$post[,4]))
hpd <- HPDinterval(as.mcmc(exp(fit.be$post[,4])),0.95)
cat("Odds Median: ",round(O,2),"\n")
cat("95% HPD c.i.: (",round(hpd[1],2),",",round(hpd[2],2),")\n")
rm(median.0,median.1,O)

##########################
## To run the simulation study change to TRUE below
## Note that the simulation is a time consuming procedure, may take a few days to run.
## We have used a supercomputing centre, and still took quite some time, because there
## are many many runs...
if (FALSE) {
  source("article_simulation.r")
}

#####################
# End

## finally, just print the total time to run it.
## Without running the article_simulation.r (which many take a very long time),
## this file should run in about 2 to 3 hours in a modern computer (as of June-2012).
run.time <- (proc.time()[3]-initial.time[3])/60
print(run.time)
