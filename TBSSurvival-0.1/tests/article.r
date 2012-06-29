#library("TBSSurvival")

# To perform convergence analysis of BE set flag.convergence "TRUE".
flag.convergence <- FALSE

initial.time <- proc.time()
# Seting initial seed to have the same results that we obtained.
set.seed(1234)

##########################################
## Section 2. Transform-both-side model ##
##########################################


i.fig <- 1 # counter figure number
# Figure (a)
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

# Figure (b)
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
# Number of Cycles (in Thousands) of Fatigue Life for 67 of 72 Alloy T7987 Speciments that #
# Failed Before 300 Thousand Cycles.                                                       #
# Meeker & Escobar, pp. 130-131 (1998)                                                     #
############################################################################################
time <- c( 94,  96,  99,  99, 104, 108, 112, 114, 117, 117, 118, 121, 121, 123, 129, 131, 133, 135, 136,
          139, 139, 140, 141, 141, 143, 144, 149, 149, 152, 153, 159, 159, 159, 159, 162, 168, 168, 169,
          170, 170, 171, 172, 173, 176, 177, 180, 180, 184, 187, 188, 189, 190, 196, 197, 203, 205, 211,
          213, 224, 226, 227, 256, 257, 269, 271, 274, 291, 300, 300, 300, 300, 300)
delta <- c( 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,
            1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,
            1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,
            1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   0,   0,   0,   0,   0)

#data(alloyT7987)
#time  <- alloyT7987$time
#delta <- alloyT7987$delta

#Summary statistic of the failure time in thousand cycles for the Alloy T7987 data set.
cat("Time AlloyT7987: Summary statistics\n")
print(summary(time))


###########
# Part MLE
cat("\n"); cat("Maximum Likelihood Estimation:\n")

# MLE Estimation for each possible error distribution
tbs.mle.norm     <- tbs.survreg.mle(Surv(time,delta) ~ 1,dist="norm")
tbs.mle.doubexp  <- tbs.survreg.mle(Surv(time,delta) ~ 1,dist="doubexp")
tbs.mle.t        <- tbs.survreg.mle(Surv(time,delta) ~ 1,dist="t")
tbs.mle.cauchy   <- tbs.survreg.mle(Surv(time,delta) ~ 1,dist="cauchy")
tbs.mle.logistic <- tbs.survreg.mle(Surv(time,delta) ~ 1,dist="logistic")

# Building the Table of Example (MLE)
text <- paste("\\begin{center}\n",
              "\\begin{table}\n",
              "\\caption{Some quantities of the TBS Model for the Alloy T7987 data set ",
                        "using MLE: AIC, BIC, parameter estimates and their standard ",
                        "deviations in parenthesis.}\n",
              "\\label{table_alloy-mle}\n",
              "\\centering\n",
              "\\begin{tabular}{c|cc|ccc}\n",
              "\\hline\n",
              "Error Distribution & AIC & BIC & $\\widehat{\\lambda$} & $\\widehat{\\beta_0}$ & $\\widehat{\\xi}$ \\\\ \n",
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

# Estimating the Kaplan-Meier
km <- survfit(formula = Surv(time, delta == 1) ~ 1)

i.fig <- i.fig+1
# Figure (a) - Reliability functions (MLE)
name <- paste("tbs_fig_",i.fig,"a",".eps",sep="")
eps(name,width=5,height=5,paper="special",colormodel="gray")
plot(km,ylab="",xlab="",xlim=c(min(time),max(time)),conf.int=FALSE,axes=FALSE,lty=1,lwd=1)
t <- seq(min(time),max(time),(max(time)-min(time)-0.01)/1000)
title(ylab="R(t)",xlab="t: number of cycles (in thousands)",main="Reliability function (MLE)",cex.lab=1.2)
axis(1,at=c(93,100,150,200,250,300),lwd=2,lwd.ticks=2,pos=0)
axis(2,lwd=2,lwd.ticks=2,pos=min(time))
legend(170,0.95,c("Kaplan-Meier",
                 expression(textstyle(paste("TBS / ",sep="")) ~ epsilon ~ textstyle(paste("~",sep="")) ~ Normal)),
                 col=c(1,"gray20"),lty=c(1,1),cex=1.1,lwd=c(1,2),bg="white")
lines(t,1-ptbs(t,lambda=tbs.mle.norm$par[1],xi=tbs.mle.norm$par[2],
               beta=tbs.mle.norm$par[3],dist=tbs.mle.norm$error.dist),type="l",lwd=2,col="gray20",lty=1)
dev.off()

# Figure (b) - Hazard functions (MLE)
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

# Residual summary statistics of TBS model with normal error distribution (MLE).
cat("\n"); cat("Error: Summary statistics\n")
print(summary(tbs.mle.norm$error))


## Quantile Estimation
cat("\n"); cat("Quatile estimates:\n")
cat("Median time: ",round(exp(tbs.mle.norm$par[3]),2),"\n")
cat("95% c.i.: (",round(exp(tbs.mle.norm$par[3]+qnorm(0.025,0,1)*tbs.mle.norm$std.error[3]),2),",",
                  round(exp(tbs.mle.norm$par[3]+qnorm(0.975,0,1)*tbs.mle.norm$std.error[3]),2),")\n")


###########
# Part BE
cat("\n"); cat("Bayesian Estimation:\n")

#Bayesian Estimation
tbs.bayes.norm     <- tbs.survreg.be(Surv(time,delta) ~ 1,dist="norm",
                                     guess.beta=tbs.mle.norm$par[3],
                                     guess.lambda=tbs.mle.norm$par[1],
                                     guess.xi=tbs.mle.norm$par[2],
                                     burn=500000,jump=2000,size=1000,scale=0.07)
tbs.bayes.t        <- tbs.survreg.be(Surv(time,delta) ~ 1,dist="t",
                                     guess.beta=tbs.mle.t$par[3],
                                     guess.lambda=tbs.mle.t$par[1],
                                     guess.xi=tbs.mle.t$par[2],
                                     burn=500000,jump=2000,size=1000,scale=0.07)
tbs.bayes.cauchy   <- tbs.survreg.be(Surv(time,delta) ~ 1,dist="cauchy",
                                     guess.beta=tbs.mle.cauchy$par[3],
                                     guess.lambda=tbs.mle.cauchy$par[1],
                                     guess.xi=tbs.mle.cauchy$par[2],
                                     burn=500000,jump=2000,size=1000,scale=0.07)
tbs.bayes.doubexp  <- tbs.survreg.be(Surv(time,delta) ~ 1,dist="doubexp",
                                     guess.beta=tbs.mle.doubexp$par[3],
                                     guess.lambda=tbs.mle.doubexp$par[1],
                                     guess.xi=tbs.mle.doubexp$par[2],
                                     burn=500000,jump=2000,size=1000,scale=0.07)
tbs.bayes.logistic <- tbs.survreg.be(Surv(time,delta) ~ 1,dist="logistic",
                                     guess.beta=tbs.mle.logistic$par[3],
                                     guess.lambda=tbs.mle.logistic$par[1],
                                     guess.xi=tbs.mle.logistic$par[2],
                                     burn=500000,jump=2000,size=1000,scale=0.07)

# Building the Table of Example (Bayes)
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
              "& & & {\\footnotesize(",round(tbs.bayes.norm$par.sd[1],4),")} & {\\footnotesize(",round(tbs.bayes.norm$par.sd[3],4),")} & {\\footnotesize(",round(tbs.bayes.norm$par.sd[2],5),")} \\\\ \n",
              "DoubExp       & ",round(tbs.bayes.doubexp$DIC,2),    " & ",
                                 round(tbs.bayes.doubexp$par[1],4), " & ",round(tbs.bayes.doubexp$par[3],4),  " & ",round(tbs.bayes.doubexp$par[2],4),  " \\\\ \n",
              "& & & {\\footnotesize(",round(tbs.bayes.doubexp$par.sd[1],4),")} & {\\footnotesize(",round(tbs.bayes.doubexp$par.sd[3],4),")} & {\\footnotesize(",round(tbs.bayes.doubexp$par.sd[2],5),")} \\\\ \n",
              "t-Student     & ",round(tbs.bayes.t$DIC,2),          " & ",
                                 round(tbs.bayes.t$par[1],4),       " & ",round(tbs.bayes.t$par[3],4),        " & ",round(tbs.bayes.t$par[2],4),        " \\\\ \n",
              "& & & {\\footnotesize(",round(tbs.bayes.t$par.sd[1],4),")} & {\\footnotesize(",round(tbs.bayes.t$par.sd[3],4),")} & {\\footnotesize(",round(tbs.bayes.t$par.sd[2],5),")} \\\\ \n",
              "Cauchy        & ",round(tbs.bayes.cauchy$DIC,2),     " & ",
                                 round(tbs.bayes.cauchy$par[1],4),  " & ",round(tbs.bayes.cauchy$par[3],4),   " & ",round(tbs.bayes.cauchy$par[2],4),   " \\\\ \n",
              "& & & {\\footnotesize(",round(tbs.bayes.cauchy$par.sd[1],4),")} & {\\footnotesize(",round(tbs.bayes.cauchy$par.sd[3],4),")} & {\\footnotesize(",round(tbs.bayes.cauchy$par.sd[2],5),")} \\\\ \n",
              "Logistic      & ",round(tbs.bayes.logistic$DIC,2),   " & ",
                                 round(tbs.bayes.logistic$par[1],4)," & ",round(tbs.bayes.logistic$par[3],4), " & ",round(tbs.bayes.logistic$par[2],4), " \\\\ \n",
              "& & & {\\footnotesize(",round(tbs.bayes.logistic$par.sd[1],4),")} & {\\footnotesize(",round(tbs.bayes.logistic$par.sd[3],4),")} & {\\footnotesize(",round(tbs.bayes.logistic$par.sd[2],5),")} \\\\ \n",
              "\\hline\n",
              "\\end{tabular}\n",
              "\\end{table}\n",
              "\\end{center}\n",
              sep="")
cat(text, file="table_alloy-bayes.tex")
rm(text)


# Evaluating the Survival and Hazard functions
aux.hazard   <- matrix(0,length(time),length(tbs.bayes.logistic$post[,1]))
aux.survival <- matrix(0,length(time),length(tbs.bayes.logistic$post[,1]))
for (j in 1:length(tbs.bayes.logistic$post[,1])) {
  aux.hazard[,j]   <- c(htbs(time, lambda=tbs.bayes.logistic$post[j,1], xi=tbs.bayes.logistic$post[j,2],
                             beta=tbs.bayes.logistic$post[j,3:length(tbs.bayes.logistic$post[1,])],
                             dist="logistic"))
  aux.survival[,j] <- c(1-ptbs(time, lambda=tbs.bayes.logistic$post[j,1], xi=tbs.bayes.logistic$post[j,2],
                               beta=tbs.bayes.logistic$post[j,3:length(tbs.bayes.logistic$post[1,])],
                               dist="logistic"))
}
be.hazard   <- matrix(0,length(time),3)
be.survival <- matrix(0,length(time),3)
for (i in 1:length(time)) {
  be.hazard[i,]   <- c(  mean(aux.hazard[i,]),HPDinterval(as.mcmc(  aux.hazard[i,]),0.95))
  be.survival[i,] <- c(mean(aux.survival[i,]),HPDinterval(as.mcmc(aux.survival[i,]),0.95))
}
rm(aux.hazard,aux.survival)


i.fig <- i.fig+1
# Figure (a) - Reliability functions / Boundary (Bayes)
name <- paste("tbs_fig_",i.fig,"a",".eps",sep="")
eps(name,width=5,height=5,paper="special",colormodel="gray")
plot(km,ylab="",xlab="",xlim=c(min(time),max(time)),conf.int=FALSE,axes=FALSE,lty=1,lwd=1)
t <- seq(min(time),max(time),(max(time)-min(time)-0.01)/1000)
title(ylab="R(t)",xlab="t: number of cycles (in thousands)",main="Reliability function (BE)",cex.lab=1.2)
axis(1,at=c(93,100,150,200,250,300),lwd=2,lwd.ticks=2,pos=0)
axis(2,lwd=2,lwd.ticks=2,pos=min(time))
legend(170,0.95,c("Kaplan-Meier",
                 expression(textstyle(paste("TBS / ",sep="")) ~ epsilon ~ textstyle(paste("~",sep="")) ~ Logistic),
                 "95% HPD Interval"),
                 col=c(1,"gray20","gray20"),lty=c(1,1,2),cex=1.1,lwd=c(1,2,2),bg="white")
lines(time,be.survival[,1],type="l",lwd=2,col="gray20",lty=1)
lines(time,be.survival[,2],type="l",lwd=2,col="gray20",lty=2)
lines(time,be.survival[,3],type="l",lwd=2,col="gray20",lty=2)
dev.off()

# Figure (b) - Hazard functions (Bayes)
name <- paste("tbs_fig_",i.fig,"b",".eps",sep="")
eps(name,width=5,height=5,paper="special",colormodel="gray")
 plot(time,be.hazard[,1],type="l",col="gray20",lwd=2,axes="FALSE",ylim=c(0,0.05),xlab="",ylab="")
title(ylab="h(t)",xlab="t: number of cycles (in thousands)",main="Hazard function (BE)",cex.lab=1.2)
legend(100,0.048,c(expression(textstyle(paste("TBS / ",sep="")) ~ epsilon ~ textstyle(paste("~",sep="")) ~ Logistic),
                   "95% HPD Interval"),
                   col=c("gray20","gray20"),lty=c(1,2),cex=1.1,lwd=2,bg="white")
axis(1,at=c(80,100,150,200,250,300),lwd=2,lwd.ticks=2,pos=0)
axis(2,lwd=2,lwd.ticks=2,pos=min(time))
lines(time,be.hazard[,2],type="l",lwd=2,col="gray20",lty=2)
lines(time,be.hazard[,3],type="l",lwd=2,col="gray20",lty=2)
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
  # Figure - ACF and Time series graphics
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
  chain.logistic1 <- tbs.survreg.be(Surv(time,delta) ~ 1,dist="logistic",
                                    guess.beta=tbs.mle.logistic$par[3],
                                    guess.lambda=tbs.mle.logistic$par[1],
                                    guess.xi=tbs.mle.logistic$par[2],
                                    burn=1,jump=1,size=250000,scale=0.06)
  chain.logistic2 <- tbs.survreg.be(Surv(time,delta) ~ 1,dist="logistic",
                                    guess.beta=6,
                                    guess.lambda=1,
                                    guess.xi=1,
                                    burn=1,jump=1,size=250000,scale=0.06)
  chain.logistic3 <- tbs.survreg.be(Surv(time,delta) ~ 1,dist="logistic",
                                    guess.beta=4,
                                    guess.lambda=2,
                                    guess.xi=2,
                                    burn=1,jump=1,size=250000,scale=0.06)
  chain.logistic4 <- tbs.survreg.be(Surv(time,delta) ~ 1,dist="logistic",
                                    guess.beta=5.5,
                                    guess.lambda=0.5,
                                    guess.xi=0.5,
                                    burn=1,jump=1,size=250000,scale=0.06)
  
  n <- length(chain.logistic1$post[,1])
  aux <- seq(1,n,500) 
  #Figure: Ergodic Means - lambda
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
  name <- paste("tbs_fig-conv_",i.fig,".eps",sep="")
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
  
  # Gelman and Rubin Statistics
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


if (FALSE) {

############################################################################################
# Data                                                                                     #
# Colon Cancer - library(survival)                                                         #
############################################################################################

cat("\n"); cat("\n"); cat("Data - Colon\n")

data <- colon

###########
# Part MLE
cat("\n"); cat("Maximum Likelihood Estimation:\n")

# Estimating TBS model with normal error distribution
fit.mle <- tbs.survreg.mle(Surv(data$time,data$status) ~ data$node4,dist="norm",method="BFGS")

# Kaplan-Meier estimates
km <- survfit(formula = Surv(data$time, data$status) ~ data$node4)

# evaluating the estimated survival functions for both groups
axis.x <- seq(0.01,3000,1)
x0 <- 1-ptbs(axis.x,lambda=fit.mle$par[1],xi=fit.mle$par[2],beta=sum(fit.mle$par[3]),dist="norm")
x1 <- 1-ptbs(axis.x,lambda=fit.mle$par[1],xi=fit.mle$par[2],beta=sum(fit.mle$par[3:4]),dist="norm")

#Figure: estimated survival functions (MLE)
i.fig <- i.fig+1
name <- paste("tbs_fig-conv_",i.fig,".eps",sep="")

name <- paste("tbs_fig_ex-surv_mle-survival",".eps",sep="")
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

exp(fit.mle$par[3])
c(exp(fit.mle$par[3]+qnorm(0.025,0,1)*fit.mle$std.error[3]),
  exp(fit.mle$par[3]+qnorm(0.975,0,1)*fit.mle$std.error[3]))
exp(sum(fit.mle$par[3:4]))
c(exp(sum(fit.mle$par[3:4])+qnorm(0.025,0,1)*sqrt(sum(fit.mle$std.error[3:4]^2))),
  exp(sum(fit.mle$par[3:4])+qnorm(0.975,0,1)*sqrt(sum(fit.mle$std.error[3:4]^2))))

exp(fit.mle$par[4])
c(exp(fit.mle$par[4]+qnorm(0.025,0,1)*fit.mle$std.error[4]),
  exp(fit.mle$par[4]+qnorm(0.9755,0,1)*fit.mle$std.error[4]))

###########
# Part BE
cat("\n"); cat("Bayesian Estimation:\n")

# Estimating TBS model with normal error distribution
fit.be <- tbs.survreg.be(Surv(data$time,data$status) ~ data$node4,dist="norm",
                         guess.lambda=fit.mle$par[1],
                         guess.xi=fit.mle$par[2],
                         guess.beta=fit.mle$par[3:4],
                         burn=50000,
                         jump=500,
                         size=1000,
                         scale=0.05)


par(mfrow=c(2,4))
plot(ts(fit.be$post[,1]),xlab="iteration",ylab=expression(lambda),main="",type="l")
plot(ts(fit.be$post[,2]),xlab="iteration",ylab=expression(xi),main="",type="l")
plot(ts(fit.be$post[,3]),xlab="iteration",ylab=expression(beta[0]),main="",type="l")
plot(ts(fit.be$post[,4]),xlab="iteration",ylab=expression(beta[1]),main="",type="l")
acf(fit.be$post[,1],main=expression(lambda),ci.col="gray40")
acf(fit.be$post[,2],main=expression(xi),ci.col="gray40")
acf(fit.be$post[,3],main=expression(beta[0]),ci.col="gray40")
acf(fit.be$post[,4],main=expression(beta[1]),ci.col="gray40")
par(mfrow=c(1,1))
dev.off()


# Figures
axis.x <- seq(0.01,3000,1)
aux.x0 <- matrix(NA,length(axis.x),length(fit.be$post[,1]))
aux.x1 <- matrix(NA,length(axis.x),length(fit.be$post[,1]))
aux.x0.h <- matrix(NA,length(axis.x),length(fit.be$post[,1]))
aux.x1.h <- matrix(NA,length(axis.x),length(fit.be$post[,1]))
for (j in 1:length(fit.be$post[,1])) {
  aux.x0[,j] <- 1-ptbs(axis.x,lambda=fit.be$post[j,1],xi=fit.be$post[j,2],beta=fit.be$post[j,3],dist="norm")
  aux.x1[,j] <- 1-ptbs(axis.x,lambda=fit.be$post[j,1],xi=fit.be$post[j,1],beta=sum(fit.be$post[j,3:4]),dist="norm")
  aux.x0.h[,j] <- htbs(axis.x,lambda=fit.be$post[j,1],xi=fit.be$post[j,2],beta=fit.be$post[j,3],dist="norm")
  aux.x1.h[,j] <- htbs(axis.x,lambda=fit.be$post[j,1],xi=fit.be$post[j,1],beta=sum(fit.be$post[j,3:4]),dist="norm")
}
x0 <- matrix(NA,length(axis.x),6)
x1 <- matrix(NA,length(axis.x),6)
for (i in 1:length(axis.x)) {
  x0[i,] <- c(mean(aux.x0[i,]),HPDinterval(as.mcmc(aux.x0[i,]),0.95),
              mean(aux.x0.h[i,]),HPDinterval(as.mcmc(aux.x0.h[i,]),0.95))
  x1[i,] <- c(mean(aux.x1[i,]),HPDinterval(as.mcmc(aux.x1[i,]),0.95),
              mean(aux.x1.h[i,]),HPDinterval(as.mcmc(aux.x1.h[i,]),0.95))
}
rm(aux.x0,aux.x1,aux.x0.h,aux.x1.h,i,j)

# Example - Figure 3: Reliability functions / Boundary (Bayes)
name <- paste("tbs_fig_ex-surv_be-survival",".eps",sep="")
eps(name,width=5,height=5,paper="special",colormodel="gray")
plot(km,xlim=c(0,3000),conf.int=FALSE,axes=FALSE,lty=1,lwd=1,mark.time=FALSE,xlab="",ylab="")
title(ylab="S(t)",xlab="time",main="Survival function (BE)",cex.lab=1.2)
lines(axis.x,x0[,1],type="l",lwd=3,col="gray30",lty=2)
lines(axis.x,x0[,2],type="l",lwd=3,col="gray70",lty=3)
lines(axis.x,x0[,3],type="l",lwd=3,col="gray70",lty=3)
lines(axis.x,x1[,1],type="l",lwd=3,col="gray30",lty=1)
lines(axis.x,x1[,2],type="l",lwd=3,col="gray70",lty=3)
lines(axis.x,x1[,3],type="l",lwd=3,col="gray70",lty=3)
axis(1,lwd=2,lwd.ticks=2,pos=0)
axis(2,lwd=2,lwd.ticks=2,pos=0)
legend(100,0.25,c("TBS | X=0","TBS | X=1","HPD CI 95%","Kaplan-Meier"),lty=c(2,1,3,1),lwd=c(3,3,3,1),
       col=c("gray30","gray30","gray70",1),cex=0.7)
dev.off()

# Example - Figure 2: Hazard functions (Bayes)
name <- paste("tbs_fig_ex-surv_be-hazard",".eps",sep="")
eps(name,width=5,height=5,paper="special",colormodel="gray")
plot(axis.x,x0[,4],type="l",lwd=3,col="gray30",lty=2,xlim=c(0,3000),axes=FALSE,
     xlab="",ylab="",ylim=c(0,0.002))
title(ylab="h(t)",xlab="time",main="Hazard function (BE)",cex.lab=1.2)
lines(axis.x,x0[,5],type="l",lwd=3,col="gray70",lty=3)
lines(axis.x,x0[,6],type="l",lwd=3,col="gray70",lty=3)
lines(axis.x,x1[,4],type="l",lwd=3,col="gray30",lty=1)
lines(axis.x,x1[,5],type="l",lwd=3,col="gray70",lty=3)
lines(axis.x,x1[,6],type="l",lwd=3,col="gray70",lty=3)
axis(1,lwd=2,lwd.ticks=2,pos=0)
axis(2,lwd=2,lwd.ticks=2,pos=0)
legend(1900,0.0019,c("TBS | X=0","TBS | X=1","HPD CI 95%"),lty=c(2,1,3),lwd=c(3,3,3),
       col=c("gray30","gray30","gray70"),cex=0.7)
dev.off()


fit.be$par
fit.be$par.std.error

median.0 <- exp(fit.be$post[,3])
median.1 <- exp(apply(fit.be$post[,3:4],1,sum))
O        <- exp(fit.be$post[,4])
c(mean(median.0),HPDinterval(as.mcmc(median.0),0.95))
c(mean(median.1),HPDinterval(as.mcmc(median.1),0.95))
c(mean(O),HPDinterval(as.mcmc(O),0.95))


}







run.time <- (proc.time()[3]-initial.time[3])/60
print(run.time)
