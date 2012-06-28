library("TBSSurvival")

# To perform convergence analysis of BE set flag.convergence "TRUE".
flag.convergence <- FALSE
# To perform Gompertz simulation set flag.gompertz "TRUE".
flag.gompertz <- FALSE

initial.time <- proc.time()

# Seting initial seed to have the same results that we obtained.
set.seed(1234)

#####################################################
## Section 2.1 - Density and Reliability Functions ##
#####################################################

# Figure 1a
name <- paste("gpt_fig_1a",".eps",sep="")
t <- seq(0.01,10,0.01)
eps(name,width=5,height=5,paper="special",colormodel="gray")
plot(t, (1-pgpt(t,l=1,param=sqrt(0.5),beta=1,dist="norm")),type="l",lwd=2,lty=1,col="gray50",
     axes="FALSE",ylim=c(0,1),xlab="",ylab="")
lines(t,(1-pgpt(t,l=1,  param=sqrt(1),beta=1,dist="norm")),type="l",lwd=2,lty=2,col="gray40")
lines(t,(1-pgpt(t,l=1,  param=sqrt(2),beta=1,dist="norm")),type="l",lwd=2,lty=4,col="gray30")
lines(t,(1-pgpt(t,l=1,  param=sqrt(4),beta=1,dist="norm")),type="l",lwd=2,lty=5,col="gray20")
lines(t,(1-pgpt(t,l=1, param=sqrt(10),beta=1,dist="norm")),type="l",lwd=2,lty=6,col="gray10")
title(xlab="t",ylab="R(t)",main=expression(lambda == 1 ~ textstyle(paste("and",sep="")) ~ epsilon ~ textstyle(paste("~",sep="")) ~ Normal(0,sigma^2)),cex.lab=1.2)
axis(1,at=seq(0,10,1),lwd=2,lwd.ticks=2,pos=0)
axis(2,at=seq(0,1,0.2),lwd=2,lwd.ticks=2,pos=0)
legend(6,0.95,c(expression(sigma^2 == 0.5),
                expression(sigma^2 == 1),
                expression(sigma^2 == 2),
                expression(sigma^2 == 4),
                expression(sigma^2 == 10)),
              col=c("gray50","gray40","gray30","gray20","gray10"),lty=c(1,2,4,5,6),cex=1.1,lwd=2,bg="white")
dev.off()

# Figure 1b
name <- paste("gpt_fig_1b",".eps",sep="")
t <- seq(0.01,10,0.01)
eps(name,width=5,height=5,paper="special",colormodel="gray")
plot(t, dgpt(t,l=1,param=sqrt(0.5),beta=1,dist="norm"),type="l",lwd=2,lty=1,col="gray50",
     axes="FALSE",ylim=c(0,1),xlab="",ylab="")
lines(t,dgpt(t,l=1,  param=sqrt(1),beta=1,dist="norm"),type="l",lwd=2,lty=2,col="gray40")
lines(t,dgpt(t,l=1,  param=sqrt(2),beta=1,dist="norm"),type="l",lwd=2,lty=4,col="gray30")
lines(t,dgpt(t,l=1,  param=sqrt(4),beta=1,dist="norm"),type="l",lwd=2,lty=5,col="gray20")
lines(t,dgpt(t,l=1, param=sqrt(10),beta=1,dist="norm"),type="l",lwd=2,lty=6,col="gray10")
title(xlab="t",ylab="f(t)",main=expression(lambda == 1 ~ textstyle(paste("and",sep="")) ~ epsilon ~ textstyle(paste("~",sep="")) ~ Normal(0,sigma^2)),cex.lab=1.2)
axis(1,at=seq(0,10,1),lwd=2,lwd.ticks=2,pos=0)
axis(2,at=seq(0,1,0.2),lwd=2,lwd.ticks=2,pos=0)
legend(6,0.95,c(expression(sigma^2 == 0.5),
                expression(sigma^2 == 1),
                expression(sigma^2 == 2),
                expression(sigma^2 == 4),
                expression(sigma^2 == 10)),
              col=c("gray10","gray20","gray30","gray40","gray50"),lty=c(1,2,4,5,6),cex=1.1,lwd=2,bg="white")
dev.off()

# Figure 2a
name <- paste("gpt_fig_2a",".eps",sep="")
t <- seq(0.01,10,0.01)
eps(name,width=5,height=5,paper="special",colormodel="gray")
plot(t, (1-pgpt(t,l=2,param=sqrt(0.5),beta=1,dist="norm")),type="l",lwd=2,lty=1,col="gray50",
     axes="FALSE",ylim=c(0,1),xlab="",ylab="")
lines(t,(1-pgpt(t,l=2,  param=sqrt(1),beta=1,dist="norm")),type="l",lwd=2,lty=2,col="gray40")
lines(t,(1-pgpt(t,l=2,  param=sqrt(2),beta=1,dist="norm")),type="l",lwd=2,lty=4,col="gray30")
lines(t,(1-pgpt(t,l=2  ,param=sqrt(4),beta=1,dist="norm")),type="l",lwd=2,lty=5,col="gray20")
lines(t,(1-pgpt(t,l=2, param=sqrt(10),beta=1,dist="norm")),type="l",lwd=2,lty=6,col="gray10")
title(xlab="t",ylab="R(t)",main=expression(lambda == 2 ~ textstyle(paste("and",sep="")) ~ epsilon ~ textstyle(paste("~",sep="")) ~ Normal(0,sigma^2)),cex.lab=1.2)
axis(1,at=seq(0,10,1),lwd=2,lwd.ticks=2,pos=0)
axis(2,at=seq(0,1,0.2),lwd=2,lwd.ticks=2,pos=0)
legend(6,0.95,c(expression(sigma^2 == 0.5),
                expression(sigma^2 == 1),
                expression(sigma^2 == 2),
                expression(sigma^2 == 4),
                expression(sigma^2 == 10)),
              col=c("gray10","gray20","gray30","gray40","gray50"),lty=c(1,2,4,5,6),cex=1.1,lwd=2,bg="white")
dev.off()

# Figure 2b
name <- paste("gpt_fig_2b",".eps",sep="")
t <- seq(0.01,10,0.01)
eps(name,width=5,height=5,paper="special",colormodel="gray")
plot(t, dgpt(t,l=2,param=sqrt(0.5),beta=1,dist="norm"),type="l",lwd=2,lty=1,col="gray50",
     axes="FALSE",ylim=c(0,1),xlab="",ylab="")
lines(t,dgpt(t,l=2,  param=sqrt(1),beta=1,dist="norm"),type="l",lwd=2,lty=2,col="gray40")
lines(t,dgpt(t,l=2,  param=sqrt(2),beta=1,dist="norm"),type="l",lwd=2,lty=4,col="gray30")
lines(t,dgpt(t,l=2,  param=sqrt(4),beta=1,dist="norm"),type="l",lwd=2,lty=5,col="gray20")
lines(t,dgpt(t,l=2, param=sqrt(10),beta=1,dist="norm"),type="l",lwd=2,lty=6,col="gray10")
title(xlab="t",ylab="f(t)",main=expression(lambda == 2 ~ textstyle(paste("and",sep="")) ~ epsilon ~ textstyle(paste("~",sep="")) ~ Normal(0,sigma^2)),cex.lab=1.2)
axis(1,at=seq(0,10,1),lwd=2,lwd.ticks=2,pos=0)
axis(2,at=seq(0,1,0.2),lwd=2,lwd.ticks=2,pos=0)
legend(6,0.95,c(expression(sigma^2 == 0.5),
                expression(sigma^2 == 1),
                expression(sigma^2 == 2),
                expression(sigma^2 == 4),
                expression(sigma^2 == 10)),
              col=c("gray10","gray20","gray30","gray40","gray50"),lty=c(1,2,4,5,6),cex=1.1,lwd=2,bg="white")
dev.off()

####################################
## Section 2.2 - Hazard functions ##
####################################

# Figure 3a
name <- paste("gpt_fig_3a",".eps",sep="")
t <- seq(0.01,10,0.01)
eps(name,width=5,height=5,paper="special",colormodel="gray")
plot(t, hazard.gpt(t,l=0.5,param=sqrt(0.5),beta=1,dist="norm"),type="l",lwd=2,lty=1,col="gray50",
     axes="FALSE",ylim=c(0,1),xlab="",ylab="")
lines(t,hazard.gpt(t,l=0.5,  param=sqrt(1),beta=1,dist="norm"),type="l",lwd=2,lty=2,col="gray40")
lines(t,hazard.gpt(t,l=0.5,  param=sqrt(2),beta=1,dist="norm"),type="l",lwd=2,lty=4,col="gray30")
lines(t,hazard.gpt(t,l=0.5,  param=sqrt(4),beta=1,dist="norm"),type="l",lwd=2,lty=5,col="gray20")
lines(t,hazard.gpt(t,l=0.5, param=sqrt(10),beta=1,dist="norm"),type="l",lwd=2,lty=6,col="gray10")
title(xlab="t",ylab="h(t)",main=expression(lambda == 0.5 ~ textstyle(paste("and",sep="")) ~ epsilon ~ textstyle(paste("~",sep="")) ~ Normal(0,sigma^2)),cex.lab=1.2)
axis(1,at=seq(0,10,1),lwd=2,lwd.ticks=2,pos=0)
axis(2,at=seq(0,1,0.2),lwd=2,lwd.ticks=2,pos=0)
legend(6,0.95,c(expression(sigma^2 == 0.5),
                expression(sigma^2 == 1),
                expression(sigma^2 == 2),
                expression(sigma^2 == 4),
                expression(sigma^2 == 10)),
              col=c("gray10","gray20","gray30","gray40","gray50"),lty=c(1,2,4,5,6),cex=1.1,lwd=2,bg="white")
dev.off()

# Figure 3b
name <- paste("gpt_fig_3b",".eps",sep="")
t <- seq(0.01,10,0.01)
eps(name,width=5,height=5,paper="special",colormodel="gray")
plot(t, hazard.gpt(t,l=1,param=sqrt(0.5),beta=1,dist="norm"),type="l",lwd=2,lty=1,col="gray50",
     axes="FALSE",ylim=c(0,1),xlab="",ylab="")
lines(t,hazard.gpt(t,l=1,  param=sqrt(1),beta=1,dist="norm"),type="l",lwd=2,lty=2,col="gray40")
lines(t,hazard.gpt(t,l=1,  param=sqrt(2),beta=1,dist="norm"),type="l",lwd=2,lty=4,col="gray30")
lines(t,hazard.gpt(t,l=1,  param=sqrt(4),beta=1,dist="norm"),type="l",lwd=2,lty=5,col="gray20")
lines(t,hazard.gpt(t,l=1, param=sqrt(10),beta=1,dist="norm"),type="l",lwd=2,lty=6,col="gray10")
title(xlab="t",ylab="h(t)",main=expression(lambda == 1 ~ textstyle(paste("and",sep="")) ~ epsilon ~ textstyle(paste("~",sep="")) ~ Normal(0,sigma^2)),cex.lab=1.2)
axis(1,at=seq(0,10,1),lwd=2,lwd.ticks=2,pos=0)
axis(2,at=seq(0,1,0.2),lwd=2,lwd.ticks=2,pos=0)
legend(6,0.95,c(expression(sigma^2 == 0.5),
                expression(sigma^2 == 1),
                expression(sigma^2 == 2),
                expression(sigma^2 == 4),
                expression(sigma^2 == 10)),
              col=c("gray10","gray20","gray30","gray40","gray50"),lty=c(1,2,4,5,6),cex=1.1,lwd=2,bg="white")
dev.off()

# Figure 3c
name <- paste("gpt_fig_3c",".eps",sep="")
t <- seq(0.01,10,0.01)
eps(name,width=5,height=5,paper="special",colormodel="gray")
plot(t, hazard.gpt(t,l=2,param=sqrt(0.5),beta=1,dist="norm"),type="l",lwd=2,lty=1,col="gray50",
     axes="FALSE",ylim=c(0,1),xlab="",ylab="")
lines(t,hazard.gpt(t,l=2,  param=sqrt(1),beta=1,dist="norm"),type="l",lwd=2,lty=2,col="gray40")
lines(t,hazard.gpt(t,l=2,  param=sqrt(2),beta=1,dist="norm"),type="l",lwd=2,lty=4,col="gray30")
lines(t,hazard.gpt(t,l=2,  param=sqrt(4),beta=1,dist="norm"),type="l",lwd=2,lty=5,col="gray20")
lines(t,hazard.gpt(t,l=2, param=sqrt(10),beta=1,dist="norm"),type="l",lwd=2,lty=6,col="gray10")
title(xlab="t",ylab="h(t)",main=expression(lambda == 2 ~ textstyle(paste("and",sep="")) ~ epsilon ~ textstyle(paste("~",sep="")) ~ Normal(0,sigma^2)),cex.lab=1.2)
axis(1,at=seq(0,10,1),lwd=2,lwd.ticks=2,pos=0)
axis(2,at=seq(0,1,0.2),lwd=2,lwd.ticks=2,pos=0)
legend(6,0.95,c(expression(sigma^2 == 0.5),
                expression(sigma^2 == 1),
                expression(sigma^2 == 2),
                expression(sigma^2 == 4),
                expression(sigma^2 == 10)),
              col=c("gray10","gray20","gray30","gray40","gray50"),lty=c(1,2,4,5,6),cex=1.1,lwd=2,bg="white")
dev.off()

# Figure 3d
name <- paste("gpt_fig_3d",".eps",sep="")
t <- seq(0.01,10,0.01)
eps(name,width=5,height=5,paper="special",colormodel="gray")
plot(t, hazard.gpt(t,l=3,param=sqrt(0.5),beta=1,dist="norm"),type="l",lwd=2,lty=1,col="gray50",
     axes="FALSE",ylim=c(0,1),xlab="",ylab="")
lines(t,hazard.gpt(t,l=3,  param=sqrt(1),beta=1,dist="norm"),type="l",lwd=2,lty=2,col="gray40")
lines(t,hazard.gpt(t,l=3,  param=sqrt(2),beta=1,dist="norm"),type="l",lwd=2,lty=4,col="gray30")
lines(t,hazard.gpt(t,l=3,  param=sqrt(4),beta=1,dist="norm"),type="l",lwd=2,lty=5,col="gray20")
lines(t,hazard.gpt(t,l=3, param=sqrt(10),beta=1,dist="norm"),type="l",lwd=2,lty=6,col="gray10")
title(xlab="t",ylab="h(t)",main=expression(lambda == 3 ~ textstyle(paste("and",sep="")) ~ epsilon ~ textstyle(paste("~",sep="")) ~ Normal(0,sigma^2)),cex.lab=1.2)
axis(1,at=seq(0,10,1),lwd=2,lwd.ticks=2,pos=0)
axis(2,at=seq(0,1,0.2),lwd=2,lwd.ticks=2,pos=0)
legend(6,0.95,c(expression(sigma^2 == 0.5),
                expression(sigma^2 == 1),
                expression(sigma^2 == 2),
                expression(sigma^2 == 4),
                expression(sigma^2 == 10)),
              col=c("gray10","gray20","gray30","gray40","gray50"),lty=c(1,2,4,5,6),cex=1.1,lwd=2,bg="white")
dev.off()

# Figure 4a
name <- paste("gpt_fig_4a",".eps",sep="")
t <- seq(0.01,10,0.01)
eps(name,width=5,height=5,paper="special",colormodel="gray")
plot(t, hazard.gpt(t,l=0.5,param=sqrt(0.5),beta=1,dist="doubexp"),type="l",lwd=2,lty=1,col="gray50",
     axes="FALSE",ylim=c(0,1),xlab="",ylab="")
lines(t,hazard.gpt(t,l=0.5,  param=sqrt(1),beta=1,dist="doubexp"),type="l",lwd=2,lty=2,col="gray40")
lines(t,hazard.gpt(t,l=0.5,  param=sqrt(2),beta=1,dist="doubexp"),type="l",lwd=2,lty=4,col="gray30")
lines(t,hazard.gpt(t,l=0.5,  param=sqrt(4),beta=1,dist="doubexp"),type="l",lwd=2,lty=5,col="gray20")
lines(t,hazard.gpt(t,l=0.5, param=sqrt(10),beta=1,dist="doubexp"),type="l",lwd=2,lty=6,col="gray10")
title(xlab="t",ylab="h(t)",main=expression(lambda == 0.5 ~ textstyle(paste("and",sep="")) ~ epsilon ~ textstyle(paste("~",sep="")) ~ DoubExp(b)),cex.lab=1.2)
axis(1,at=seq(0,10,1),lwd=2,lwd.ticks=2,pos=0)
axis(2,at=seq(0,1,0.2),lwd=2,lwd.ticks=2,pos=0)
legend(6,0.95,c(expression(b == 0.5),
                expression(b == 1),
                expression(b == 2),
                expression(b == 4),
                expression(b == 10)),
              col=c("gray10","gray20","gray30","gray40","gray50"),lty=c(1,2,4,5,6),cex=1.1,lwd=2,bg="white")
dev.off()

# Figure 4b
name <- paste("gpt_fig_4b",".eps",sep="")
t <- seq(0.01,10,0.01)
eps(name,width=5,height=5,paper="special",colormodel="gray")
plot(t, hazard.gpt(t,l=1,param=sqrt(0.5),beta=1,dist="doubexp"),type="l",lwd=2,lty=1,col="gray50",
     axes="FALSE",ylim=c(0,1),xlab="",ylab="")
lines(t,hazard.gpt(t,l=1,  param=sqrt(1),beta=1,dist="doubexp"),type="l",lwd=2,lty=2,col="gray40")
lines(t,hazard.gpt(t,l=1,  param=sqrt(2),beta=1,dist="doubexp"),type="l",lwd=2,lty=4,col="gray30")
lines(t,hazard.gpt(t,l=1,  param=sqrt(4),beta=1,dist="doubexp"),type="l",lwd=2,lty=5,col="gray20")
lines(t,hazard.gpt(t,l=1, param=sqrt(10),beta=1,dist="doubexp"),type="l",lwd=2,lty=6,col="gray10")
title(xlab="t",ylab="h(t)",main=expression(lambda == 1 ~ textstyle(paste("and",sep="")) ~ epsilon ~ textstyle(paste("~",sep="")) ~ DoubExp(b)),cex.lab=1.2)
axis(1,at=seq(0,10,1),lwd=2,lwd.ticks=2,pos=0)
axis(2,at=seq(0,1,0.2),lwd=2,lwd.ticks=2,pos=0)
legend(6,0.95,c(expression(b == 0.5),
                expression(b == 1),
                expression(b == 2),
                expression(b == 4),
                expression(b == 10)),
              col=c("gray10","gray20","gray30","gray40","gray50"),lty=c(1,2,4,5,6),cex=1.1,lwd=2,bg="white")
dev.off()

# Figure 4c
name <- paste("gpt_fig_4c",".eps",sep="")
t <- seq(0.01,10,0.01)
eps(name,width=5,height=5,paper="special",colormodel="gray")
plot(t, hazard.gpt(t,l=2,param=sqrt(0.5),beta=1,dist="doubexp"),type="l",lwd=2,lty=1,col="gray50",
     axes="FALSE",ylim=c(0,1),xlab="",ylab="")
lines(t,hazard.gpt(t,l=2,  param=sqrt(1),beta=1,dist="doubexp"),type="l",lwd=2,lty=2,col="gray40")
lines(t,hazard.gpt(t,l=2,  param=sqrt(2),beta=1,dist="doubexp"),type="l",lwd=2,lty=4,col="gray30")
lines(t,hazard.gpt(t,l=2,  param=sqrt(4),beta=1,dist="doubexp"),type="l",lwd=2,lty=5,col="gray20")
lines(t,hazard.gpt(t,l=2, param=sqrt(10),beta=1,dist="doubexp"),type="l",lwd=2,lty=6,col="gray10")
title(xlab="t",ylab="h(t)",main=expression(lambda == 2 ~ textstyle(paste("and",sep="")) ~ epsilon ~ textstyle(paste("~",sep="")) ~ DoubExp(b)),cex.lab=1.2)
axis(1,at=seq(0,10,1),lwd=2,lwd.ticks=2,pos=0)
axis(2,at=seq(0,1,0.2),lwd=2,lwd.ticks=2,pos=0)
legend(6,0.95,c(expression(b == 0.5),
                expression(b == 1),
                expression(b == 2),
                expression(b == 4),
                expression(b == 10)),
              col=c("gray10","gray20","gray30","gray40","gray50"),lty=c(1,2,4,5,6),cex=1.1,lwd=2,bg="white")
dev.off()

# Figure 4d
name <- paste("gpt_fig_4d",".eps",sep="")
t <- seq(0.01,10,0.01)
eps(name,width=5,height=5,paper="special",colormodel="gray")
plot(t, hazard.gpt(t,l=3,param=sqrt(0.5),beta=1,dist="doubexp"),type="l",lwd=2,lty=1,col="gray50",
     axes="FALSE",ylim=c(0,1),xlab="",ylab="")
lines(t,hazard.gpt(t,l=3,  param=sqrt(1),beta=1,dist="doubexp"),type="l",lwd=2,lty=2,col="gray40")
lines(t,hazard.gpt(t,l=3,  param=sqrt(2),beta=1,dist="doubexp"),type="l",lwd=2,lty=4,col="gray30")
lines(t,hazard.gpt(t,l=3,  param=sqrt(4),beta=1,dist="doubexp"),type="l",lwd=2,lty=5,col="gray20")
lines(t,hazard.gpt(t,l=3, param=sqrt(10),beta=1,dist="doubexp"),type="l",lwd=2,lty=6,col="gray10")
title(xlab="t",ylab="h(t)",main=expression(lambda == 3 ~ textstyle(paste("and",sep="")) ~ epsilon ~ textstyle(paste("~",sep="")) ~ DoubExp(b)),cex.lab=1.2)
axis(1,at=seq(0,10,1),lwd=2,lwd.ticks=2,pos=0)
axis(2,at=seq(0,1,0.2),lwd=2,lwd.ticks=2,pos=0)
legend(6,0.95,c(expression(b == 0.5),
                expression(b == 1),
                expression(b == 2),
                expression(b == 4),
                expression(b == 10)),
              col=c("gray10","gray20","gray30","gray40","gray50"),lty=c(1,2,4,5,6),cex=1.1,lwd=2,bg="white")
dev.off()


##############################
## Section 4 - Data Example ##
##############################

############################################################################################
# Data                                                                                     #
# Number of Cycles (in Thousands) of Fatigue Life for 67 of 72 Alloy T7987 Speciments that #
# Failed Before 300 Thousand Cycles.                                                       #
# Meeker & Escobar, pp. 130-131 (1998)                                                     #
############################################################################################
d <- NULL
d$time <- c( 94,  96,  99,  99, 104, 108, 112, 114, 117, 117, 118, 121, 121, 123, 129, 131, 133, 135, 136,
            139, 139, 140, 141, 141, 143, 144, 149, 149, 152, 153, 159, 159, 159, 159, 162, 168, 168, 169,
            170, 170, 171, 172, 173, 176, 177, 180, 180, 184, 187, 188, 189, 190, 196, 197, 203, 205, 211,
            213, 224, 226, 227, 256, 257, 269, 271, 274, 291, 300, 300, 300, 300, 300)
d$delta <- c( 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,
              1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,
              1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,
              1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   0,   0,   0,   0,   0)

# MLE Estimation
# error: Normal
gpt.mle.norm     <- estimation.mle(time=d$time,delta=d$delta,dist="norm"    ,method="Nelder-Mead")
# error: Double-Exponential
gpt.mle.doubexp  <- estimation.mle(time=d$time,delta=d$delta,dist="doubexp" ,method="Nelder-Mead")
# error: t-Student
gpt.mle.t        <- estimation.mle(time=d$time,delta=d$delta,dist="t"       ,method="Nelder-Mead")
# error: Cauchy
gpt.mle.cauchy   <- estimation.mle(time=d$time,delta=d$delta,dist="cauchy"  ,method="Nelder-Mead")
# error: Logistic
gpt.mle.logistic <- estimation.mle(time=d$time,delta=d$delta,dist="logistic",method="Nelder-Mead")

###############################
## Benchmark - Model Weibull ##
weib.initial.time <- proc.time()
weib.model <- NULL
mle.weib             <- weibreg(Surv(d$time,d$delta)~1)
weib.model$par       <- c(exp(mle.weib$coefficients[1]),exp(mle.weib$coefficients[2]))
names(weib.model$par) <- NULL
weib.model$std.error <- diag(sqrt(solve(-hessian(lweib,weib.model$par,time=d$time,delta=d$delta))))
weib.model$log.lik   <- mle.weib$loglik[2]
weib.model$nparam    <- 2
weib.model$AIC       <- 2*weib.model$nparam-2*weib.model$log.lik
weib.model$AICc      <- weib.model$AIC + 2*weib.model$nparam*(weib.model$nparam+1)/(length(d$time)-weib.model$nparam-1)
weib.model$BIC       <- -2*weib.model$log.lik+weib.model$nparam*log(length(d$time))
km <- survfit(formula = Surv(d$time, d$delta == 1) ~ 1)
weib.model$KS        <- max(abs(exp(-(km$time/weib.model$par[1])^weib.model$par[2])-km$surv))
weib.model$run.time  <- (proc.time()[3] - weib.initial.time[3])/60
rm(weib.initial.time,mle.weib,km)

# Building the Table of Example (MLE)
text <- paste("\\begin{center}\n",
              "\\begin{table}\n",
              "\\caption{GPT Model (MLE)}\n",
              "\\label{table_mle}\n",
              "\\centering\n",
              "\\begin{tabular}{c|ccc|ccc}\n",
              "\\hline\n",
              "Error Distribution & KS & AIC & BIC & $\\lambda$ & $\\beta_0$ & $\\xi$ \\\\ \n",
              "\\hline \n",
              "\\hline \n",
              "Normal        & ",round(gpt.mle.norm$KS,4),        " & ",round(gpt.mle.norm$AIC,2),        " & ",round(gpt.mle.norm$BIC,2),        " & ",
                                 round(gpt.mle.norm$par[1],4),    " & ",round(gpt.mle.norm$par[3],4),     " & ",round(gpt.mle.norm$par[2],4),     " \\\\ \n",
              "& & & & {\\footnotesize(",round(gpt.mle.norm$std.error[1],4),")} & {\\footnotesize(",round(gpt.mle.norm$std.error[3],4),")} & {\\footnotesize(",round(gpt.mle.norm$std.error[2],5),")} \\\\ \n",
              "DoubExp       & ",round(gpt.mle.doubexp$KS,4),     " & ",round(gpt.mle.doubexp$AIC,2),     " & ",round(gpt.mle.doubexp$BIC,2),     " & ",
                                 round(gpt.mle.doubexp$par[1],4), " & ",round(gpt.mle.doubexp$par[3],4),  " & ",round(gpt.mle.doubexp$par[2],4),  " \\\\ \n",
              "& & & & {\\footnotesize(",round(gpt.mle.doubexp$std.error[1],4),")} & {\\footnotesize(",round(gpt.mle.doubexp$std.error[3],4),")} & {\\footnotesize(",round(gpt.mle.doubexp$std.error[2],5),")} \\\\ \n",
              "t-Student     & ",round(gpt.mle.t$KS,4),           " & ",round(gpt.mle.t$AIC,2),           " & ",round(gpt.mle.t$BIC,2),           " & ",
                                 round(gpt.mle.t$par[1],4),       " & ",round(gpt.mle.t$par[3],4),        " & ",round(gpt.mle.t$par[2],4),        " \\\\ \n",
              "& & & & {\\footnotesize(",round(gpt.mle.t$std.error[1],4),")} & {\\footnotesize(",round(gpt.mle.t$std.error[3],4),")} & {\\footnotesize(",round(gpt.mle.t$std.error[2],5),")} \\\\ \n",
              "Cauchy        & ",round(gpt.mle.cauchy$KS,4),      " & ",round(gpt.mle.cauchy$AIC,2),      " & ",round(gpt.mle.cauchy$BIC,2),      " & ",
                                 round(gpt.mle.cauchy$par[1],4),  " & ",round(gpt.mle.cauchy$par[3],4),   " & ",round(gpt.mle.cauchy$par[2],4),   " \\\\ \n",
              "& & & & {\\footnotesize(",round(gpt.mle.cauchy$std.error[1],4),")} & {\\footnotesize(",round(gpt.mle.cauchy$std.error[3],4),")} & {\\footnotesize(",round(gpt.mle.cauchy$std.error[2],5),")} \\\\ \n",
              "Logistic      & ",round(gpt.mle.logistic$KS,4),    " & ",round(gpt.mle.logistic$AIC,2),    " & ",round(gpt.mle.logistic$BIC,2),    " & ",
                                 round(gpt.mle.logistic$par[1],4)," & ",round(gpt.mle.logistic$par[3],4), " & ",round(gpt.mle.logistic$par[2],4), " \\\\ \n",
              "& & & & {\\footnotesize(",round(gpt.mle.logistic$std.error[1],4),")} & {\\footnotesize(",round(gpt.mle.logistic$std.error[3],4),")} & {\\footnotesize(",round(gpt.mle.logistic$std.error[2],5),")} \\\\ \n",
              "\\hline\n",
              "\\hline\n",
              "& & & & Scale & Shape & \\\\ \n",
              "\\hline\n",
              "Weibull Model & ",round(weib.model$KS,4),      " & ",round(weib.model$AIC,2),      " & ",round(weib.model$BIC,2),      " & ",
                                 round(weib.model$par[1],4)," & ",round(weib.model$par[2],4)," & \\\\ \n",
              "& & & & {\\footnotesize(",round(weib.model$std.error[1],4),")} & {\\footnotesize(",round(weib.model$std.error[2],4),")} &  \\\\ \n",
              "\\hline\n",
              "\\end{tabular}\n",
              "\\end{table}\n",
              "\\end{center}\n",
              sep="")
cat(text, file="table_ex-mle.tex")
rm(text)

# Example - Figure 1: Reliability functions (MLE)
name <- paste("gpt_fig_ex-mle1",".eps",sep="")
eps(name,width=5,height=5,paper="special",colormodel="gray")
km <- survfit(formula = Surv(d$time, d$delta == 1) ~ 1)
plot(km,ylab="",xlab="",xlim=c(min(d$time),max(d$time)),conf.int=FALSE,axes=FALSE,lty=1,lwd=1)
t <- seq(min(d$time),max(d$time),(max(d$time)-min(d$time)-0.01)/1000)
title(ylab="R(t)",xlab="t: number of cycles (in thousands)",main="Reliability function (MLE)",cex.lab=1.2)
axis(1,at=c(93,100,150,200,250,300),lwd=2,lwd.ticks=2,pos=0)
axis(2,lwd=2,lwd.ticks=2,pos=min(d$time))
legend(170,0.95,c("Kaplan-Meier",
                 expression(textstyle(paste("GPT / ",sep="")) ~ epsilon ~ textstyle(paste("~",sep="")) ~ Logistic),
                 "Weibull Model"),
                 col=c(1,"gray20","gray10"),lty=c(1,1,2),cex=1.1,lwd=c(1,2,2),bg="white")
lines(t,1-pgpt(t,l=gpt.mle.logistic$par[1],param=gpt.mle.logistic$par[2],
               beta=gpt.mle.logistic$par[3],dist=gpt.mle.logistic$error.dist),type="l",lwd=2,col="gray20",lty=1)
lines(t,exp(-(t/weib.model$par[1])^weib.model$par[2]) ,type="l",lwd=2,col="gray10",lty=2)
dev.off()

# Example - Figure 2: Hazard functions (MLE)
name <- paste("gpt_fig_ex-mle2",".eps",sep="")
eps(name,width=5,height=5,paper="special",colormodel="gray")
t <- seq(0.01,350,(350-0.01)/1000)
 plot(t,hazard.gpt(t,l=gpt.mle.logistic$par[1],param=gpt.mle.logistic$par[2],
                   beta=gpt.mle.logistic$par[3],dist=gpt.mle.logistic$error.dist),
      type="l",col="gray20",lwd=2,axes="FALSE",ylim=c(0,0.05),xlab="",ylab="")
title(ylab="h(t)",xlab="t: number of cycles (in thousands)",main="Hazard function (MLE)",cex.lab=1.2)
legend(30,0.048,c(expression(textstyle(paste("GPT / ",sep="")) ~ epsilon ~ textstyle(paste("~",sep="")) ~ Logistic),
                   "Weibull Model"),
                   col=c("gray20","gray10"),lty=c(1,2),cex=1.1,lwd=2,bg="white")
axis(1,lwd=2,lwd.ticks=2,pos=0)
axis(2,lwd=2,lwd.ticks=2,pos=0)
lines(t,(weib.model$par[2]/weib.model$par[1])*((t/weib.model$par[1])^(weib.model$par[2]-1)),type="l",col="gray10",lwd=2,lty=2)
lines(t,hazard.gpt(t,l=gpt.mle.logistic$par[1],param=gpt.mle.logistic$par[2],
                   beta=gpt.mle.logistic$par[3],dist=gpt.mle.logistic$error.dist),type="l",col="gray20",lwd=2)
dev.off()

# Residual analysis of GPT model with logistic error. (MLE)
name <- paste("gpt_fig_ex-mle_error1",".eps",sep="")
eps(name,width=5,height=5,paper="special",colormodel="gray")
error <- (g.l(log(km$time),gpt.mle.logistic$par[1]) - g.l(gpt.mle.logistic$par[3],gpt.mle.logistic$par[1]))
pplot(error,qdist=qlogistic,s=gpt.mle.logistic$par[2],xlab=expression(epsilon),main="Quantil Plot - Logistic (MLE)")
dev.off()

name <- paste("gpt_fig_ex-mle_error2",".eps",sep="")
eps(name,width=5,height=5,paper="special",colormodel="gray")
hist(error,freq = FALSE, main=expression(paste("Histogram of ",epsilon," (MLE)",sep="")),
                         xlab=expression(epsilon))
curve(dlogistic(x,s=gpt.mle.logistic$par[2]),col = "gray10", lty = 1, lwd = 2, add = TRUE)
dev.off()

round(summary(error),5)

# Kolmogorov-Smirnov test
ks.test(error, "plogistic", s=gpt.mle.logistic$par[2])


#########################
## Quantile Estimation ##

# Building the Table of Example
text <- paste("\\begin{center}\n",
              "\\begin{table}\n",
              "\\caption{Quantile estimation (MLE)}\n",
              "\\label{table_mle}\n",
              "\\centering\n",
              "\\begin{tabular}{c|c|cc}\n",
              "\\hline\n",
              "Error Distribution & Median & $L_B$ & $U_B$ \\\\ \n",
              "\\hline \n",
              "\\hline \n",
#              "Normal    & ",round(exp(gpt.mle.norm$par[3]),2)," & ",
#                round(exp(gpt.mle.norm$par[3]+qnorm(0.025,0,1)*gpt.mle.norm$std.error[3]),2)," & ",
#                round(exp(gpt.mle.norm$par[3]+qnorm(0.975,0,1)*gpt.mle.norm$std.error[3]),2)," \\\\ \n",
#              "DoubExp   & ",round(exp(gpt.mle.doubexp$par[3]),2)," & ",
#                round(exp(gpt.mle.doubexp$par[3]+qnorm(0.025,0,1)*gpt.mle.doubexp$std.error[3]),2)," & ",
#                round(exp(gpt.mle.doubexp$par[3]+qnorm(0.975,0,1)*gpt.mle.doubexp$std.error[3]),2)," \\\\ \n",
#              "t-Student & ",round(exp(gpt.mle.t$par[3]),2)," & ",
#                round(exp(gpt.mle.t$par[3]+qnorm(0.025,0,1)*gpt.mle.t$std.error[3]),2)," & ",
#                round(exp(gpt.mle.t$par[3]+qnorm(0.975,0,1)*gpt.mle.t$std.error[3]),2)," \\\\ \n",
#              "Cauchy    & ",round(exp(gpt.mle.cauchy$par[3]),2)," & ",
#                round(exp(gpt.mle.cauchy$par[3]+qnorm(0.025,0,1)*gpt.mle.cauchy$std.error[3]),2)," & ",
#                round(exp(gpt.mle.cauchy$par[3]+qnorm(0.975,0,1)*gpt.mle.cauchy$std.error[3]),2)," \\\\ \n",
              "Logistic  & ",round(exp(gpt.mle.logistic$par[3]),2)," & ",
                round(exp(gpt.mle.logistic$par[3]+qnorm(0.025,0,1)*gpt.mle.logistic$std.error[3]),2)," & ",
                round(exp(gpt.mle.logistic$par[3]+qnorm(0.975,0,1)*gpt.mle.logistic$std.error[3]),2)," \\\\ \n",
              "\\hline\n",
              "\\end{tabular}\n",
              "\\end{table}\n",
              "\\end{center}\n",
              sep="")
cat(text, file="table_ex-mle-quantile.tex")
rm(text)


#Bayesian Estimation
gpt.bayes.norm     <- estimation.bayes(    guess=gpt.mle.norm$par,time=d$time,delta=d$delta,dist="norm",
                                       burn=500000,jump=2000,size=1000,scale=0.07)
gpt.bayes.t        <- estimation.bayes(       guess=gpt.mle.t$par,time=d$time,delta=d$delta,dist="t",
                                       burn=500000,jump=2000,size=1000,scale=0.1)
gpt.bayes.cauchy   <- estimation.bayes(  guess=gpt.mle.cauchy$par,time=d$time,delta=d$delta,dist="cauchy",
                                       burn=500000,jump=2000,size=1000,scale=0.1)
gpt.bayes.doubexp  <- estimation.bayes( guess=gpt.mle.doubexp$par,time=d$time,delta=d$delta,dist="doubexp",
                                       burn=500000,jump=2000,size=1000,scale=0.1)
gpt.bayes.logistic <- estimation.bayes(guess=gpt.mle.logistic$par,time=d$time,delta=d$delta,dist="logistic",
                                       burn=500000,jump=2000,size=1000,scale=0.06)

#bayes.weib.model <- metrop(lweib,weib.model$par,time=d$time,delta=d$delta,nbatch=200000,scale=1)
bayes.weib.model <- bweib(guess=weib.model$par,time=d$time,delta=d$delta,
                          burn=20000,jump=180,size=1000,scale=1)

# Building the Table of Example (Bayes)
text <- paste("\\begin{center}\n",
              "\\begin{table}\n",
              "\\caption{GPT Model (BE)}\n",
              "\\label{table_bayes}\n",
              "\\centering\n",
              "\\begin{tabular}{c|cc|ccc}\n",
              "\\hline\n",
              "Error Distribution & KS & DIC & $\\lambda$ & $\\beta_0$ & $\\xi$ \\\\ \n",
              "\\hline \n",
              "\\hline \n",
              "Normal        & ",round(gpt.bayes.norm$KS,4),        " & ",round(gpt.bayes.norm$DIC,2),        " & ",
                                 round(gpt.bayes.norm$par[1],4),    " & ",round(gpt.bayes.norm$par[3],4),     " & ",round(gpt.bayes.norm$par[2],4),     " \\\\ \n",
              "& & & {\\footnotesize(",round(gpt.bayes.norm$par.std.error[1],4),")} & {\\footnotesize(",round(gpt.bayes.norm$par.std.error[3],4),")} & {\\footnotesize(",round(gpt.bayes.norm$par.std.error[2],5),")} \\\\ \n",
              "DoubExp       & ",round(gpt.bayes.doubexp$KS,4),     " & ",round(gpt.bayes.doubexp$DIC,2),     " & ",
                                 round(gpt.bayes.doubexp$par[1],4), " & ",round(gpt.bayes.doubexp$par[3],4),  " & ",round(gpt.bayes.doubexp$par[2],4),  " \\\\ \n",
              "& & & {\\footnotesize(",round(gpt.bayes.doubexp$par.std.error[1],4),")} & {\\footnotesize(",round(gpt.bayes.doubexp$par.std.error[3],4),")} & {\\footnotesize(",round(gpt.bayes.doubexp$par.std.error[2],5),")} \\\\ \n",
              "t-Student     & ",round(gpt.bayes.t$KS,4),           " & ",round(gpt.bayes.t$DIC,2),           " & ",
                                 round(gpt.bayes.t$par[1],4),       " & ",round(gpt.bayes.t$par[3],4),        " & ",round(gpt.bayes.t$par[2],4),        " \\\\ \n",
              "& & & {\\footnotesize(",round(gpt.bayes.t$par.std.error[1],4),")} & {\\footnotesize(",round(gpt.bayes.t$par.std.error[3],4),")} & {\\footnotesize(",round(gpt.bayes.t$par.std.error[2],5),")} \\\\ \n",
              "Cauchy        & ",round(gpt.bayes.cauchy$KS,4),      " & ",round(gpt.bayes.cauchy$DIC,2),      " & ",
                                 round(gpt.bayes.cauchy$par[1],4),  " & ",round(gpt.bayes.cauchy$par[3],4),   " & ",round(gpt.bayes.cauchy$par[2],4),   " \\\\ \n",
              "& & & {\\footnotesize(",round(gpt.bayes.cauchy$par.std.error[1],4),")} & {\\footnotesize(",round(gpt.bayes.cauchy$par.std.error[3],4),")} & {\\footnotesize(",round(gpt.bayes.cauchy$par.std.error[2],5),")} \\\\ \n",
              "Logistic      & ",round(gpt.bayes.logistic$KS,4),    " & ",round(gpt.bayes.logistic$DIC,2),    " & ",
                                 round(gpt.bayes.logistic$par[1],4)," & ",round(gpt.bayes.logistic$par[3],4), " & ",round(gpt.bayes.logistic$par[2],4), " \\\\ \n",
              "& & & {\\footnotesize(",round(gpt.bayes.logistic$par.std.error[1],4),")} & {\\footnotesize(",round(gpt.bayes.logistic$par.std.error[3],4),")} & {\\footnotesize(",round(gpt.bayes.logistic$par.std.error[2],5),")} \\\\ \n",
              "\\hline\n",
              "\\hline\n",
              "& & & Scale & Shape & \\\\ \n",
              "\\hline\n",
              "Weibull Model & ",round(bayes.weib.model$KS,4),      " & ",round(bayes.weib.model$DIC,2),      " & ",
                                 round(bayes.weib.model$par[1],4),  " & ",round(bayes.weib.model$par[2],4)," & \\\\ \n",
              "& & & {\\footnotesize(",round(bayes.weib.model$par.std.error[1],4),")} & {\\footnotesize(",round(bayes.weib.model$par.std.error[2],4),")} &  \\\\ \n",
              "\\hline\n",
              "\\end{tabular}\n",
              "\\end{table}\n",
              "\\end{center}\n",
              sep="")
cat(text, file="table_ex-bayes.tex")
rm(text)

# Example - Figure 1: Reliability functions (Bayes)
name <- paste("gpt_fig_ex-bayes1",".eps",sep="")
eps(name,width=5,height=5,paper="special",colormodel="gray")
km <- survfit(formula = Surv(d$time, d$delta == 1) ~ 1)
plot(km,ylab="",xlab="",xlim=c(min(d$time),max(d$time)),conf.int=FALSE,axes=FALSE,lty=1,lwd=1)
t <- seq(min(d$time),max(d$time),(max(d$time)-min(d$time)-0.01)/1000)
title(ylab="R(t)",xlab="t: number of cycles (in thousands)",main="Reliability function (BE)",cex.lab=1.2)
axis(1,at=c(93,100,150,200,250,300),lwd=2,lwd.ticks=2,pos=0)
axis(2,lwd=2,lwd.ticks=2,pos=min(d$time))
legend(170,0.95,c("Kaplan-Meier",
                 expression(textstyle(paste("GPT / ",sep="")) ~ epsilon ~ textstyle(paste("~",sep="")) ~ Logistic),
                 "Weibull Model"),
                 col=c(1,"gray20","gray10"),lty=c(1,1,2),cex=1.1,lwd=c(1,2,2),bg="white")
lines(gpt.bayes.logistic$time,gpt.bayes.logistic$survival[,3],type="l",lwd=2,col="gray20",lty=1)
lines(  bayes.weib.model$time,  bayes.weib.model$survival[,3],type="l",lwd=2,col="gray10",lty=2)
dev.off()

# Example - Figure 2: Hazard functions (Bayes)
name <- paste("gpt_fig_ex-bayes2",".eps",sep="")
eps(name,width=5,height=5,paper="special",colormodel="gray")
 plot(gpt.bayes.logistic$time,gpt.bayes.logistic$hazard[,3],
      type="l",col="gray20",lwd=2,axes="FALSE",ylim=c(0,0.05),xlab="",ylab="")
title(ylab="h(t)",xlab="t: number of cycles (in thousands)",main="Hazard function (BE)",cex.lab=1.2)
legend(100,0.048,c(expression(textstyle(paste("GPT / ",sep="")) ~ epsilon ~ textstyle(paste("~",sep="")) ~ Logistic),
                   "Weibull Model"),
                   col=c("gray20","gray10"),lty=c(1,2),cex=1.1,lwd=2,bg="white")
axis(1,at=c(80,100,150,200,250,300),lwd=2,lwd.ticks=2,pos=0)
axis(2,lwd=2,lwd.ticks=2,pos=min(gpt.bayes.logistic$time))
lines(  bayes.weib.model$time,  bayes.weib.model$hazard[,3],type="l",col="gray10",lwd=2,lty=2)
lines(gpt.bayes.logistic$time,gpt.bayes.logistic$hazard[,3],type="l",col="gray20",lwd=2,lty=1)
dev.off()

# Example - Figure 3: Reliability functions / Boundary (Bayes)
name <- paste("gpt_fig_ex-bayes-boundary",".eps",sep="")
eps(name,width=5,height=5,paper="special",colormodel="gray")
km <- survfit(formula = Surv(d$time, d$delta == 1) ~ 1)
plot(km,ylab="",xlab="",xlim=c(min(d$time),max(d$time)),conf.int=FALSE,axes=FALSE,lty=1,lwd=1)
t <- seq(min(d$time),max(d$time),(max(d$time)-min(d$time)-0.01)/1000)
title(ylab="R(t)",xlab="t: number of cycles (in thousands)",main="Reliability function (BE)",cex.lab=1.2)
axis(1,at=c(93,100,150,200,250,300),lwd=2,lwd.ticks=2,pos=0)
axis(2,lwd=2,lwd.ticks=2,pos=min(d$time))
legend(170,0.95,c("Kaplan-Meier",
                 expression(textstyle(paste("GPT / ",sep="")) ~ epsilon ~ textstyle(paste("~",sep="")) ~ Logistic),
                 "95% HPD Interval"),
                 col=c(1,"gray20","gray20"),lty=c(1,1,2),cex=1.1,lwd=c(1,2,2),bg="white")
lines(gpt.bayes.logistic$time,gpt.bayes.logistic$survival[,3],type="l",lwd=2,col="gray20",lty=1)
lines(gpt.bayes.logistic$time,gpt.bayes.logistic$survival[,7],type="l",lwd=2,col="gray20",lty=2)
lines(gpt.bayes.logistic$time,gpt.bayes.logistic$survival[,8],type="l",lwd=2,col="gray20",lty=2)
dev.off()


#########################
## Quantile Estimation ##

# Building the Table of Example (Bayes)
text <- paste("\\begin{center}\n",
              "\\begin{table}\n",
              "\\caption{Quantile estimation (BE)}\n",
              "\\label{table_quantile-bayes}\n",
              "\\centering\n",
              "\\begin{tabular}{c|c|cc}\n",
              "\\hline\n",
              "Error Distribution & Median & $L_B$ & $U_B$ \\\\ \n",
              "\\hline \n",
              "\\hline \n",
#              "Normal    & ",    round(gpt.bayes.norm$t.median[3],2)," & ",    round(gpt.bayes.norm$t.median[9],2)," & ",    round(gpt.bayes.norm$t.median[10],2)," \\\\ \n",
#              "Doubexp   & ", round(gpt.bayes.doubexp$t.median[3],2)," & ", round(gpt.bayes.doubexp$t.median[9],2)," & ", round(gpt.bayes.doubexp$t.median[10],2)," \\\\ \n",
#              "t         & ",       round(gpt.bayes.t$t.median[3],2)," & ",       round(gpt.bayes.t$t.median[9],2)," & ",       round(gpt.bayes.t$t.median[10],2)," \\\\ \n",
#              "Cauchy    & ",  round(gpt.bayes.cauchy$t.median[3],2)," & ",  round(gpt.bayes.cauchy$t.median[9],2)," & ",  round(gpt.bayes.cauchy$t.median[10],2)," \\\\ \n",
              "Logistic  & ",round(gpt.bayes.logistic$t.median[3],2)," & ",round(gpt.bayes.logistic$t.median[9],2)," & ",round(gpt.bayes.logistic$t.median[10],2)," \\\\ \n",
              "\\hline\n",
              "\\end{tabular}\n",
              "\\end{table}\n",
              "\\end{center}\n",
              sep="")
cat(text, file="table_ex-bayes-quantile.tex")
rm(text)

# Figures: Prior .vs. Posterior
name <- paste("estimation_prior_x_posterior-lambda",".eps",sep="")
eps(name,width=5,height=5,paper="special",colormodel="gray")
#Density Lambda
aux.dens <- density(gpt.bayes.logistic$post[,1])
aux.x <- seq(0,5,0.01) 
aux.y <- dunif.exp(aux.x,a=0.00001,b=3,p=0.8)
plot(aux.x,aux.y,type="l",lwd=2,ylim=c(0,max(aux.y,aux.dens$y)),ylab="density",xlab=expression(lambda))
lines(aux.dens,lwd=2,col="gray30")
legend(2.5,1,c("prior","posterior"),col=c(1,"gray30"),lty=c(1,1),lwd=c(2,2),cex=1)
dev.off()
rm(aux.dens,aux.x,aux.y)
#Density xi
name <- paste("estimation_prior_x_posterior-xi",".eps",sep="")
eps(name,width=5,height=5,paper="special",colormodel="gray")
aux.dens <- density(gpt.bayes.logistic$post[,2])
aux.x <- seq(0,4,0.01) 
aux.y <- dunif.exp(aux.x,a=0.00001,b=2,p=0.9)
plot(aux.x,aux.y,type="l",lwd=2,ylim=c(0,max(aux.y,aux.dens$y)),ylab="density",xlab=expression(xi))
lines(aux.dens,lwd=2,col="gray30")
legend(2.5,0.55,c("prior","posterior"),col=c(1,"gray30"),lty=c(1,1),lwd=c(2,2),cex=1)
dev.off()
rm(aux.dens,aux.x,aux.y)
#Density beta_0
name <- paste("estimation_prior_x_posterior-beta",".eps",sep="")
eps(name,width=5,height=5,paper="special",colormodel="gray")
aux.dens <- density(gpt.bayes.logistic$post[,3])
aux.x <- seq(4.7,5.5,0.1) 
aux.y <- dnorm(aux.x,5,5)
plot(aux.x,aux.y,type="l",lwd=2,ylim=c(0,max(aux.y,aux.dens$y)),ylab="density",xlab=expression(beta[0]))
lines(aux.dens,lwd=2,col="gray30")
legend(5.2,9,c("prior","posterior"),col=c(1,"gray30"),lty=c(1,1),lwd=c(2,2),cex=1)
dev.off()
rm(aux.dens,aux.x,aux.y)
#end


## Convergence analysis of the chain
if (flag.convergence) {
  # To perform convergence analysis set "TRUE" in the if above.
  # See the first lines of this code.
  
  name <- paste("gpt_fig_ex-bayes-converge-post",".eps",sep="")
  eps(name,width=5,height=5,paper="special",colormodel="gray")
#  par(mfrow=c(3,3))
  par(mfrow=c(2,3))
  plot(ts(gpt.bayes.logistic$post[,1]),xlab="iteration",ylab=expression(lambda),main="",type="l")
  plot(ts(gpt.bayes.logistic$post[,2]),xlab="iteration",ylab=expression(beta[0]),main="",type="l")
  plot(ts(gpt.bayes.logistic$post[,3]),xlab="iteration",ylab=expression(xi),main="",type="l")
  acf(gpt.bayes.logistic$post[,1],main=expression(lambda),ci.col="gray40")
  acf(gpt.bayes.logistic$post[,2],main=expression(beta[0]),ci.col="gray40")
  acf(gpt.bayes.logistic$post[,3],main=expression(xi),ci.col="gray40")
#  plot(seq(1,length(gpt.bayes.logistic$post[,1]),1),
#       cumsum(gpt.bayes.logistic$post[,1])/seq(1,length(gpt.bayes.logistic$post[,1]),1),
#       xlab="iteration",ylab=expression(lambda),main="",type="l")
#  plot(seq(1,length(gpt.bayes.logistic$post[,2]),1),
#       cumsum(gpt.bayes.logistic$post[,2])/seq(1,length(gpt.bayes.logistic$post[,2]),1),
#       xlab="iteration",ylab=expression(beta[0]),main="",type="l")
#  plot(seq(1,length(gpt.bayes.logistic$post[,3]),1),
#       cumsum(gpt.bayes.logistic$post[,3])/seq(1,length(gpt.bayes.logistic$post[,3]),1),
#       xlab="iteration",ylab=expression(xi),main="",type="l")
  par(mfrow=c(1,1))
  dev.off()
  
  # Analysis with 4 chains
  chain.logistic1 <- metrop(log.post,gpt.mle.logistic$par,time=d$time,
                            delta=d$delta,dist="logistic",nbatch=250000,blen=1,scale=0.06)
  chain.logistic2 <- metrop(log.post,c(1,1,6),            time=d$time,
                            delta=d$delta,dist="logistic",nbatch=250000,blen=1,scale=0.06)
  chain.logistic3 <- metrop(log.post,c(2,2,4),            time=d$time,
                            delta=d$delta,dist="logistic",nbatch=250000,blen=1,scale=0.06)
  chain.logistic4 <- metrop(log.post,c(0.5,0.5,5.5),      time=d$time,
                            delta=d$delta,dist="logistic",nbatch=250000,blen=1,scale=0.06)
  
  #Figure: Ergodic Means
  n <- length(chain.logistic1$batch[,1])
  aux <- seq(1,n,500) 
  name <- paste("gpt_fig_ex-bayes-converge-chain_lambda",".eps",sep="")
  eps(name,width=5,height=5,paper="special",colormodel="gray")
  aux.1 <- cumsum(chain.logistic1$batch[,1])/seq(1,n,1)
  aux.2 <- cumsum(chain.logistic2$batch[,1])/seq(1,n,1)
  aux.3 <- cumsum(chain.logistic3$batch[,1])/seq(1,n,1)
  aux.4 <- cumsum(chain.logistic4$batch[,1])/seq(1,n,1)
   plot(aux,aux.1[aux],type="l",xlab="iteration",ylab=expression(lambda),main="",
        ylim=c(min(aux.1,aux.2,aux.3,aux.4),max(aux.1,aux.2,aux.3,aux.4)))
  lines(aux,aux.2[aux],type="l",col="gray30")
  lines(aux,aux.3[aux],type="l",col="gray20")
  lines(aux,aux.4[aux],type="l",col="gray10")
  dev.off()
  name <- paste("gpt_fig_ex-bayes-converge-chain_xi",".eps",sep="")
  eps(name,width=5,height=5,paper="special",colormodel="gray")
  aux.1 <- cumsum(chain.logistic1$batch[,2])/seq(1,n,1)
  aux.2 <- cumsum(chain.logistic2$batch[,2])/seq(1,n,1)
  aux.3 <- cumsum(chain.logistic3$batch[,2])/seq(1,n,1)
  aux.4 <- cumsum(chain.logistic4$batch[,2])/seq(1,n,1)
   plot(aux,aux.1[aux],type="l",xlab="iteration",ylab=expression(xi),main="",
        ylim=c(min(aux.1,aux.2,aux.3,aux.4),max(aux.1,aux.2,aux.3,aux.4)))
  lines(aux,aux.2[aux],type="l",col="gray30")
  lines(aux,aux.3[aux],type="l",col="gray20")
  lines(aux,aux.4[aux],type="l",col="gray10")
  dev.off()
  name <- paste("gpt_fig_ex-bayes-converge-chain_beta0",".eps",sep="")
  eps(name,width=5,height=5,paper="special",colormodel="gray")
  aux.1 <- cumsum(chain.logistic1$batch[,3])/seq(1,n,1)
  aux.2 <- cumsum(chain.logistic2$batch[,3])/seq(1,n,1)
  aux.3 <- cumsum(chain.logistic3$batch[,3])/seq(1,n,1)
  aux.4 <- cumsum(chain.logistic4$batch[,3])/seq(1,n,1)
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
    W <- (var(chain.logistic1$batch[,j])+var(chain.logistic2$batch[,j])+var(chain.logistic3$batch[,j])+var(chain.logistic4$batch[,j]))/4
    mean.theta1 <- mean(chain.logistic1$batch[,j])
    mean.theta2 <- mean(chain.logistic2$batch[,j])
    mean.theta3 <- mean(chain.logistic3$batch[,j])
    mean.theta4 <- mean(chain.logistic4$batch[,j])
    mean.theta <- mean(c(chain.logistic1$batch[,j],chain.logistic2$batch[,j],chain.logistic3$batch[,j],chain.logistic4$batch[,j]))
    B <- (n/3)*((mean.theta1-mean.theta)^2+(mean.theta2-mean.theta)^2+(mean.theta3-mean.theta)^2+(mean.theta4-mean.theta)^2)
    R[j] <- (((1-1/n)*W+B/n)/W)
  }
  text <- paste("\\begin{center}\n",
                "\\begin{table}\n",
                "\\caption{Gelman-Rubin Statistics (BE)}\n",
                "\\label{table_gelman}\n",
                "\\centering\n",
                "\\begin{tabular}{cccc}\n",
                "\\hline\n",
                " & $\\lambda$ & $\\xi$ & $\\beta_0$ \\\\ \n",
                "\\hline \n",
                "\\hline \n",
                "Logistic  & ",round(R[1],5)," & ",round(R[2],5)," & ",round(R[3],5)," \\\\ \n",
                "\\hline\n",
                "\\end{tabular}\n",
                "\\end{table}\n",
                "\\end{center}\n",
                sep="")
  cat(text, file="table_ex-bayes-gelman.tex")
  rm(text)

}

####################################
## Section 5 - Concluding Remarks ##
####################################
if (flag.gompertz) {
  # To perform simulation of Section 5 set "TRUE" in the if above.
  # See the first lines of this code.

sample.size <- c(  30, 100,30,100,30,100,30,100)
gomp.gamma  <- c(3/10,3/10, 1,  1, 1,  1, 2,  2)
gomp.k      <- c(10/3,10/3, 1,  1, 2,  2, 1,  1)
n.rep <- 100

text <- paste("\\begin{center}\n","\\begin{table}\n",
              "\\caption{KS - Comparison between Gompertz distribution and GPT model with normal error distribution (MLE)}\n",
              "\\label{table_gompertz}\n","\\centering\n","\\begin{tabular}{cc|ccccc}\n","\\hline\n",
              " n & Parameters & Min. & 1st Qu. & Median & 3rd. Qu. & Max. \\\\ \n","\\hline \n","\\hline \n",sep="")
for (j in 1:length(sample.size)) {
  KS <- rep(0,n.rep)
  sec5 <- NULL
  for (i in 1:n.rep) {
    sec5$time <- rgompertz(sample.size[j],gomp.gamma[j],1/gomp.k[j])
    sec5$delta <- rep(1,length(sec5$time))
    gpt.mle.sec5 <- estimation.mle(time=sec5$time,delta=sec5$delta,dist="norm")
    KS[i] <- gpt.mle.sec5$KS
  }
  aux <- summary(KS)
  text <- paste(text,sample.size[j]," & $\\gamma = ",round(gomp.gamma[j],2),"$, $\\kappa = ",round(gomp.k[j],2),"$ & ",
                round(aux[1],4)," & ",round(aux[2],4)," & ",round(aux[3],4)," & ",round(aux[5],4)," & ",round(aux[6],4)," \\\\ \n",sep="")
}
text <- paste(text,"\\hline\n","\\end{tabular}\n","\\end{table}\n","\\end{center}\n",sep="")
cat(text, file="table_sec5-gomp.tex")
rm(text)

sec5$time <- rgompertz(30,0.3,0.3)
sec5$delta <- rep(1,length(sec5$time))
gpt.mle.sec5 <- estimation.mle(time=sec5$time,delta=sec5$delta,dist="norm")
name <- paste("gpt_fig_sec5-gompertz-30",".eps",sep="")
eps(name,width=5,height=5,paper="special",colormodel="gray")
km <- survfit(formula = Surv(sec5$time, sec5$delta == 1) ~ 1)
plot(km,ylab="",xlab="",xlim=c(0,max(sec5$time)),conf.int=FALSE,axes=FALSE,lty=1,lwd=2)
t <- seq(min(sec5$time),max(sec5$time),(max(sec5$time)-min(sec5$time)-0.01)/1000)
title(ylab="R(t)",xlab="time (n=30)",main="Reliability function (MLE)",cex.lab=1.2)
axis(1,at=c(-0.1,seq(0,max(sec5$time)+0.19,0.2)),lwd=2,lwd.ticks=2,pos=0)
axis(2,lwd=2,lwd.ticks=2,pos=0)
legend(0.05,0.25,c("Kaplan-Meier",
                 expression(textstyle(paste("GPT / ",sep="")) ~ epsilon ~ textstyle(paste("~",sep="")) ~ Normal)),
                 col=c(1,"gray20"),lty=c(1,1),cex=1.1,lwd=c(2,2),bg="white")
lines(t,1-pgpt(t,l=gpt.mle.sec5$par[1],param=gpt.mle.sec5$par[2],
               beta=gpt.mle.sec5$par[3],dist=gpt.mle.sec5$error.dist),type="l",lwd=2,col="gray20",lty=1)
dev.off()

sec5$time <- rgompertz(100,0.3,0.3)
sec5$delta <- rep(1,length(sec5$time))
gpt.mle.sec5 <- estimation.mle(time=sec5$time,delta=sec5$delta,dist="norm")
name <- paste("gpt_fig_sec5-gompertz-100",".eps",sep="")
eps(name,width=5,height=5,paper="special",colormodel="gray")
km <- survfit(formula = Surv(sec5$time, sec5$delta == 1) ~ 1)
plot(km,ylab="",xlab="",xlim=c(0,max(sec5$time)),conf.int=FALSE,axes=FALSE,lty=1,lwd=2)
t <- seq(min(sec5$time),max(sec5$time),(max(sec5$time)-min(sec5$time)-0.01)/1000)
title(ylab="R(t)",xlab="time (n=100)",main="Reliability function (MLE)",cex.lab=1.2)
axis(1,at=c(-0.1,seq(0,max(sec5$time)+0.19,0.2)),lwd=2,lwd.ticks=2,pos=0)
axis(2,lwd=2,lwd.ticks=2,pos=0)
legend(0.05,0.25,c("Kaplan-Meier",
                 expression(textstyle(paste("GPT / ",sep="")) ~ epsilon ~ textstyle(paste("~",sep="")) ~ Normal)),
                 col=c(1,"gray20"),lty=c(1,1),cex=1.1,lwd=c(2,2),bg="white")
lines(t,1-pgpt(t,l=gpt.mle.sec5$par[1],param=gpt.mle.sec5$par[2],
               beta=gpt.mle.sec5$par[3],dist=gpt.mle.sec5$error.dist),type="l",lwd=2,col="gray20",lty=1)
dev.off()

}

run.time <- (proc.time()[3]-initial.time[3])/60
print(run.time)
