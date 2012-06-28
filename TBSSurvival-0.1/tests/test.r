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

## This code is used for testing purposes. The TBSSurvival library does not
## depend on it for any of its functionalities
source('loadlibs.r')
loadlibs()

######################
# simulation: function used to generate data from Weibull distribution and
#             then compare the Weibull model with TBSmodel.
#             Also, this function is used to analyse if TBSmodel MLE is 
#             working well.
sim.weib <- function(gen=list(n=100,cens=0,shape=2,scale=3),est=list(dist="norm"),prefix="",
                     n.copies=1,initial.seed=1234,method="Rsolnp",verbose=TRUE) {
  # gen: is a list with the parameters to generate from a Weibull distribution
  #          n: sample size
  #       cens: censor rate
  #      shape: shape Weibull parameter
  #      scale: scale Weibull parameter
  # est: is a list with the parameters used to estimate the model
  #      dist: error distribution of TBS model
  # n.copies: how many copies will be performed in the simulation
  # initial.seed: to control the simulation

  result.1  <- matrix(est$dist,n.copies,1)
  result.2  <- matrix(NA,n.copies,23)

  for (i in 1:n.copies) {
    # Generating data
    initial.time <- .gettime()
    seed <- initial.seed+i-1
    set.seed(seed)
    d <- NULL
    d$time  <- rweibull(gen$n,shape=1/gen$shape,scale=exp(gen$scale))
    d$delta <- rep(1,gen$n)
    censor  <- quantile(d$time,probs=(1-gen$cens))
    d$delta <- ifelse(d$time > censor, 0, 1)
    d$time  <- ifelse(d$time > censor, censor, d$time)

    # Original Survival function
    orig.y   <- pweibull(d$time[d$delta == 1],shape=1/gen$shape,scale=exp(gen$scale))

    # Weibull fit
    weib.fit <- survreg(Surv(d$time,d$delta) ~ 1, dist="weibull")
    weib.y   <- pweibull(d$time[d$delta == 1],shape=1/weib.fit$scale,
                                              scale=exp(weib.fit$coefficients))

    # Tbs fit
    tbs.fit  <- tbs.survreg.mle(Surv(d$time,d$delta) ~ 1,dist=est$dist,method=method,verbose=verbose)

    # Result matrix
    if (tbs.fit$convergence) {
      tbs.y <- ptbs(d$time[d$delta == 1],lambda=tbs.fit$par[1],xi=tbs.fit$par[2],
                    beta=tbs.fit$par[3],dist=est$dist)
      run.time <- .gettime() - initial.time
      result.2[i,] <- c(seed,gen$n,gen$cens,gen$shape,gen$scale,tbs.fit$par,
                        max(abs(orig.y-tbs.y)),mean(abs(orig.y-tbs.y)),
                        1/weib.fit$scale,exp(weib.fit$coefficients),
                        max(abs(orig.y-weib.y)),mean(abs(orig.y-weib.y)),
                        ifelse(max(abs(orig.y-tbs.y)) > max(abs(orig.y-weib.y)), 0, 1),
                        ifelse(mean(abs(orig.y-tbs.y)) > mean(abs(orig.y-weib.y)), 0, 1),
                        exp(gen$scale)*log(2)^gen$shape,median(d$time),
                        exp(tbs.fit$par[3]),
                        (exp(tbs.fit$par[3])-exp(gen$scale)*log(2)^gen$shape)^2,
                        exp(weib.fit$coefficients)*log(2)^weib.fit$scale,
                        (exp(weib.fit$coefficients)*log(2)^weib.fit$scale-
                         exp(gen$scale)*log(2)^gen$shape)^2,run.time)
    } else {
      run.time <- .gettime() - initial.time
      result.2[i,] <- c(seed,gen$n,gen$cens,gen$shape,gen$scale,rep(NA,5),
                        1/weib.fit$scale,exp(weib.fit$coefficients),
                        max(abs(orig.y-weib.y)),mean(abs(orig.y-weib.y)),
                        rep(NA,2),exp(gen$scale)*log(2)^gen$shape,median(d$time),NA,NA,
                        exp(weib.fit$coefficients)*log(2)^weib.fit$scale,
                        (exp(weib.fit$coefficients)*log(2)^weib.fit$scale-
                         exp(gen$scale)*log(2)^gen$shape)^2,run.time)
    }
    rm(d)
  }
  
  colnames(result.1) <- "dist"
  colnames(result.2) <- c("gen.seed","n","censor","orig.par1","orig.par2",
                          "tbs.lambda","tbs.xi","tbs.beta","tbs.Max.AE",
                          "tbs.Mean.AE","weib.par1","weib.par2",
                          "weib.Max.AE","weib.Mean.AE","comp.Max.AE","comp.Mean.AE",
                          "orig.median","data.median","tbs.median","tbs.sqe",
                          "weib.median","weib.sqe","run.time")

  out <- data.frame(result.1,result.2)
  name <- paste(prefix,"Sim_Weib-",est$dist,"_",gen$n,"_",gen$cens,"_",gen$shape,"_",gen$scale,".csv",sep="")
  write.csv(out,file=name)
#  return(out)
}

######################
# sim.tbs: function used to generate data from TBS model and 
#          then check the estimation procedure.  
sim.tbs <- function(gen=list(n=100,cens=0,lambda=2,xi=3,beta=1,dist="norm"),prefix="",
                    n.copies=1,initial.seed=1234,method="Rsolnp",verbose=TRUE) {
  # gen: is a list with the parameters to generate from a TBS model
  #          n: sample size
  #       cens: censor rate
  #     lambda, xi, beta: TBS parameters
  # n.copies: how many copies will be performed in the simulation
  # initial.seed: to control the simulation

  result.1  <- matrix(gen$dist,n.copies,1)
  result.2  <- matrix(NA,n.copies,19)

  for (i in 1:n.copies) {
    # Generating data
    initial.time <- .gettime()
    seed <- initial.seed+i-1
    set.seed(seed)
    d <- NULL
    d$time  <- rtbs(gen$n,gen$lambda,gen$xi,gen$beta,dist=gen$dist)
    d$delta <- rep(1,gen$n)
    censor  <- quantile(d$time,probs=(1-gen$cens))
    d$delta <- ifelse(d$time > censor, 0, 1)
    d$time  <- ifelse(d$time > censor, censor, d$time)

    # Original Survival function
    orig.y  <- ptbs(d$time[d$delta == 1],gen$n,gen$lambda,gen$xi,gen$beta,dist=gen$dist)

    # Tbs fit
    tbs.fit  <- tbs.survreg.mle(Surv(d$time,d$delta) ~ 1,dist=gen$dist,method=method,verbose=verbose)

    # Result matrix
    if (tbs.fit$convergence) {
      tbs.y <- ptbs(d$time[d$delta == 1],lambda=tbs.fit$par[1],xi=tbs.fit$par[2],
                    beta=tbs.fit$par[3],dist=gen$dist)
      run.time <- .gettime() - initial.time
      result.2[i,] <- c(seed,gen$n,gen$cens,gen$lambda,gen$xi,gen$beta,tbs.fit$par,
                        (gen$lambda-tbs.fit$par[1])^2,(gen$xi-tbs.fit$par[2])^2,
                        (gen$beta-tbs.fit$par[3])^2,
                        max(abs(orig.y-tbs.y)),mean(abs(orig.y-tbs.y)),median(d$time),
                        exp(gen$beta),exp(tbs.fit$par[3]),(exp(tbs.fit$par[3])-exp(gen$beta))^2,
                        run.time)
    } else {
      run.time <- .gettime() - initial.time
      result.2[i,] <- c(seed,gen$n,gen$cens,gen$lambda,gen$xi,gen$beta,rep(NA,8),
                        median(d$time),rep(NA,3),run.time)
    }
    rm(d)
  }
  
  colnames(result.1) <- "dist"
  colnames(result.2) <- c("gen.seed","n","censor","orig.lambda","orig.xi","orig.beta",
                          "tbs.lambda","tbs.xi","tbs.beta","lambda.sqe","xi.sqe","beta.sqe",
                          "tbs.Max.AE","tbs.Mean.AE","orig.median","data.median",
                          "tbs.median","tbs.sqe","run.time")
  out <- data.frame(result.1,result.2)

  name <- paste(prefix,"Sim_TBS-",gen$dist,"_",gen$n,"_",gen$cens,"_",gen$lambda,
                "_",gen$xi,"_",gen$beta,".csv",sep="")
  write.csv(out,file=name)
#  return(out)
}


if(FALSE) {


  ###############################################################################
## The following function is based on the work of Evgenia Dimitriadou,
## Kurt Hornik, Friedrich Leisch, David Meyer, and Andreas Weingessel.
## This function is a version of the "probplot" function from package "e1071".
## The changes are in the color (to gray scale) and the labels of the graphic.
pplot <- function (x, qdist = qnorm, probs = NULL, line = TRUE, xlab = NULL, 
    ylab = "Probability in %", main="", ...) 
{
    DOTARGS <- as.list(substitute(list(...)))[-1]
    DOTARGS <- paste(names(DOTARGS), DOTARGS, sep = "=", collapse = ", ")
    x <- sort(x)
    QNAME <- deparse(substitute(qdist))
    DOTS <- list(...)
    qdist <- match.fun(qdist)
    QFUN <- function(p) {
        args = DOTS
        args$p = p
        do.call("qdist", args)
    }
    y <- QFUN(ppoints(length(x)))
    if (is.null(probs)) {
        probs <- c(0.01, 0.05, seq(0.1, 0.9, by = 0.1), 0.95, 
            0.99)
        if (length(x) >= 1000) 
            probs <- c(0.001, probs, 0.999)
    }
    qprobs <- QFUN(probs)
    plot(x, y, axes = FALSE, type = "n", ylim = range(c(y, qprobs)), 
        xlab = xlab, ylab = ylab, main=main)
    box()
    abline(h = qprobs, col = "grey")
    axis(1)
    axis(2, at = qprobs, labels = 100 * probs)
    points(x, y)
    xl <- quantile(x, c(0.25, 0.75))
    yl <- qdist(c(0.25, 0.75), ...)
    slope <- diff(yl)/diff(xl)
    int <- yl[1] - slope * xl[1]
    if (line) {
        abline(int, slope, col = "gray10")
    }
    z <- list(qdist = QFUN, int = int, slope = slope)
    class(z) <- "pplot"
    invisible(z)
}



if (FALSE) {

d <- NULL
d$time <- c( 94,  96,  99,  99, 104, 108, 112, 114, 117, 117, 118, 121, 121, 123, 129, 131, 133, 135, 136,
            139, 139, 140, 141, 141, 143, 144, 149, 149, 152, 153, 159, 159, 159, 159, 162, 168, 168, 169,
            170, 170, 171, 172, 173, 176, 177, 180, 180, 184, 187, 188, 189, 190, 196, 197, 203, 205, 211,
            213, 224, 226, 227, 256, 257, 269, 271, 274, 291, 300, 300, 300, 300, 300)
d$delta <- c( 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,
              1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,
              1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,
              1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   0,   0,   0,   0,   0)

mle.norm   <- tbs.survreg.mle(Surv(d$time,d$delta) ~ 1, dist="norm", method="Nelder-Mead",verbose=TRUE)
bayes.norm <- tbs.survreg.bayes(Surv(d$time,d$delta) ~ 1,dist="norm",
                               kick.beta=mle.norm$par[3],kick.lambda=mle.norm$par[1],kick.xi=mle.norm$par[2],
                               burn=500000,jump=2000,size=1000,scale=0.07)

bayes.norm <- tbs.survreg.bayes(Surv(d$time,d$delta) ~ 1,dist="norm",
                               kick.beta=mle.norm$par[3],kick.lambda=mle.norm$par[1],kick.xi=mle.norm$par[2],
                               burn=5,jump=2,size=10,scale=0.07)

####################################################
}

if (FALSE) {
#########################
## Gerando do modelo TBS

n.sample <- 10
dist     <- c("norm","t","cauchy","logistic","doubexp")
ni       <- 30 #c(30,100,1000)
cens     <- c(0,0.2,0.4,0.6)
par      <- matrix(NA,2,6)
par[1,]  <- c(1.5,0.5,1,1.5,2,-0.5)
par[2,]  <- c(0.5,2,0.5,-1.5,0.7,-0.9)

n.total <- length(dist)*length(cens)*length(par[,1])*length(ni)*n.sample
result.1  <- matrix(NA,n.total,1)
result.2  <- matrix(NA,n.total,20)

i.1 <- 1
initial.seed <- 1234
initial.time <- .gettime()
names(initial.time) <- NULL
run.time     <- 0
for (i.c in 1:length(cens)) {
  for (i.p in 1:length(par[,1])) {
    for (i.n in 1:length(ni)) {
      for (i.s in 1:n.sample) {
        seed <- (initial.seed+i.1)
        set.seed(seed)
        d    <- NULL
        x1   <- rnorm(ni[i.n],0,1)
        x2   <- sample(c(1,0,0),replace=TRUE,ni[i.n])
        x3   <- sample(c(1,1,1,0),replace=TRUE,ni[i.n])
        beta <- par[i.p,3:6]
        d$x  <- cbind(rep(1,ni[i.n]),x1,x2,x3)
        d$time  <- rtbs(ni[i.n],lambda=par[i.p,1],xi=par[i.p,2],beta=beta,x=d$x,dist="norm")
        d$delta <- rep(1,ni[i.n])
        censor  <- quantile(d$time,probs=(1-cens[i.c]))
        for (i in 1:ni[i.n]) {
          if (d$time[i] > censor) {
            d$delta[i] <- 0
            d$time[i]  <- censor
          }
        }

        for (i.d in 1:length(dist)) {
          mle.norm <- tbs.survreg.mle(Surv(d$time,d$delta) ~ d$x[,2] + d$x[,3] + d$x[,4],
                                      dist=dist[i.d], method="BFGS",verbose=TRUE)
          fit <- survreg(Surv(d$time,d$delta) ~ d$x[,2] + d$x[,3] + d$x[,4], dist="weibull")

          y1 <- ptbs(d$time[d$delta == 1],lambda=par[i.p,1],xi=par[i.p,2],
                     beta=par[i.p,3:6],x=d$x[d$delta == 1,],dist="norm")
          y3 <- pweibull(d$time[d$delta == 1],shape=1/fit$scale,
                         scale=exp(fit$linear.predictors[d$delta == 1]))
          result.1[i.1,1] <- dist[i.d]
          if (mle.norm$convergence) {
            y2 <- ptbs(d$time[d$delta == 1],lambda=mle.norm$par[1],xi=mle.norm$par[2],
                       beta=mle.norm$par[3:6],x=d$x[d$delta == 1,],dist=dist[i.d])
            result.2[i.1,] <- c(seed,cens[i.c],par[i.p,],ni[i.n],i.s,mle.norm$par,
                                max(abs(y1-y2)),mean(abs(y1-y2)),
                                max(abs(y1-y3)),mean(abs(y1-y3)))
          } else {
            result.2[i.1,] <- c(seed,cens[i.c],par[i.p,],ni[i.n],i.s,rep(NA,length(par[i.p,])),
                                NA,NA,max(abs(y1-y3)),mean(abs(y1-y3)))
          }
           
          if ((i.1 %% 1) == 0) {
            run.time2 <- run.time
            run.time  <- .gettime() - initial.time
            names(run.time) <- NULL
            print(paste("i: ",i.1,
                        "| resta: ",n.total-i.1,
                        "| g1: ",round(run.time,2),
                        "| g2: ",round(run.time-run.time2,2),
                        "| ",date(),sep=""))
          }
          i.1 <- i.1+1
        }
      }
    }
  }
}
colnames(result.1) <- "dist"
result.2 <- cbind(result.2[,1:20],ifelse(result.2[,17] < result.2[,19],1,0),
                  ifelse(result.2[,18] < result.2[,20],1,0),
                  result.2[,17]-result.2[,19],result.2[,18]-result.2[,20])
colnames(result.2) <- c("seed","censor","lambda","xi","beta0","beta1","beta2",
                        "beta3","n.sample","sample","est.lambda","est.xi",
                        "est.beta0","est.beta1","est.beta2","est.beta3","max",
                        "mean","max.weib","mean.weib","comp.max","comp.mean",
                        "dif.max","dif.mean")

result <- data.frame(cbind(result.1,result.2))

print(table(result$comp.mean,result$dist,result$censor,result$n.sample))
boxplot(result.2[result$dist=="norm",23:24])

#result[((result$dist == "norm") & (result$n.sample==1000) &
#        (result$censor == 0) & (result$lambda == par[1,1])),c(18,20)]
}

###################################################
if (TRUE) {


#########################
## Gerando da weibull

n.sample <- 2
dist     <- c("norm","t","cauchy","logistic","doubexp")
ni       <- c(30,100,1000)
cens     <- c(0,0.2,0.4,0.6)
par      <- matrix(NA,4,2)
par[1,]  <- c(1/0.5,log(1.0))
par[2,]  <- c(1/3.5,log(0.5))
par[3,]  <- c(1/3.5,log(5.0))
par[4,]  <- c(1/5.0,log(1.5))

#median   <- c(exp(par[1,2])*log(2)^(par[1,1]),
#              exp(par[2,2])*log(2)^(par[2,1]),
#              exp(par[3,2])*log(2)^(par[3,1]),
#              exp(par[4,2])*log(2)^(par[4,1]))

n.total <- length(dist)*length(cens)*length(par[,1])*length(ni)*n.sample
result.1  <- matrix(NA,n.total,1)
result.2  <- matrix(NA,n.total,15)

i.1 <- 1
initial.seed <- 1234
initial.time <- .gettime()
names(initial.time) <- NULL
run.time     <- 0
for (i.c in 1:length(cens)) {
  for (i.p in 1:length(par[,1])) {
    for (i.n in 1:length(ni)) {
      for (i.s in 1:n.sample) {
        seed <- (initial.seed+i.1)
        set.seed(seed)
        d    <- NULL
        beta <- par[i.p,2]
        d$time  <- rweibull(ni[i.n],shape=1/par[i.p,1],scale=exp(beta))
        d$delta <- rep(1,ni[i.n])
        censor  <- quantile(d$time,probs=(1-cens[i.c]))
        for (i in 1:ni[i.n]) {
          if (d$time[i] > censor) {
            d$delta[i] <- 0
            d$time[i]  <- censor
          }
        }

        fit <- survreg(Surv(d$time,d$delta) ~ 1, dist="weibull")
        y3 <- pweibull(d$time[d$delta == 1],shape=1/fit$scale, scale=exp(fit$coefficients))
        for (i.d in 1:length(dist)) {
          mle.norm <- tbs.survreg.mle(Surv(d$time,d$delta) ~ 1,dist=dist[i.d], method="Rsolnp",verbose=TRUE)

          y1 <- pweibull(d$time[d$delta == 1],shape=1/par[i.p,1],scale=exp(par[i.p,2]))
          result.1[i.1,1] <- dist[i.d]
          if (mle.norm$convergence) {
            y2 <- ptbs(d$time[d$delta == 1],lambda=mle.norm$par[1],xi=mle.norm$par[2],
                       beta=mle.norm$par[3],dist=dist[i.d])
            result.2[i.1,] <- c(seed,cens[i.c],par[i.p,],ni[i.n],i.s,mle.norm$par,
                                max(abs(y1-y2)),mean(abs(y1-y2)),
                                max(abs(y1-y3)),mean(abs(y1-y3)),
                                exp(mle.norm$par[3]),exp(fit$coefficients)*log(2)^fit$scale)
          } else {
            result.2[i.1,] <- c(seed,cens[i.c],par[i.p,],ni[i.n],i.s,rep(NA,(length(par[i.p,])+1)),
                                NA,NA,max(abs(y1-y3)),mean(abs(y1-y3)),
                                NA,exp(fit$coefficients)*log(2)^fit$scale)
          }
           
          if ((i.1 %% 1) == 0) {
            run.time2 <- run.time
            run.time  <- .gettime() - initial.time
            names(run.time) <- NULL
            cat("step: ",i.1,
                "| remains: ",n.total-i.1,
                "| ll: ",mle.norm$log.lik,
                "| g1: ",round(run.time,2),
                "| g2: ",round(run.time-run.time2,2),
                "| ",date(),'\n',sep="")
          }
          i.1 <- i.1+1
        }
      }
    }
  }
}
colnames(result.1) <- "dist"
result.2 <- cbind(result.2[,1:15],ifelse(result.2[,10] < result.2[,12],1,0),
                  ifelse(result.2[,11] < result.2[,13],1,0),
                  result.2[,10]-result.2[,12],result.2[,11]-result.2[,13])
colnames(result.2) <- c("seed","censor","gamma","beta","n.sample","sample",
                        "est.lambda","est.xi","est.beta","max","mean",
                        "max.weib","mean.weib","median","weib.median",
                        "comp.max","comp.mean","dif.max","dif.mean")
result <- data.frame(cbind(result.1,result.2))


n.total <- length(dist)*length(cens)*length(par[,1])*length(ni)
result.median1  <- matrix(NA,n.total,1)
result.median2  <- matrix(NA,n.total,9)
i <- 1
for (i.c in 1:length(cens)) {
  for (i.p in 1:length(par[,1])) {
    for (i.n in 1:length(ni)) {
      median <- exp(par[i.p,2])*log(2)^(par[i.p,1])
      for (i.d in 1:length(dist)) {
        p1 <- apply(result.2[((result.1[,1] == dist[i.d]) & (result.2[,2] == cens[i.c]) &
                         (result.2[,4] == par[i.p,2]) & (result.2[,5] == ni[i.n])),14:15],2,mean)
        p2 <- apply((result.2[((result.1[,1] == dist[i.d]) & (result.2[,2] == cens[i.c]) &
                         (result.2[,4] == par[i.p,2]) & (result.2[,5] == ni[i.n])),14:15]-
                     exp(par[i.p,2])*log(2)^(par[i.p,1]))^2,2,mean)
        result.median1[i,] <- dist[i.d]
        result.median2[i,] <- c(cens[i.c],par[i.p,],ni[i.n],median,p1,p2)
        i <- i+1
      }
    }
  }
}
colnames(result.median1) <- "dist"
colnames(result.median2) <- c("censor","gamma","beta","n.sample","orig.median",
                              "est.median","weib.median","mse.est","weib.mse")
result.median <- data.frame(cbind(result.median1,result.median2))

table(result$comp.mean,result$dist,result$censor,result$n.sample)
boxplot(result.2[result$dist=="norm",18:19])


#result[((result$dist == "norm") & (result$n.sample==1000) &
#        (result$censor == 0) & (result$lambda == par[1,1])),c(17,19)]

###################################################
}

if (FALSE) {

n.sample <- 1
dist     <- c("cauchy")
ni       <- c(30)
cens     <- c(0)
par      <- matrix(NA,1,6)
par[1,]  <- c(0.5,2,0.5,-1.5,0.7,-0.9)

n.total <- length(dist)*length(cens)*length(par[,1])*length(ni)*n.sample
result.1  <- matrix(NA,n.total,1)
result.2  <- matrix(NA,n.total,18)

i.1 <- 1
for (i.c in 1:length(cens)) {
  for (i.p in 1:length(par[,1])) {
    for (i.n in 1:length(ni)) {
      for (i.s in 1:n.sample) {
        seed <- (7670)
        set.seed(seed)
        d    <- NULL
        x1   <- rnorm(ni[i.n],0,1)
        x2   <- sample(c(1,0,0),replace=TRUE,ni[i.n])
        x3   <- sample(c(1,1,1,0),replace=TRUE,ni[i.n])
        beta <- par[i.p,3:6]
        d$x  <- cbind(rep(1,ni[i.n]),x1,x2,x3) 
        d$time  <- rtbs(ni[i.n],lambda=par[i.p,1],xi=par[i.p,2],beta=beta,x=d$x,dist="norm")
        d$delta <- rep(1,ni[i.n])
        censor  <- quantile(d$time,probs=(1-cens[i.c]))
        for (i in 1:ni[i.n]) {
          if (d$time[i] > censor) {
            d$delta[i] <- 0
            d$time[i]  <- censor
          }
        }

        for (i.d in 1:length(dist)) {
          mle.norm <- tbs.survreg.mle(Surv(d$time,d$delta) ~ d$x[,2] + d$x[,3] + d$x[,4],
                                      dist=dist[i.d], method="BFGS",verbose=TRUE)

          y1 <- ptbs(d$time[d$delta == 1],lambda=par[i.p,1],xi=par[i.p,2],
                     beta=par[i.p,3:6],x=d$x[d$delta == 1,],dist="norm")
          y2 <- ptbs(d$time[d$delta == 1],lambda=mle.norm$par[1],xi=mle.norm$par[2],
                     beta=mle.norm$par[3:6],x=d$x[d$delta == 1,],dist=dist[i.d])
          result.1[i.1,1] <- dist[i.d]
          result.2[i.1,] <- c(seed,cens[i.c],par[i.p,],ni[i.n],i.s,mle.norm$par,
                             max(abs(y1-y2)),mean(abs(y1-y2)))
 
#          if ((i.1 %% 100) == 0)
#            print(n.total-i.1)
            print(i.1)
          i.1 <- i.1+1
        }
      }
    }
  }
}
colnames(result.1) <- "dist"
colnames(result.2) <- c("seed","censor","lambda","xi","beta0","beta1",
                        "beta2","beta3","n.sample","sample","est.lambda","est.xi","est.beta0",
                        "est.beta1","est.beta2","est.beta3","max","mean")
result <- data.frame(cbind(result.1,result.2))

dtbs(d$time[2],lambda=mle.norm$par[1],xi=mle.norm$par[2],beta=mle.norm$par[3:6],x=d$x[2,],dist=dist[i.d])
dtbs(d$time[2],lambda=347.356925,xi=260.00285,
     beta=c(14.315272,-43.463057,3.291133,12.8141),x=d$x[2,],dist=dist[i.d])


###################################################



bayes.norm <- tbs.survreg.bayes(Surv(d$time,d$delta) ~ x1 + x2 + x3,dist="norm",
                               kick.beta=mle.norm$par[3:6],kick.lambda=mle.norm$par[1],kick.xi=mle.norm$par[2],
                               burn=500000,jump=2000,size=1000,scale=0.07)



set.seed(1234)
n    <- 1000
d    <- NULL
beta <- 0.5
d$x  <- rep(1,n) 
d$time  <- sort(rtbs(n,lambda=1.5,xi=0.5,beta=beta,dist="norm"))
d$delta <- rep(1,n)
d$time[(n-floor(n*0.2)):n]  <- rep(d$time[(n-floor(n*0.2))-1],length(d$time[(n-floor(n*0.2)):n]))
d$delta[(n-floor(n*0.2)):n] <- rep(0,length(d$time[(n-floor(n*0.2)):n]))
print(sum(d$delta)/n)

mle.norm2  <- tbs.survreg.mle(Surv(d$time,d$delta) ~ 1, dist="norm", method="Nelder-Mead",verbose=TRUE)
mle.norm2$par

eixo.y1 <- ptbs(d$time,lambda=1.5,xi=0.5,beta=beta,dist="norm")
eixo.y2 <- ptbs(d$time,lambda=mle.norm2$par[1],xi=mle.norm2$par[2],beta=mle.norm2$par[3],dist="norm")
plot(d$time,eixo.y2,type="l",ylim=c(0,1))
lines(d$time,eixo.y1,type="l",col=2,lty=2)
}

##################################################
if (FALSE) {

n.sample <- 1
dist     <- "norm"
ni       <- 1000
cens     <- 0.6
par      <- matrix(NA,1,2)
par[1,]  <- c(1/3.5,log(0.5))

n.total <- length(dist)*length(cens)*length(par[,1])*length(ni)*n.sample
result.1  <- matrix(NA,n.total,1)
result.2  <- matrix(NA,n.total,15)

i.1 <- 1
initial.seed <- 1644
initial.time <- .gettime()
names(initial.time) <- NULL
run.time     <- 0
for (i.c in 1:length(cens)) {
  for (i.p in 1:length(par[,1])) {
    for (i.n in 1:length(ni)) {
      for (i.s in 1:n.sample) {
        seed <- (initial.seed+i.1)
        set.seed(seed)
        d    <- NULL
        beta <- par[i.p,2]
        d$time  <- rweibull(ni[i.n],shape=1/par[i.p,1],scale=exp(beta))
        d$delta <- rep(1,ni[i.n])
        censor  <- quantile(d$time,probs=(1-cens[i.c]))
        for (i in 1:ni[i.n]) {
          if (d$time[i] > censor) {
            d$delta[i] <- 0
            d$time[i]  <- censor
          }
        }

        fit <- survreg(Surv(d$time,d$delta) ~ 1, dist="weibull")
        y3 <- pweibull(d$time[d$delta == 1],shape=1/fit$scale, scale=exp(fit$coefficients))
        for (i.d in 1:length(dist)) {
          mle.norm <- tbs.survreg.mle(Surv(d$time,d$delta) ~ 1,dist=dist[i.d], method="BFGS",verbose=TRUE)

          y1 <- pweibull(d$time[d$delta == 1],shape=1/par[i.p,1],scale=exp(par[i.p,2]))
          result.1[i.1,1] <- dist[i.d]
          if (mle.norm$convergence) {
            y2 <- ptbs(d$time[d$delta == 1],lambda=mle.norm$par[1],xi=mle.norm$par[2],
                       beta=mle.norm$par[3],dist=dist[i.d])
            result.2[i.1,] <- c(seed,cens[i.c],par[i.p,],ni[i.n],i.s,mle.norm$par,
                                max(abs(y1-y2)),mean(abs(y1-y2)),
                                max(abs(y1-y3)),mean(abs(y1-y3)),
                                exp(mle.norm$par[3]),exp(fit$coefficients)*log(2)^fit$scale)
          } else {
            result.2[i.1,] <- c(seed,cens[i.c],par[i.p,],ni[i.n],i.s,rep(NA,(length(par[i.p,])+1)),
                                NA,NA,max(abs(y1-y3)),mean(abs(y1-y3)),
                                NA,exp(fit$coefficients)*log(2)^fit$scale)
          }
           
          if ((i.1 %% 1) == 0) {
            run.time2 <- run.time
            run.time  <- .gettime() - initial.time
            names(run.time) <- NULL
            print(paste("i: ",i.1,
                        "| resta: ",n.total-i.1,
                        "| g1: ",round(run.time,2),
                        "| g2: ",round(run.time-run.time2,2),
                        "| ",date(),sep=""))
          }
          i.1 <- i.1+1
        }
      }
    }
  }
}
#colnames(result.1) <- "dist"
#result.2 <- cbind(result.2[,1:15],ifelse(result.2[,10] < result.2[,12],1,0),
#                  ifelse(result.2[,11] < result.2[,13],1,0),
#                  result.2[,10]-result.2[,12],result.2[,11]-result.2[,13])
#colnames(result.2) <- c("seed","censor","gamma","beta","n.sample","sample",
#                        "est.lambda","est.xi","est.beta","max","mean",
#                        "max.weib","mean.weib","median","weib.median",
#                        "comp.max","comp.mean","dif.max","dif.mean")
#result <- data.frame(cbind(result.1,result.2))





}




}
