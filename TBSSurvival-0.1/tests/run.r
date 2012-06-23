dist <- c("norm","t","logistic","cauchy","doubexp")
n    <- c(30,100,1000)
cens <- c(0,0.2,0.4,0.6)
par  <- matrix(NA,4,2)
par[1,] <- c(1,1)
par[2,] <- c(2,1)
par[3,] <- c(1,2)
par[4,] <- c(2,2)
nn = 1
text <- ""
##text <- "initial.time <- proc.time()[3]/60\ncat('starting at ',initial.time,'\\n')\n"
for (i.d in 1:length(dist)) {
  for (i.n in 1:length(n)) {
    for (i.c in 1:length(cens)) {
      for (i.p in 1:length(par[,1])) {
        metodo <- paste("sim.weib(gen=list(n=",n[i.n],",cens=",cens[i.c],",scale=",par[i.p,1],",shape=",par[i.p,2],
                      "),est=list(dist=\"",dist[i.d],"\"),n.copies=1000,method=\"Rsolnp\",prefix=\"/project/s221/cassio/tbs/oe/\")",sep="")
        text <- paste(text,"cat(",nn,",'",metodo," '); ",metodo,"\n",sep="")
        nn <- nn +1
      }
    }
  }
}

par  <- matrix(NA,18,3)
par[1,] <- c(0.5,0.5,1)
par[2,] <- c(0.5,1.0,1)
par[3,] <- c(0.5,2.0,1)
par[4,] <- c(1.0,0.5,1)
par[5,] <- c(1.0,1.0,1)
par[6,] <- c(1.0,2.0,1)
par[7,] <- c(2.0,0.5,1)
par[8,] <- c(2.0,1.0,1)
par[9,] <- c(2.0,2.0,1)
par[10,] <- c(0.5,0.5,5)
par[11,] <- c(0.5,1.0,5)
par[12,] <- c(0.5,2.0,5)
par[13,] <- c(1.0,0.5,5)
par[14,] <- c(1.0,1.0,5)
par[15,] <- c(1.0,2.0,5)
par[16,] <- c(2.0,0.5,5)
par[17,] <- c(2.0,1.0,5)
par[18,] <- c(2.0,2.0,5)
nn=1
for (i.d in 1:length(dist)) {
  for (i.n in 1:length(n)) {
    for (i.c in 1:length(cens)) {
      for (i.p in 1:length(par[,1])) {
        metodo <- paste("sim.tbs(gen=list(n=",n[i.n],",cens=",cens[i.c],",lambda=",par[i.p,1],",xi=",par[i.p,2],
                        ",beta=",par[i.p,3],",dist=\"",dist[i.d],"\"),n.copies=1000,method=\"Rsolnp\",prefix=\"/project/s221/cassio/tbs/oe/\")",sep="")
        text <- paste(text,"cat(",nn,",'",metodo," '); ",metodo,"\n",sep="")
        nn=nn+1
      }
    }
  }
}

##text <- paste(text,"run.time <- proc.time()[3]/60-initial.time\ncat('end at ',run.time,'\\n')\n",sep="")

cat(text,file="run-all.r",append=FALSE,sep="")
