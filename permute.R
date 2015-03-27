source("aibp2.R")
source("countdown.R")
library(doMC) 
registerDoMC(system("nproc",intern=TRUE))

calc.settings <- function(Zs) {
  Z11 <- lapply(Zs,function(z) if (sum(z[1,])==2) z)
  Z11 <- Z11[!unlist(lapply(Z11, is.null))]
  EZ11 <- sum.matrices(Z11)/length(Z11)
  u.z11 <- unique.matrix(Z11)

  Z11.0011 <- lapply(Z11,function(z) if (all(z[2,1:2]==c(0,0)) && sum(z[2,])==2) z)
  Z11.0011 <- Z11.0011[!unlist(lapply(Z11.0011,is.null))]
  EZ11.0011 <- sum.matrices(Z11.0011)/length(Z11.0011)
  print(length(Z11.0011))

  EZ11.0011
}

B <- 1e4
D <- matrix(9,5,5)
diag(D) <- 0
D[1,3] <- D[3,1] <- D[2,5] <- D[5,2] <- 1

Zo <- foreach(i=1:B) %dopar% raibp(a=2,N=5,l=function(s,t,d) exp(-d[s,t]))
Zs <- foreach(i=1:B) %dopar% raibp(a=2,D=D,l=function(s,t,d) exp(-d[s,t]))
Zp <- foreach(i=1:B) %dopar% raibp(a=2,D=D,l=function(s,t,d) exp(-d[s,t]),perm=TRUE)

EO <- sum.matrices(Zo) / B
EZ <- sum.matrices(Zs) / B
EP <- sum.matrices(Zp) / B

par(mfrow=c(3,1))
  a.image(EO,numbers=TRUE,main="IBP")
  a.image(EZ,numbers=TRUE,main="Not Permuted")
  a.image(EP,numbers=TRUE,main="Permuted")
par(mfrow=c(1,1))

EOc <- calc.settings(Zo) #128
EZc <- calc.settings(Zs) #128
EPc <- calc.settings(Zp) #165

par(mfrow=c(3,1))
  a.image(round(EOc,4),numbers=TRUE,main="IBP")
  a.image(round(EZc,4),numbers=TRUE,main="Not Permuted")
  a.image(round(EPc,4),numbers=TRUE,main="Permuted")
par(mfrow=c(1,1))

apply(EZc,1,sum)
apply(EPc,1,sum)
apply(EZp,1,sum)
