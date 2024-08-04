library(MASS)

uncertaintysim <- function(samples, rho, k, psi, iter){
  vehwlist = c()
  gammalist = c()
  ehwlist = c()
  thetalist1 = c()
  thetalist2 = c()
  thetalist3 = c()
  vdescrlist=c()
  vcaussamplist = c()
  for (j in 1:iter){
  #Create matrix to fill with attributes
  z= matrix(, nrow= samples, ncol= k+1 )
  #Create column of 1s
  z[,1] <- 1
  #Generate other attributes from a standard normal
  for (i in 1:k){
    z[, i+1] <- rnorm(samples, 0,1) 
  }
  theta = rnorm(samples, z[, 2]*psi,1)
  #generate xi
  xi = rnorm(samples, 0,1)
  #generate u; note that in this setup, u =x  
  u = rnorm(samples, 0,1)
  #generate potential outcomes
  y= u*theta + xi
  sampleindex <- sample(samples, samples*rho)
  
  zsample = z[sampleindex,]
  usample = u[sampleindex]
  thetasample = theta[sampleindex]
  lambdahat =  (t(usample)%*%zsample)%*%solve(t(zsample)%*%zsample)
  uresid = usample -  lambdahat%*%t(zsample)
  gammahat= (1/(samples*rho))*(uresid)%*%t(uresid)
  bigdat <- cbind(z, y)
  datsample = bigdat[sampleindex,]
  datsample = cbind(datsample, t(uresid))
  dfsample = data.frame(datsample)
  
  reg1 <- lm(y ~ ., data=subset(dfsample, select = -c(V1)))
  eresid = reg1$residuals
  e2= eresid*eresid
  ehw = (1/(samples*rho))*(e2*uresid)%*%t(uresid)
 
  
  vehw = (1/gammahat)*ehw*(1/gammahat)
  vehwlist[j] <- (vehw/(rho*samples))^{1/2}
  gammalist[j] <- gammahat
  ehwlist[j] <- ehw
  
  vdescrlist[j] <-((1- rho)*vehw/(rho*samples))^{1/2}
  
  bigdat <- cbind(bigdat, u)
  bigdf <- data.frame(bigdat)
  reg2 <- lm(y ~ ., data=subset(bigdf, select = -c(V1)))
  thetahat = reg1$coefficients[k+2]
  thetades = reg2$coefficients[k+2]
  thetalist1[j] <- thetahat- thetades
  
  thetacausal = (1/samples)*sum(theta)
  thetalist2[j] <- thetahat- thetacausal
  
  thetacausalsamp <- (1/(samples*rho))*sum(thetasample)
  thetalist3[j] <- thetahat - thetacausalsamp
  
  g = (uresid*eresid)%*%zsample%*%solve(t(zsample)%*%zsample)
  deltazhat = (1/(rho*samples))*(uresid*eresid - g%*%t(zsample))%*%t(uresid*eresid - g%*%t(zsample))
  vcaussamp = (1/gammahat)*deltazhat*(1/gammahat)
  vcaussamplist[j] <- (vcaussamp/(rho*samples))^{1/2}
  }
  return(c(mean(vehwlist), mean(gammalist), mean(ehwlist), mean(vdescrlist),mean(vcaussamplist), sd(thetalist1),sd(thetalist2), sd(thetalist3)))
}



check <- uncertaintysim(100000, .01, 1, 2, iter= 1000)

##check succesful for k =1

check <- uncertaintysim(100000, .01, 10, 2, iter= 1000)

##check succesful for k = 10

check <- uncertaintysim(100000, .01, 1, 2, iter= 50000)
check <- uncertaintysim(100000, .01, 10, 2, iter= 50000)

n= 10000
k= 10
z= matrix(, nrow= n, ncol= k+1 )
z[,1] <- 1

for(i in 1:n*rho){
  altehw = 0
  altehw = altehw + e2[i]*uresid[i]*uresid[i]
}

for (i in 1:k){
  z[, i+1] <- rnorm(n, 0,1) 
}
psi= 2
theta= z[,2]*psi
xi = rnorm(n, 0,1)
u = rnorm(n, 0,1)
y= u*theta + xi

#sample 
rho= .1
sampleindex <- sample(n, n*rho)

zsample = z[sampleindex,]
usample = u[sampleindex]
lambdahat =  (t(usample)%*%zsample)%*%solve(t(zsample)%*%zsample)
uresid = usample -  lambdahat%*%t(zsample)
gammahat= (1/(n*rho))*(uresid)%*%t(uresid)
bigdat <- cbind(z, y)
datsample = bigdat[sampleindex,]
datsample = cbind(datsample, t(uresid))
dfsample = data.frame(datsample)

reg1 <- lm(y ~ ., data=subset(dfsample, select = -c(V1)))
eresid = reg1$residuals
e2= eresid*eresid
ehw = (1/(n*rho))*(e2*uresid)%*%t(uresid)
vehw = (1/gammahat)*ehw*(1/gammahat)
finalehw = (vehw/(rho*n))^{1/2}
bigdat <- cbind(bigdat, u)
bigdf <- data.frame(bigdat)
reg2 <- lm(y ~ ., data=subset(bigdf, select = -c(V1)))
thetahat = reg$coefficients[12]
thetadesc = reg2$coefficients[12]
thetacausal = (1/n)*sum(theta)

g = (uresid*eresid)%*%zsample%*%solve(t(zsample)%*%zsample)
deltazhat = (1/(rho*n))*(uresid*eresid - g%*%t(zsample))%*%t(uresid*eresid - g%*%t(zsample))
vcaussamp = (1/gammahat)*deltazhat*(1/gammahat)

##newsim

covsim <- function(samples, rho, k, psi, iter){
  vehwlist = c()
  gammalist = c()
  ehwlist = c()
  thetalist1 = c()
  thetalist2 = c()
  thetalist3 = c()
  vdescrlist=c()
  vcaussamplist = c()
  for (j in 1:iter){
    #Create matrix to fill with attributes
    z= matrix(, nrow= samples, ncol= k+1 )
    #Create column of 1s
    z[,1] <- 1
    #Generate other attributes from a standard normal
    for (i in 1:k){
      z[, i+1] <- rnorm(samples, 0,1) 
    }
    theta = rnorm(samples, z[, 2]*psi,1)
    #generate xi
    xi = rnorm(samples, 0,1)
    #generate u; note that in this setup, u =x  
    u = rnorm(samples, 0,1)
    #generate potential outcomes
    y= u*theta + xi
    sampleindex <- sample(samples, samples*rho)
    
    zsample = z[sampleindex,]
    usample = u[sampleindex]
    thetasample = theta[sampleindex]
    lambdahat =  (t(usample)%*%zsample)%*%solve(t(zsample)%*%zsample)
    uresid = usample -  lambdahat%*%t(zsample)
    gammahat= (1/(samples*rho))*(uresid)%*%t(uresid)
    bigdat <- cbind(z, y)
    datsample = bigdat[sampleindex,]
    datsample = cbind(datsample, t(uresid))
    dfsample = data.frame(datsample)
    
    reg1 <- lm(y ~ ., data=subset(dfsample, select = -c(V1)))
    eresid = reg1$residuals
    e2= eresid*eresid
    ehw = (1/(samples*rho))*(e2*uresid)%*%t(uresid)
    
    
    vehw = (1/gammahat)*ehw*(1/gammahat)
    vehwlist[j] <- (vehw/(rho*samples))^{1/2}
    gammalist[j] <- gammahat
    ehwlist[j] <- ehw
    
    vdescrlist[j] <-((1- rho)*vehw/(rho*samples))^{1/2}
    
    bigdat <- cbind(bigdat, u)
    bigdf <- data.frame(bigdat)
    reg2 <- lm(y ~ ., data=subset(bigdf, select = -c(V1)))
    thetahat = reg1$coefficients[k+2]
    thetades = reg2$coefficients[k+2]
    thetalist1[j] <- thetahat- thetades
    
    thetacausal = (1/samples)*sum(theta)
    thetalist2[j] <- thetahat- thetacausal
    
    thetacausalsamp <- (1/(samples*rho))*sum(thetasample)
    thetalist3[j] <- thetahat - thetacausalsamp
    
    g = (uresid*eresid)%*%zsample%*%solve(t(zsample)%*%zsample)
    deltazhat = (1/(rho*samples))*(uresid*eresid - g%*%t(zsample))%*%t(uresid*eresid - g%*%t(zsample))
    vcaussamp = (1/gammahat)*deltazhat*(1/gammahat)
    vcaussamplist[j] <- (vcaussamp/(rho*samples))^{1/2}
  }
  return(c(mean(vehwlist), mean(gammalist), mean(ehwlist), mean(vdescrlist),mean(vcaussamplist), sd(thetalist1),sd(thetalist2), sd(thetalist3)))
}

samples = 10000

data = mvrnorm(n =samples, mu = c(0,0), Sigma = matrix(c(1, u1u2, u1u2, 1), nrow=2), empirical = TRUE)

  #Create matrix to fill with attributes
  z= matrix(, nrow= samples, ncol= k+1 )
  #Create column of 1s
  z[,1] <- 1
  #Generate other attributes from a standard normal
  for (i in 1:k){
    z[, i+1] <- rnorm(samples, 0,1) 
  }
  theta = rnorm(samples, z[, 2]*psi,1)
  #generate xi
  xi = rnorm(samples, 0,1)
  #generate u; note that in this setup, u =x  
  u = rnorm(samples, 0,1)
  #generate potential outcomes
  y= u*theta + xi
  sampleindex <- sample(samples, samples*rho)
  
  zsample = z[sampleindex,]
  usample = u[sampleindex]
  thetasample = theta[sampleindex]
  lambdahat =  (t(usample)%*%zsample)%*%solve(t(zsample)%*%zsample)
  uresid = usample -  lambdahat%*%t(zsample)
  gammahat= (1/(samples*rho))*(uresid)%*%t(uresid)
  bigdat <- cbind(z, y)
  datsample = bigdat[sampleindex,]
  datsample = cbind(datsample, t(uresid))
  dfsample = data.frame(datsample)
  
  reg1 <- lm(y ~ ., data=subset(dfsample, select = -c(V1)))
  eresid = reg1$residuals
  e2= eresid*eresid
  ehw = (1/(samples*rho))*(e2*uresid)%*%t(uresid)
  
  
  vehw = (1/gammahat)*ehw*(1/gammahat)
  vehwlist[j] <- (vehw/(rho*samples))^{1/2}
  gammalist[j] <- gammahat
  ehwlist[j] <- ehw
  
  vdescrlist[j] <-((1- rho)*vehw/(rho*samples))^{1/2}
  
  bigdat <- cbind(bigdat, u)
  bigdf <- data.frame(bigdat)
  reg2 <- lm(y ~ ., data=subset(bigdf, select = -c(V1)))
  thetahat = reg1$coefficients[k+2]
  thetades = reg2$coefficients[k+2]
  thetalist1[j] <- thetahat- thetades
  
  thetacausal = (1/samples)*sum(theta)
  thetalist2[j] <- thetahat- thetacausal
  
  thetacausalsamp <- (1/(samples*rho))*sum(thetasample)
  thetalist3[j] <- thetahat - thetacausalsamp
  
  g = (uresid*eresid)%*%zsample%*%solve(t(zsample)%*%zsample)
  deltazhat = (1/(rho*samples))*(uresid*eresid - g%*%t(zsample))%*%t(uresid*eresid - g%*%t(zsample))
  vcaussamp = (1/gammahat)*deltazhat*(1/gammahat)
