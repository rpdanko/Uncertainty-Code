setwd("C:/Users/chris/Documents/Uncertainty Code")

library(MASS)
library(AER)
library(brew)

set.seed(48824)

uncertaintysim <- function(samples, rho, k, psi, iter, nu, delta, sigma2=1){
  vehwlist = c()
  gammalist = c()
  ehwlist = c()
  thetalist1 = c()
  thetalist2 = c()
  thetalist3 = c()
  vdescrlist=c()
  vcaussamplist = c()
  vcauslist = c()
  thetacausalcoverage = c()
  thetaehwcoverage = c()
  olscausalcoverage = c()
  olsehwcoverage = c()
  #Create matrix to fill with attributes
  x= matrix(, nrow= samples, ncol= k+1 )
  #Create column of 1s
  x[,1] <- 1
  #Generate other attributes from a standard normal
  for (i in 1:k){
    x[, i+1] <- rnorm(samples, 0,1) 
  }
  theta = rnorm(samples,x[,2:(k+1)]%*%matrix(psi),sigma2)
  #generate xi
  xi = rnorm(samples, 0,1)
  #generate w; note that in this setup, w =u  
  z = rnorm(samples, 0, 1)
  w = delta*z + nu*xi
  #generate potential outcomes
  y= w*theta + xi
   for (j in 1:iter){
    
    sampleindex <- sample(samples, samples*rho)
    
    xsample = x[sampleindex,]
    wsample = w[sampleindex]
    zsample = z[sampleindex]
    ysample = y[sampleindex]
    thetasample = theta[sampleindex]
    lambdahat =  (t(wsample)%*%xsample)%*%solve(t(xsample)%*%xsample)
    wresid = wsample -  lambdahat%*%t(xsample)
    bigxihat =  (t(zsample)%*%xsample)%*%solve(t(xsample)%*%xsample)
    zresid = zsample -  lambdahat%*%t(xsample)
    gammahat= (1/(samples*rho))*(zresid)%*%t(wresid)
    bigdat <- cbind(y,x)
    datsample = bigdat[sampleindex,]
    dfsample = data.frame(datsample)
    
    reg1 <- lm(ysample ~ t(wresid) + xsample, data=subset(dfsample, select = -c(V2)))
    ivreg1 <- ivreg(ysample ~ t(wresid) + xsample | t(zresid) + xsample)
    eresid = ivreg1$residuals
    e2= eresid*eresid
    ehw = (1/(samples*rho))*(e2*zresid)%*%t(zresid)
    olstheta <- reg1$coefficients[2]
    
    vehw = (1/gammahat)*ehw*(1/gammahat)
    
    sdehw <- (vehw/(rho*samples))^{1/2}
    
    vehwlist[j] <- (vehw/(rho*samples))^{1/2}
    gammalist[j] <- gammahat
    ehwlist[j] <- ehw
    
    vdescrlist[j] <-((1- rho)*vehw/(rho*samples))^{1/2}
    
    
    
    
    bigdf <- data.frame(bigdat)
    ##calculate thetades
    regdes <- ivreg(y ~ w + x | z+ x)
    thetahat = ivreg1$coefficients[2]
    thetades = regdes$coefficients[2]
    thetalist1[j] <- thetahat- thetades
    
    
    
    thetacausal = (1/samples)*sum(theta)
    thetalist2[j] <- thetahat- thetacausal
    
    thetacausalsamp <- (1/(samples*rho))*sum(thetasample)
    thetalist3[j] <- thetahat - thetacausalsamp
    
    g = (zresid*eresid)%*%xsample%*%solve(t(xsample)%*%xsample)
    deltazhat = (1/(rho*samples))*(zresid*eresid - g%*%t(xsample))%*%t(zresid*eresid - g%*%t(xsample))
    vcaussamp = (1/gammahat)*deltazhat*(1/gammahat)
    vcaussamplist[j] <- (vcaussamp/(rho*samples))^(1/2)
    
    vcaus = rho*vcaussamp + (1-rho)*vehw
    
    ##calculate standard errors
    
    sdcaus <- (vcaus/(rho*samples))^(1/2) 
    
    # Collect coverages
    
    if(thetacausal < thetahat + 1.96*sdcaus & thetacausal > thetahat - 1.96*sdcaus ){
      thetacausalcoverage[j] <- 1
    } else {
      thetacausalcoverage[j] <- 0
    }
    
    if(thetacausal < thetahat + 1.96*sdehw & thetacausal > thetahat - 1.96*sdehw ){
      thetaehwcoverage[j] <- 1
    } else {
      thetaehwcoverage[j] <- 0
    }
    
    if(olstheta < thetahat + 1.96*sdcaus & olstheta > thetahat - 1.96*sdcaus ){
      olscausalcoverage[j] <- 1
    } else {
      olscausalcoverage[j] <- 0
    }
    
    if(olstheta < thetahat + 1.96*sdehw & olstheta > thetahat - 1.96*sdehw ){
      olsehwcoverage[j] <- 1
    } else { 
      olsehwcoverage[j] <- 0
    }
    
    
    
    
    
    
    #Collect standard errors
    
    vcauslist[j] <- sdcaus
    
    
    
    
  }
  return(c(mean(vehwlist), mean(gammalist), mean(ehwlist), mean(vdescrlist),
           mean(vcaussamplist), sd(thetalist1),sd(thetalist2), sd(thetalist3),
           mean(vcauslist)/mean(vehwlist),
           mean(vcauslist), mean(thetacausalcoverage), mean(thetaehwcoverage),
           mean(olscausalcoverage), mean(olsehwcoverage)))}

result1 <- uncertaintysim(1000, rho=.5, k = 1, psi = c(2),iter = 10000, nu = 1, delta =1)

result2 <- uncertaintysim(1000, rho=.5, k = 2, psi = c(2, 2),iter = 10000, nu = 1, delta =1)

result3 <- uncertaintysim(1000, rho=.5, k = 3, psi = c(2,2,2),iter = 10000, nu = 1, delta =1)

result4 <- uncertaintysim(1000, rho=.5, k = 4, psi = c(2,2,2,2),iter = 10000, nu = 1, delta =1)

result5 <- uncertaintysim(1000, rho=.5, k = 10, psi = c(2,2,2,2,2,2,2,2,2,2),iter = 10000, nu = 1, delta =1)

result6 <- uncertaintysim(10000, rho=.5, k = 3, psi = c(2,2,2),iter = 10000, nu = 1, delta =1)

simput <- c(result1, result2, result3, result4,result5, result6)

brew(file = "Text Input/6simformat.brew", output = "Text Output/Table 1.txt")

uncertaintysim(1000, rho=.5, k = 1, psi = c(2),iter = 1000, nu = 1, delta =.01)

uncertaintysim(1000, rho=.5, k = 1, psi = c(2),iter = 1000, nu = 1, delta =.1)

uncertaintysim(1000, rho=.5, k = 1, psi = c(2),iter = 1000, nu = 1, delta =.25)

uncertaintysim(1000, rho=.5, k = 1, psi = c(2),iter = 1000, nu = 1, delta =.5)

uncertaintysim(1000, rho=.5, k = 1, psi = c(2),iter = 1000, nu = 1, delta =.75)

uncertaintysim(1000, rho=.5, k = 1, psi = c(2),iter = 1000, nu = 1, delta =1)

uncertaintysim(1000, rho=.5, k = 1, psi = c(2),iter = 1000, nu = 1, delta =10)

uncertaintysim(1000, rho=.5, k = 1, psi = c(0),iter = 1000, nu = 1, delta =1)

uncertaintysim(1000, rho=.5, k = 1, psi = c(1),iter = 1000, nu = 1, delta =1)

uncertaintysim(1000, rho=.5, k = 1, psi = c(5),iter = 1000, nu = 1, delta =1)

uncertaintysim(1000, rho=.5, k = 1, psi = c(10),iter = 1000, nu = 1, delta =1)

uncertaintysim(1000, rho=.5, k = 1, psi = c(100),iter = 1000, nu = 1, delta =1)

uncertaintysim(1000, rho=.1, k = 1, psi = c(2),iter = 1000, nu = 1, delta =1)

uncertaintysim(1000, rho=.25, k = 1, psi = c(2),iter = 1000, nu = 1, delta =1)

uncertaintysim(1000, rho=.5, k = 1, psi = c(2),iter = 1000, nu = 1, delta =1)

uncertaintysim(1000, rho=.75, k = 1, psi = c(2),iter = 1000, nu = 1, delta =1)

uncertaintysim(1000, rho=1, k = 1, psi = c(2),iter = 1000, nu = 1, delta =1)

uncertaintysim(1000, rho=.1, k = 1, psi = c(2),iter = 1000, nu = 1, delta =1, sigma2 = .5)

uncertaintysim(1000, rho=.1, k = 1, psi = c(2),iter = 1000, nu = 1, delta =1, sigma2 = 1)

uncertaintysim(1000, rho=.1, k = 1, psi = c(2),iter = 1000, nu = 1, delta =1, sigma2 = 2)

uncertaintysim(1000, rho=.1, k = 1, psi = c(2),iter = 1000, nu = 1, delta =1, sigma2 = 5)

uncertaintysim(1000, rho=.1, k = 1, psi = c(2),iter = 1000, nu = 1, delta =1, sigma2 = 10)

uncertaintysim(1000, rho=.1, k = 1, psi = c(2),iter = 1000, nu = 1, delta =1, sigma2 = 50)

brew(file = "Text Input/brewtest.brew", output = "Text Output/test.txt")

sink()