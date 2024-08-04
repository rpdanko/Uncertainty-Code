setwd("C:/Users/chris/Documents/Uncertainty Code")

library(MASS)
library(AER)
library(brew)
library(simstudy)
library(data.table)
set.seed(48824)

uncertaintysim <- function(samples, rho, k, psi, iter, sigma2=1, g, gsig){
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
    xi <- c()
  for (i in 1:g){
    
  }
  xi = rnorm(samples, 0,1)
  #generate w; note that in this setup, w =u  
  w = rnorm(samples, 0, 1)
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

###playing around

check <- matrix(1, nrow= 1000, ncol =1000) 

check[col(check) != row(check)] <- .5

n =3

mu <- c()
mu[1:n] <- 0
final <- mvrnorm(1, mu = mu, Sigma = check)



final1 <- final[1:500]
final2 <- final[501:1000]
cor(final1, final2)
cov(final1, final2)

t(final1) %*% final2 /n

## try something else

Sigma = matrix(c(1,.5,.5, 1), ncol = 2)
R = chol(Sigma)
n =1000
X = t(R) %*%matrix(rnorm(n*2), 2)

X %*% t(X)/n


R = chol(check)

### Continuing to mess around


defC <- defData(varname = "xi", formula = "0", dist= "normal")
dc <- genData(n = 1000, dtDefs = defC, id = "unit")
dc <- genCluster(dtClust = dc, cLevelVar = "unit", numIndsVar = 1000)