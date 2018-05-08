



# rm(list=ls())
# gc()
#library(invgamma)

#library(MCMCpack)
library(nimble)

#' A Cat Function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' cat_function()

bsvEst=function (x){
  y=x
  N = 100
  n=N/2
  vol = rep(0, length(y))
  
  aa = matrix(rep(0, N*3),  nrow=N, ncol=3)
  
  #  number of iterations
  T1 = length(y)
  #  number of observations
  z = 1:T1
  
  # the initial values for the start of the Gibbs
  
  phi = 0.7
  sigma = 1.5
  mu = 0.8
  
  alpha_mu = 0.0
  sig_mu = 10.0
  alpha_phi = 0.0
  sig_phi = 10.0
  alpha_psi = 0.0
  sig_psi = 10.0
  alpha_tau = 2
  sig_tau = 0.0025
  
  
  muEst=0
  phiEst=0
  sigmaEst=0
  
  d1 <- data.frame(1, mu, phi,  sigma)                
  write.table(d1, "bsv-mcmc-R.txt", row.names = FALSE, col.names = FALSE)
  
  norm_temp = rnorm(1,0,1)
  z[1] = mu + sigma/(1.0-phi*phi)^0.5 * norm_temp
  
  for (j in 2:T1)
  {
    norm_temp = rnorm(1, 0, 1)
    z[j] = mu + phi * (z[j-1]-mu) + sigma*norm_temp
  }
  
  for (i in 1: N)
  {  
    #print(i)
    u1 = mu + phi*(z[2] - mu) - sigma*sigma/2.0
    rnd = rnorm(1, u1, sigma)
    
    a1 = exp(-y[1] * y[1] / 2.0 * exp(-rnd) + y[1] * y[1] / 2.0 * exp(-z[1]))
    
    u = runif(1, min = 0, max = 1)
    if (u <= min(1.0, a1))
    {
      z[1] = rnd
    }
    for (j in 2: (T1-1))
    {
      sigma_temp = sigma*sigma/(1.0 + phi*phi)
      s1 = sigma_temp ^ 0.5
      u1 = mu + phi*(z[j-1] - mu + z[j+1] - mu)/(1.0 + phi*phi) - sigma_temp / 2.0
      rnd = rnorm(1, u1, s1)
      a1 = exp(-y[j] * y[j] / 2.0 * exp(-rnd) + y[j] * y[j] / 2.0 * exp(-z[j]))
      
      u = runif(1, min = 0, max = 1)
      if (u <= min(1.0, a1))
      {
        z[j] = rnd
      }
    }
    
    u1 = mu + phi * (z[T1-1] - mu) -sigma*sigma/2
    rnd = rnorm(1, u1, sigma)
    a1 = exp(-y[T1] * y[T1] / 2.0 * exp(-rnd) + y[T1] * y[T1] / 2.0 * exp(-z[T1]))
    
    u = runif(1, 0,1)
    if (u <= min(1.0, a1))
    {
      z[T1] = rnd
    }
    # generate the random number for the variable mu
    
    a1 = (T1 - 1) * (1.0 - phi)^2.0 + (1.0 - phi * phi)
    a1 = a1 / sigma/sigma
    # a1 = a1 + 1.0 / sig_mu**2.0
    
    b1 = 0.0
    for (h in 1: (T1-1))
    {
      b1 = b1 + (1.0 - phi) * (z[h + 1] - phi * z[h])
    }
    
    b1 = b1+(1.0-phi*phi)*z[1]
    b1 = b1 / sigma/sigma
    # b1 = b1 + alpha_mu / sig_mu / sig_mu
    
    temp1 = rnorm(1, 0, 1)
    mu = b1 / a1 + 1.0 / a1^0.5 * temp1
    
    ######################
    
    # u=runif(1, 0,1)
    # 
    # latent=sqrt(1.0-phi*phi)*u
    # bounder=min(1.0, sqrt(1.0-latent*latent))
    # 
    a1 = 0.0
    for (h in 1: (T1-1))
    {
      a1 = a1 + (z[h] - mu)^2
    }
    
    a1 = a1 - (z[1] - mu)^2
    a1 = a1 / sigma^2
    # a1 = a1 + 1.0 / sig_phi / sig_phi
    
    b1 = 0.0
    for ( h in 1: (T1-1))
    {
      b1 = b1 + (z[h] - mu) * (z[h + 1] - mu)
    }
    b1 = b1 / sigma/sigma
    # b1 = b1 + alpha_phi / sig_phi / sig_phi
    
    mu1 = b1 / a1
    sig =sqrt(1.0 / a1)
    
    
    # u=runif(1, 0,1)
    # 
    # u=u*exp(-(phi-mu1)^2/2.0/sig^2);
    # 
    # c1_left= mu1-sqrt(-2.0*sig^2*log(u))
    # 
    # 
    # c1_right=mu1+sqrt(-2.0*sig^2*log(u))
    # 
    # left=max(-bounder,c1_left);
    # 
    # 
    # right=min(bounder, c1_right);
    # 
    # u=runif(1, 0,1)
    # phi=left + (right-left)*u
    # 
    # 
    rnd = rnorm(1, mu1, sig)
    mm=1
    while(rnd<=0 | rnd>=1){
      rnd = rnorm(1, mu1, sig)
      mm=mm+1
      if (mm == 10) break; 
    }
    
    if(mm !=10){
      phi=rnd  }
    
    # 
    # a11 = sqrt((1-rnd*rnd)/(1-phi*phi))
    # 
    # u = runif(1, 0, 1)
    # if (u <= min(1.0, a11))
    # {
    #   phi = rnd
    # }
    # 
    
    # p = dnorm((-bounder - mu1) / sig) + (dnorm((bounder - mu1) / sig) - dnorm((-bounder - mu1) / sig)) * u
    #   phi = mu1 + sig * qnorm(p)
    
    # generate the random  numbers for sigma^2
    
    a1 = 0.0
    for (h in 1: (T1-1))
    {
      a1 = a1 + (z[h + 1] - mu - phi * (z[h] - mu))^2.0
    }
    
    a1 = a1 + (z[1] - mu)^2.0 * (1.0 - phi * phi)
    a1 = a1 / 2.0
    # a1 = a1 + sig_tau
    
    # sigma =1/ (rinvgamma(1, (T1 / 2-1), a1))^0.5
    sigma= rinvgamma(1,(T1 / 2-1), a1) 
    
    sigma=sqrt(sigma)
    #sigma=0.2
    #print(c(mu, phi, sigma))
    
    if(i >n){
      muEst= muEst + mu 
      phiEst=phiEst +phi
      sigmaEst = sigmaEst +sigma
      vol=vol +z
    }
    
    # aa[i,1] = mu
    # aa[i,2]=phi
    # aa[i,3]=sigma
    # 
  }
  
  # df=data.frame(aa)
  # names(df) <- c("mu", "phi", "sigma")
  # write.table(df, "bsv-mcmc-R.txt", row.names = FALSE,
  #             col.names = FALSE, append = TRUE)
  
  muEst= muEst/(N-n) 
  phiEst=phiEst/(N-n) 
  sigmaEst = sigmaEst/(N-n) 
  vol=vol/(N-n)
  
  coef =c(muEst, phiEst, sigmaEst)
  
  # resi=y*exp(-vol/2)
  
  list(coefficients = coef,
       vol=vol
  )
  
  
}


###############################################
#' A Cat Function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' cat_function()

bsvSimu=function(mu1, phi1, sigma1, length1){
  
  
  mu=mu1
  phi=phi1
  sigma=sigma1
  N = length1
  
  
  print(mu)
  
  print(phi)
  print(sigma)
  
  print(N)
  
  n = 20000
  z = 1:N
  x = z
  y=z
  
  
  z[1] = rnorm(1, mu, sigma/sqrt(1-phi*phi))
  for (i in 2:n)
    
  {
    z[2] = mu + phi*(z[1] -mu) + sigma*rnorm(1, 0,1)
    y[1] = exp(z[2]/2) * rnorm(1,0,1)
    z[1] = z[2]
    
  }  
  
  for (i in 1:N)
    
  {
    z[2] = mu + phi*(z[1] -mu) + sigma*rnorm(1, 0,1)
    y[i] = exp(z[2]/2) * rnorm(1,0,1)
    x[i] = z[2]
    z[1] = z[2]
    
  }  
  
  list(ySimu = y,
       volSimu=x)
}

# simu=bsvSimu(-5.5, 0.95, 0.15, 1200)
# 
# simu=bsvSimu(-0.68, 0.96, 0.19, 1500)
# 
# 
# print(ks.test(simu$ySimu*exp(-simu$volSimu /2), "pnorm", 0, 1))
# 
# #plot(simu$volSimu, fit$vol)
# 
# fit=bsvEst(simu$ySimu)
# 
# print(fit$coefficients)
# 
# print(ks.test(simu$ySimu*exp(-fit$vol/2), "pnorm", 0, 1))

 
# mydata = read.table("simulated-returns.txt") 

# mydata = read.table("aus2005.txt") 
# y=mydata$V1
# fit=bsvEst(y)
# fit$coefficients 
# 
# print(fit$coefficients)
# print(ks.test(y*exp(-fit$vol/2), "pnorm", 0, 1))
# 
# zz=y*exp(-fit$vol/2)
# par(mar=c(2,2,2,2))
# par(mfrow=c(1,1))
# qqnorm(zz, xlim=c(-3,3), ylim=c(-3,3))

# 
# 
# nn=1:length(simu$ySimu)
# 
# par(mar=c(1,1,1,1))
# par(mfrow=c(1,1))
# 
# plot(nn,abs(simu$ySimu), type="l", col="blue")
# 
# lines(nn, simu$volSimu, type="l", col="red")
# 
# 
# zz=simu$ySimu*exp(-fit$vol/2)
# 
# par(mar=c(2,2,2,2))
# par(mfrow=c(1,1))
# qqnorm(zz, xlim=c(-3,3), ylim=c(-3,3))
# 
# print(ks.test((simu$ySimu)*exp(-fit$vol/2), "pnorm", 0, 1))
# 
# print(ks.test(zz, "pnorm", 0, 1))
# 
# par(mar=c(1,1,1,1))
# par(mfrow=c(2,1))
# 
# plot(nn,simu$ySimu, type="l", col="blue", ylim=c(0, max(simu$ySimu)) )
# 
# plot(nn, exp(fit$vol/2), type="l", col="red", ylim=c(min(exp(fit$vol/2)), max(exp(fit$vol/2))))
# 





