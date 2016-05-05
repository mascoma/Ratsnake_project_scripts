#misfits for one tree
###The functions for the coalescent models were written by Hélène Morlon and publshed in the PLoSB (2010). I have modified them to some degree, but the original descriptions as she wrote them appear below:

library(ape)

library(picante)

library(laser)

library(geiger)

library(qpcR)

#This code computes the likelihood of a given phylogeny under the constant diversity/exponentially varying turnover rate model (Model 2 from the PloSB 2010 paper), with turnover rate at present tau0, exponential variation in turnover rate gamma, and diversity N0



getLikelihood.coalMoranEXP<-function(Vtimes,Ntips,tau0,gamma,N0)
  
  
  
{
  
  Ttimes <- diff(Vtimes)
  
  Vtimes<-Vtimes[2:length(Vtimes)]
  
  nbint<-length(Ttimes)
  
  samp<-seq((Ntips-2),(Ntips-nbint-1),by=-1)
  
  indLikelihood<-samp*(samp+1)/2*2*tau0/N0*exp(gamma*Vtimes)*exp(-samp*(samp+1)/2*2*tau0/N0*1/gamma*exp(gamma*Vtimes)*(1-exp(-gamma*Ttimes)))
  
  res<-sum(log(indLikelihood))
  
  return(list("res"=res,"all"=indLikelihood))
  
}



#This code computes the likelihood of a given phylogeny under the constant diversity/constant turnover rate model (Model 1 from the PloSB 2010 paper), with turnover rate tau0 and diversity N0



getLikelihood.coalMoranCST<-function(Vtimes,Ntips,tau0,N0)
  
{
  
  Ttimes <- diff(Vtimes)
  
  nbint<-length(Ttimes)
  
  samp<-seq((Ntips-2),(Ntips-nbint-1),by=-1)
  
  indLikelihood<-samp*(samp+1)/2*2*tau0/N0*exp(-samp*(samp+1)/2*2*tau0/N0*Ttimes)
  
  res<-sum(log(indLikelihood))
  
  return(list("res"=res,"all"=indLikelihood))
  
}

#This code computes the likelihood of a given phylogeny under various flavors of the birth-death model (Models 2 to 6 from the PloSB 2010 paper), with parameters lamb0,alpha,mu0 and beta, and diversity N0

getLikelihood.coalBD <- function(Vtimes,Ntips,lamb0,alpha,mu0,beta,N0,pos=TRUE)
  
  # The extinction rate is forced to be less than the speciation rate over the history of the clade
  
  
  
{
  
  
  Ttimes <- diff(Vtimes)
  
  Vtimes <- Vtimes[2:length(Vtimes)]
  
  nbint<-length(Ttimes)
  
  samp<-seq((Ntips-2),(Ntips-nbint-1),by=-1)
  
  
  times<-c(0,sort(Vtimes))
  
  
  if (min(abs(lamb0)*exp(alpha*times)-abs(mu0)*exp(beta*times))<=0)
    
  {
    
    indLikelihood<-0*vector(length=length(samp))
    
    res<-sum(log(indLikelihood))
    
  }
  
  
  else if ((alpha==0) & (beta==0))
    
  {  if (pos==FALSE)
    
  {r<-lamb0-mu0
   
   indLikelihood<-samp*(samp+1)/2*2*lamb0/N0*1/(exp(-r*Vtimes))*exp(-samp*(samp+1)/2*2*lamb0/N0/r*exp(r*Vtimes)*(1-exp(-r*Ttimes)))}
  
  else
    
  {r<-abs(abs(lamb0)-abs(mu0))
   
   indLikelihood<-samp*(samp+1)/2*2*abs(lamb0)/N0*1/(exp(-r*Vtimes))*exp(-samp*(samp+1)/2*2*abs(lamb0)/N0/r*exp(r*Vtimes)*(1-exp(-r*Ttimes)))}
  
  res<-sum(log(indLikelihood))
  
  }
  
  
  else
    
    
  {
    
    if ((beta==0) & !(alpha==0))
      
    {	if (pos==FALSE)
      
    {demfun<-function(x){2*lamb0*exp(alpha*x)/(N0*exp(lamb0/alpha*(1-exp(alpha*x))+mu0*x))}}
    
    else {demfun<-function(x){2*abs(lamb0)*exp(alpha*x)/(N0*exp(abs(lamb0)/alpha*(1-exp(alpha*x))+abs(mu0)*x))}}
    
    }
    
    
    else if ((alpha==0) & !(beta==0))
      
    {	if (pos==FALSE)
      
    {demfun<-function(x){2*lamb0/(N0*exp(-lamb0*x-mu0/beta*(1-exp(beta*x))))}}
    
    else
      
    {demfun<-function(x){2*abs(lamb0)/(N0*exp(-abs(lamb0)*x-abs(mu0)/beta*(1-exp(beta*x))))}}
    
    }
    
    
    else {	if (pos==FALSE)
      
    {demfun<-function(x){2*lamb0*exp(alpha*x)/(N0*exp(lamb0/alpha*(1-exp(alpha*x))-mu0/beta*(1-exp(beta*x))))}}
    
    else {demfun<-function(x){2*abs(lamb0)*exp(alpha*x)/(N0*exp(abs(lamb0)/alpha*(1-exp(alpha*x))-abs(mu0)/beta*(1-exp(beta*x))))}}
    
    }
    
    
    
    if (FALSE %in% is.finite(demfun(Vtimes)))
      
    {
      
      indLikelihood<-0*vector(length=length(samp))
      
      res<-sum(log(indLikelihood))}
    
    
    else
      
    {
      
      integrals<-c()
      
      demfunval<-c()
      
      
      for (i in 1:length(Vtimes))
        
      {
        
        demfunvali<-demfun(Vtimes[i])
        
        integrali<-integrate(demfun,(Vtimes[i]-Ttimes[i]),Vtimes[i],stop.on.error=FALSE)$value
        
        
        demfunval<-c(demfunval,demfunvali)
        
        integrals<-c(integrals,integrali)}
      
      
      indLikelihood<-samp*(samp+1)/2*demfunval*exp(-samp*(samp+1)/2*integrals)
      
      res<-sum(log(indLikelihood))
      
    }
    
  }
  
  
  return(list("res"=res,"all"=indLikelihood))}



#This code fits the constant diversity/exponentially varying turnover rate model (Model 2 from the PloSB 2010 paper) to a given phylogeny, by maximum likelihood, using the Nelder-Mead algorithm

#Outputs are the log-likelihood, the second order Akaike's Information Criterion, and the maximum likelihood estimates of the turnover rate at present (tau0) and the exponential variation in turnover rate (gamma). See notations in the PloSB 2010 paper.

#The code uses the ape package

#The code uses the getLikelihood.coalMoranEXP code. use source("getLikelihood.coalMoranEXP.r")





fitcoalMoranEXP<-function (phylo, tau0=10^-2, gamma=1, meth = "Nelder-Mead", N0=0, Vtimes=FALSE)
  
  #Assuming we know N0 (the total number of species at present), we estimate tau0 and gamma. If N0=0 (default), N0 is set to the #number of tips in the phylogeny (i.e. the phylogeny is assumed to be 100% complete). Otherwise, enter the number of species at #present.
  
  #the input in tau0 and gamma are initial parameter values. Try several inital values to make sure you are not stuck in a local optimum.
  
  #The code can take as input either a phylogeny (default, Vtimes=FALSE), or branching times (if Vtimes=TRUE)
  
  
  
  
  
{
  
  if (Vtimes==TRUE) {
    
    Vtimes<-sort(phylo)
    
    Ntips<-length(phylo)+1
    
    if (N0==0) {N0<-Ntips}
    
  }
  
  
  else{
    
    Vtimes <- sort(branching.times(phylo))
    
    Ntips<-Ntip(phylo)
    
    if (N0==0) {N0<-Ntips}
    
  }
  
  
  
  init<-c(tau0,gamma)
  
  nbpar<-length(init)
  
  nbobs<-length(Vtimes)-1
  
  
  
  optimLH.MoranEXP <- function(init) {
    
    tau0 <- init[1]
    
    gamma <- init[2]
    
    LH <- getLikelihood.coalMoranEXP(Vtimes,Ntips,tau0,gamma,N0)$res
    
    return(-LH)
    
  }
  
  
  
  temp<-suppressWarnings(optim(init, optimLH.MoranEXP, method = meth))
  
  
  
  res <- list(model = "MoranEXP", LH = -temp$value, aicc = 2 *
                
                temp$value + 2*nbpar + 2*nbpar*(nbpar+1)/(nbobs-nbpar-1), tau0 = temp$par[1], gamma=temp$par[2])
  
  return(res)
  
}





#This code fits the constant diversity/constant turnover rate model (Model 1 from the PloSB 2010 paper) to a given phylogeny, by maximum likelihood, using the Nelder-Mead algorithm

#Outputs are the log-likelihood, the second order Akaike's Information Criterion, and the maximum likelihood estimate of the turnover rate (tau0). See notations in the PloSB 2010 paper.

#The code uses the ape package

#The code uses the getLikelihood.coalMoranCST code. use source("getLikelihood.coalMoranCST.r")



fitcoalMoranCST<-function (phylo, tau0=1, meth = "Nelder-Mead", N0=0, Vtimes=FALSE)
  
  #Assuming we know N0 (the total number of species at present), we estimate tau0. If N0=0 (default), N0 is set to the number of tips in the phylogeny (i.e. the phylogeny is assumed to be 100% complete). Otherwise, enter the number of species at present.
  
  #the input in tau0 is an initial value for the turnover rate. Try several inital values to make sure you are not stuck in a local optimum.
  
  #The code can take as input either a phylogeny (default, Vtimes=FALSE), or branching times (if Vtimes=TRUE)
  
  
  
{
  
  if (Vtimes==TRUE) {
    
    Vtimes<-sort(phylo)
    
    Ntips<-length(phylo)+1
    
    if (N0==0) {N0<-Ntips}
    
  }
  
  
  else{
    
    Vtimes <- sort(branching.times(phylo))
    
    Ntips<-Ntip(phylo)
    
    if (N0==0) {N0<-Ntips}
    
  }
  
  
  init<-c(tau0)
  
  
  
  nbpar<-length(init)
  
  nbobs<-length(Vtimes)-1
  
  
  
  optimLH.MoranCST <- function(init) {
    
    tau0 <- init[1]
    
    LH <- getLikelihood.coalMoranCST(Vtimes,Ntips,tau0,N0)$res
    
    return(-LH)
    
  }
  
  
  
  temp<-suppressWarnings(optim(init, optimLH.MoranCST, method = meth))
  
  
  
  res <- list(model = "MoranCST", LH = -temp$value, aicc = 2 *
                
                temp$value + 2*nbpar + 2*nbpar*(nbpar+1)/(nbobs-nbpar-1), tau0 = temp$par[1])
  
  
  
  return(res)
  
}

#This code fits various flavors of the birth-death model (Models 2 to 6 from the PloSB 2010 paper) to a given phylogeny, by maximum #likelihood, using the Nelder-Mead algorithm

#Outputs are the log-likelihood, the second order Akaike's Information Criterion, and the maximum likelihood estimates of the parameters of diversification. Depending on the model, these parameters include a combination of the speciation rate at present #(lamb0), the exponential variation in speciation rate (alpha), the extinction rate at present (mu0), the exponential variation in extinction rate (beta) and the extinction fraction (extinction rate/speciation rate, eps). See notations in the PloSB 2010 paper.

#The code uses the ape package

#The code uses the getLikelihood.coalBD code. use source("getLikelihood.coalBD.r")



fitcoalBD<-function (phylo,lamb0=0.1,alpha=1,mu0=0.01,beta=0,meth = "Nelder-Mead",N0=0,cst.lamb=FALSE,cst.mu=FALSE,fix.eps=FALSE,mu.0=FALSE,pos=TRUE,Vtimes=FALSE)
  
  
  
  #The default settings allow to fit the most general model where the rates of speciation and extinction vary over time without a fixed extinction fraction (Model 4d from the PloSB 2010 paper). cst.lamb=TRUE forces the speciation rate to be constant over time (used to #fit Models 3, 5 and 4b). cst.mu=TRUE forces the extinction rate to be constant over time (used to fit Models 3, and 4a). fix.eps forces the extinction fraction to be constant over time (used to fit Model 4c). mu.0=TRUE forces the extinction rate to 0 (used to fit Models 5 #and 6).
  
  #pos=TRUE (the default) forces the rates of speciation and extinction to be positive. pos=FALSE removes this forcing.
  
  
  
{
  
  if (Vtimes==TRUE) {
    
    Vtimes<-sort(phylo)
    
    Ntips<-length(phylo)+1
    
    if (N0==0) {N0<-Ntips}
    
  }
  
  
  else{
    
    Vtimes <- sort(branching.times(phylo))
    
    Ntips<-Ntip(phylo)
    
    if (N0==0) {N0<-Ntips}
    
  }
  
  
  nbobs<-length(Vtimes)-1
  
  
  
  
  #pure birth	constant rates (Model 5)
  
  if (mu.0==TRUE & cst.lamb==TRUE)
    
  {init<-c(lamb0)}
  
  
  #birth-death constant rates	(Model 3)
  
  else if (cst.mu==TRUE & cst.lamb==TRUE)
    
  {init <- c(lamb0,mu0)}
  
  
  #pure birth varying speciation rate (Model 6)
  
  else if (mu.0==TRUE & cst.lamb==FALSE)
    
  {init<-c(lamb0,alpha)}
  
  
  #birth-death varying speciation rate (Model 4a)
  
  else if (cst.mu==TRUE & cst.lamb==FALSE)
    
  {init <- c(lamb0,alpha,mu0)}
  
  
  
  #birth-death varying extinction rate (Model 4b)
  
  else if (cst.mu==FALSE & cst.lamb==TRUE)
    
  {init <- c(lamb0,mu0,beta)}
  
  
  
  #birth-death varying speciation rate and constant extinction fraction (Model 4c)
  
  else if (fix.eps==TRUE)
    
  {init <- c(lamb0,alpha,mu0/lamb0)}
  
  
  
  #birth-death varying speciation and extinction rates (Model 4d)
  
  else
    
  {init = c(lamb0,alpha,mu0,beta)}
  
  
  
  nbpar<-length(init)
  
  
  
  ############################################################
  
  
  
  if (mu.0==TRUE & cst.lamb==TRUE)
    
  {optimLH.coalBD <- function(init) {
    
    lamb0 <- init[1]
    
    LH <- getLikelihood.coalBD(Vtimes,Ntips,lamb0,alpha=0,mu0=0,beta=0,N0,pos=pos)$res
    
    return(-LH)}
   
  }
  
  
  
  else if (cst.mu==TRUE & cst.lamb==TRUE)
    
  {optimLH.coalBD <- function(init) {
    
    lamb0 <- init[1]
    
    mu0 <- init[2]
    
    LH <- getLikelihood.coalBD(Vtimes,Ntips,lamb0,alpha=0,mu0,beta=0,N0,pos=pos)$res
    
    return(-LH)}
   
  }
  
  
  
  else if (mu.0==TRUE & cst.lamb==FALSE)
    
  {optimLH.coalBD <- function(init) {
    
    lamb0 <- init[1]
    
    alpha <- init[2]
    
    LH <- getLikelihood.coalBD(Vtimes,Ntips,lamb0,alpha,mu0=0,beta=0,N0,pos=pos)$res
    
    return(-LH)}
   
  }
  
  
  
  else if (cst.mu==TRUE & cst.lamb==FALSE)
    
  {optimLH.coalBD <- function(init) {
    
    lamb0 <- init[1]
    
    alpha <- init[2]
    
    mu0 <- init[3]
    
    LH <- getLikelihood.coalBD(Vtimes,Ntips,lamb0,alpha,mu0,beta=0,N0,pos=pos)$res
    
    return(-LH)}
   
  }
  
  
  else if (cst.mu==FALSE & cst.lamb==TRUE)
    
  {optimLH.coalBD <- function(init) {
    
    lamb0 <- init[1]
    
    mu0 <- init[2]
    
    beta<- init[3]
    
    LH <- getLikelihood.coalBD(Vtimes,Ntips,lamb0,alpha=0,mu0,beta,N0,pos=pos)$res
    
    return(-LH)}
   
  }
  
  
  
  else if (fix.eps==TRUE)
    
  {optimLH.coalBD <- function(init) {
    
    lamb0 <- init[1]
    
    alpha <- init[2]
    
    eps <- init[3]
    
    LH <- getLikelihood.coalBD(Vtimes,Ntips,lamb0,alpha=0,mu0,beta,N0,pos=pos)$res
    
    return(-LH)}
   
  }
  
  
  
  else
    
    
    
  {optimLH.coalBD <- function(init) {
    
    lamb0 <- init[1]
    
    alpha <- init[2]
    
    mu0 <- init[3]
    
    beta <- init[4]
    
    LH <- getLikelihood.coalBD(Vtimes,Ntips,lamb0,alpha,mu0,beta,N0,pos=pos)$res
    
    return(-LH)}
   
  }
  
  
  
  #######################################################################################   
  
  
  
  temp <- optim(init, optimLH.coalBD, method = meth,control=list(ndeps=10^(-4)))
  
  
  
  if (mu.0==TRUE & cst.lamb==TRUE)
    
  {
    
    if (pos==FALSE)
      
    {res <- list(model = "Pure birth constant speciation", LH = -temp$value, aicc = 2 *
                   
                   temp$value + 2*nbpar + 2*nbpar*(nbpar+1)/(nbobs-nbpar-1), lamb0 = temp$par[1])}
    
    else
      
    {res <- list(model = "Pure birth constant speciation", LH = -temp$value, aicc = 2 *
                   
                   temp$value + 2*nbpar + 2*nbpar*(nbpar+1)/(nbobs-nbpar-1), lamb0 = abs(temp$par[1]))}}
  
  
  
  else if (cst.mu==TRUE & cst.lamb==TRUE)
    
  {
    
    if (pos==FALSE)
      
    {res <- list(model = "Birth-death constant rates", LH = -temp$value, aicc = 2 *
                   
                   temp$value + 2*nbpar + 2*nbpar*(nbpar+1)/(nbobs-nbpar-1), lamb0 = temp$par[1],mu0 = temp$par[2])}
    
    else
      
    {res <- list(model = "Birth-death constant rates", LH = -temp$value, aicc = 2 *
                   
                   temp$value + 2*nbpar + 2*nbpar*(nbpar+1)/(nbobs-nbpar-1), lamb0 = abs(temp$par[1]),mu0 = abs(temp$par[2]))}}
  
  
  
  else if (mu.0==TRUE & cst.lamb==FALSE)
    
  {
    
    if (pos==FALSE)
      
    {res <- list(model = "Pure birth varying speciation", LH = -temp$value, aicc = 2 *
                   
                   temp$value + 2*nbpar + 2*nbpar*(nbpar+1)/(nbobs-nbpar-1), lamb0 = temp$par[1], alpha = temp$par[2])}
    
    else
      
    {res <- list(model = "Pure birth varying speciation", LH = -temp$value, aicc = 2 *
                   
                   temp$value + 2*nbpar + 2*nbpar*(nbpar+1)/(nbobs-nbpar-1), lamb0 = abs(temp$par[1]), alpha = temp$par[2])}}
  
  
  
  else if (cst.mu==TRUE & cst.lamb==FALSE)
    
  {
    
    if (pos==FALSE)
      
    {res <- list(model = "Birth-death varying speciation constant extinction", LH = -temp$value, aicc = 2 *
                   
                   temp$value + 2*nbpar + 2*nbpar*(nbpar+1)/(nbobs-nbpar-1), lamb0 = temp$par[1], alpha = temp$par[2], mu0 = temp$par[3])}
    
    else
      
    {res <- list(model = "Birth-death varying speciation constant extinction", LH = -temp$value, aicc = 2 *
                   
                   temp$value + 2*nbpar + 2*nbpar*(nbpar+1)/(nbobs-nbpar-1), lamb0 = abs(temp$par[1]), alpha = temp$par[2], mu0 = abs(temp$par[3]))}}
  
  
  
  else if (cst.mu==FALSE & cst.lamb==TRUE)
    
  {
    
    if (pos==FALSE)
      
    {res <- list(model = "Birth-death constant speciation varying extinction", LH = -temp$value, aicc = 2 *
                   
                   temp$value + 2*nbpar + 2*nbpar*(nbpar+1)/(nbobs-nbpar-1), lamb0 = temp$par[1], mu0 = temp$par[2], beta = temp$par[3])}
    
    else
      
    {res <- list(model = "Birth-death constant speciation varying extinction", LH = -temp$value, aicc = 2 *
                   
                   temp$value + 2*nbpar + 2*nbpar*(nbpar+1)/(nbobs-nbpar-1), lamb0 = abs(temp$par[1]), mu0 = abs(temp$par[2]), beta = temp$par[3])}}
  
  
  
  else if (fix.eps==TRUE)
    
  {
    
    if (pos==FALSE)
      
    {res <- list(model = "Birth-death constant extinction fraction", LH = -temp$value, aicc = 2 *
                   
                   temp$value + 2*nbpar + 2*nbpar*(nbpar+1)/(nbobs-nbpar-1), lamb0 = temp$par[1], alpha = temp$par[2], eps = temp$par[3])}
    
    else
      
    {res <- list(model = "Birth-death constant extinction fraction", LH = -temp$value, aicc = 2 *
                   
                   temp$value + 2*nbpar + 2*nbpar*(nbpar+1)/(nbobs-nbpar-1), lamb0 = abs(temp$par[1]), alpha = temp$par[2], eps = abs(temp$par[3]))}}
  
  
  
  else
    
  {
    
    if (pos==FALSE)
      
    {res <- list(model = "Birth-death varying speciation and extinction", LH = -temp$value, aicc = 2 *
                   
                   temp$value + 2*nbpar + 2*nbpar*(nbpar+1)/(nbobs-nbpar-1), lamb0 = temp$par[1], alpha = temp$par[2], mu0 = temp$par[3],beta = temp$par[4])}
    
    else
      
    {res <- list(model = "Birth-death varying speciation and extinction", LH = -temp$value, aicc = 2 *
                   
                   temp$value + 2*nbpar + 2*nbpar*(nbpar+1)/(nbobs-nbpar-1), lamb0 = abs(temp$par[1]), alpha = temp$par[2], mu0 = abs(temp$par[3]),beta = temp$par[4])}}
  
  
  
  
  
  return(res)
  
}



#Monotonic decay of speciation rate, sensu Rabosky and Lovette 2008       


fit.linear <- function(phy)
{
  N <- length(phy$tip.label)
  x <- c(0, branching.times(phy))
  T <- max(x)
  t <- T - x
  lambda.t <- function(lambda, mT) lambda * (1 - t[3:N] / mT)
  chi.t <- function(lambda, mT, time) exp(-((lambda * (time - T) * (-2 * mT + time + T)) / (2 * mT)))
  linear.lik <- function(lambda, mT) lfactorial(N - 1) + sum(log(lambda.t(lambda, mT))) + sum(log(chi.t(lambda, mT, t[3:N]))) + 2 * log(chi.t(lambda, mT, t[2]))
  lnL <- function(p)
  {
    lambda <- exp(p[1])
    mT <- exp(p[2]) + T
    -linear.lik(lambda, mT)
  }
  res <- NULL
  temp <- optim(c(log(log(N) / T), log(T)), lnL)
  res$LH <- -temp$value
  res$lambda0 <- exp(temp$par[1])
  res$mT <- exp(temp$par[2]) + T
  res$AIC <- -2 * res$LH + 4
  res
}


#Yule model with incomplete sampling

fit.yule.f <- function(phy, sampling.f = 1)
  
{
  
  N <- length(phy$tip.label)
  
  x <- c(0, branching.times(phy))
  
  tau <- x[2]
  
  t <- tau - x
  
  ptt.f <- function(lambda, time, tau) sampling.f / (sampling.f - (sampling.f - 1) * exp(lambda * (time - tau)))
  
  beta <- function(lambda, time) 1 - ptt.f(lambda, 0, time) * exp(-lambda * time)
  
  yule.f.lik <- function(lambda) lfactorial(N - 1) + (N - 2) * log(lambda) + sum(log(ptt.f(lambda, t[3:N], tau))) + 2 * log(1 - beta(lambda, tau)) + sum(log(1 - beta(lambda, x[3:N])))
  
  lnL <- function(p)
    
  {
    
    lambda <- exp(p[1])
    
    -yule.f.lik(lambda)
    
  }
  
  res <- NULL
  
  temp <- suppressWarnings(optim(log(0.5), lnL))
  
  res$LH <- -temp$value
  
  res$r <- exp(temp$par[1])
  
  res$AIC <- -2 * res$LH + 2
  
  res
  
}



inv.logit <- function(x) return(1 / (exp(-x) + 1));

logit <- function(p) return(log(p) - log(1-p));

#Birth-death model with incomplete sampling

fit.bd.f <- function(phy, sampling.f = 1)
  
{
  
  N <- length(phy$tip.label)
  
  x <- c(0, branching.times(phy))
  
  tau <- x[2]
  
  t <- tau - x
  
  ptt.f <- function(r, eps, time, tau) ((eps - 1) * sampling.f * exp(r * (tau - time))) / ((-sampling.f * exp(r * (tau - time))) + eps + sampling.f - 1)
  
  beta <- function(r, eps, time) 1 - ptt.f(r, eps, 0, time) * exp(-r * time)
  
  bd.f.lik <- function(r, eps) lfactorial(N - 1) - (N - 2) * log(1 - eps) + (N - 2) * log(r) + sum(log(ptt.f(r, eps, t[3:N], tau))) + 2 * log(1 - beta(r, eps, tau)) + sum(log(1 - beta(r, eps, x[3:N])))
  
  lnL <- function(p)
    
  {
    
    r <- exp(p[1])
    
    eps <- inv.logit(p[2])
    
    -bd.f.lik(r, eps)
    
  }
  
  res <- NULL
  
  temp <- optim(c(log(0.1), logit(0.5)), lnL)
  
  res$LH <- -temp$value
  
  res$r <- exp(temp$par[1])
  
  res$eps <- inv.logit(temp$par[2])
  
  res$AIC <- -2 * res$LH + 4
  
  res
  
}



#Model of Strathmann and Slatkin 1983, with hyperbolic decay of high initial speciation rate and constant extinction

fit.ss83 <- function(phy)
  
{
  
  N <- length(phy$tip.label)
  
  x <- c(0, branching.times(phy))
  
  tau <- x[2]
  
  t <- tau - x
  
  lambda.t <- function(lambda0, mu, eps) lambda0 / (1 + eps * t[3:N])
  
  rho.t <- function(k, mu, eps, tau, t) mu * (((k * log((eps * t + 1) / (eps * tau + 1))) / eps) - t + tau)
  
  ptt.t <- function(k, mu, eps, t, tau)
    
  {
    
    foo <- function(x) exp(mu * (((k * log((eps * t + 1) / (eps * x + 1))) / eps) - t + x))
    
    integ <- integrate(foo, t, tau)
    
    res <- (1 / (1 + mu * integ$value))
    
    res
    
  }
  
  chi.t <- function(k, mu, eps, t, tau) ptt.t(k, mu, eps, t, tau) * exp(rho.t(k, mu, eps, tau, t))
  
  ss83.lik <- function(lambda0, mu, eps, k)
    
  {
    
    p1 <- 0
    
    p2 <- 0
    
    for(i in 3:N) p1 <- p1 + log(ptt.t(k, mu, eps, t[i], tau))
    
    for(i in 3:N) p2 <- p2 + log(chi.t(k, mu, eps, t[i], tau))
    
    -(lfactorial(N - 1) + sum(log(lambda.t(lambda0, mu, eps))) + p1 + log(chi.t(k, mu, eps, 0, tau) ^ 2) + p2)
    
  }
  
  lnL.ss83 <- function(p)
    
  {
    
    lambda0 <- exp(p[1])
    
    mu <- exp(p[2])
    
    k <- lambda0 / mu
    
    eps <- mu / (lambda0 + mu)
    
    ss83.lik(lambda0, mu, eps, k)
    
  }
  
  temp <- optim(c(log(0.5), log(0.1)), lnL.ss83)
  
  res <- NULL
  
  res$LH <- -temp$value
  
  res$lambda0 <- exp(temp$par[1]) + exp(temp$par[2])
  
  res$mu <- exp(temp$par[2])
  
  res$eps <- res$mu / res$lambda0
  
  res$k <- exp(temp$par[1]) / exp(temp$par[2])
  
  res$AIC <- 2 * -res$LH + 4
  
  res
  
}

######################################## load all the models before running misfitssingle
#####################################################################
##################################################################
##############################################################


misfitssingle<-function(tree,f,N0, alpha, lamb0,filename){
 
 if(missing(f)) {f <- 1
                  
  } else {
    
    f <- f }  

	
	sampling.f<-f

if(missing(N0)) {N0 <- 0
                   
  } else {
    
    N0 <- N0 }  
  
  if(missing(alpha)) {alpha <- 0.001
                      
  } else {
    
    alpha<-alpha}
	
if(missing(lamb0)) {lamb0 <- 0.1
                      
  } else {
    
    lamb0<-lamb0}
	
 
	 #basic estimate of birthdeath from Ape
  
  birthdeath(tree)->birthdeath
  
  

  
  
  #Model 1 constant turnover code---UNSTABLE
  
  lamb<-birthdeath$para[[2]]
  
  fitcoalMoranCST(tree, tau0=0.0001, meth = "Nelder-Mead", N0, Vtimes=FALSE)->Morlon_Model_1_CST
  
  Morlon_Model_1_CST
  
  
  
  #Model 2 exponetial turnover code--UNSTABLE
  
  fitcoalMoranEXP(tree, tau0=.0000001, gamma=.1, meth = "Nelder-Mead", N0, Vtimes=FALSE)->Morlon_Model_2_EXP
  
  
Morlon_Model_2_EXP
  
  
  
  #Model 3 constant sp, constant ext vary over time
  
  
fitcoalBD(tree,lamb0,alpha,mu0=0.01,beta=0,meth="Nelder-Mead",N0,cst.lamb=TRUE,cst.mu=TRUE,fix.eps=FALSE,mu.0=FALSE,pos=TRUE,Vtimes=FALSE)->Morlon_Model_3_BD
  
  Morlon_Model_3_BD
  
  
  
  #Model 4a rates of sp vary and ext constant over time
  
  fitcoalBD(tree,lamb0,alpha,mu0=0.01,beta=0,meth="Nelder-Mead",N0,cst.lamb=FALSE,cst.mu=TRUE,fix.eps=FALSE,mu.0=FALSE,pos=TRUE,Vtimes=FALSE)->Morlon_Model_4a_BD
  
  Morlon_Model_4a_BD
  
  
  
  
  
  #Model 4b constant sp, variable ext vary over time
  
  fitcoalBD(tree,lamb0,alpha,mu0=0.01,beta=0,meth="Nelder-Mead",N0,cst.lamb=TRUE,cst.mu=FALSE,fix.eps=FALSE,mu.0=FALSE,pos=TRUE,Vtimes=FALSE)->Morlon_Model_4b_BD
  
  Morlon_Model_4b_BD
  
  
  
  #Model 4c extinction fraction constant over time
  
  fitcoalBD(tree,lamb0,alpha,mu0=0.01,beta=0,meth="Nelder-Mead",N0,cst.lamb=FALSE,cst.mu=FALSE,fix.eps=TRUE,mu.0=FALSE,pos=TRUE,Vtimes=FALSE)->Morlon_Model_4c_BD
  
  Morlon_Model_4c_BD
  
  
  
  #Model 4d rates of sp/ext vary over time
  
  fitcoalBD(tree,lamb0,alpha,mu0=0.01,beta=0,meth="Nelder-Mead",N0,cst.lamb=FALSE,cst.mu=FALSE,fix.eps=FALSE,mu.0=FALSE,pos=TRUE,Vtimes=FALSE)->Morlon_Model_4d_BD
  
  Morlon_Model_4d_BD
  
  
  
  #Model 5 Constant and no extinction
  
  fitcoalBD(tree,lamb0,alpha,mu0=0.01,beta=0,meth="Nelder-Mead",N0,cst.lamb=TRUE,cst.mu=FALSE,fix.eps=FALSE,mu.0=TRUE,pos=TRUE,Vtimes=FALSE)->Morlon_Model_5
  
  Morlon_Model_5
  
  
  
  #Model 6 Varying speciation and no extinction
  
  
fitcoalBD(tree,lamb0,alpha,mu0=0.01,beta=0,meth="Nelder-Mead",N0,cst.lamb=FALSE,cst.mu=FALSE,fix.eps=FALSE,mu.0=TRUE,pos=TRUE,Vtimes=FALSE)->Morlon_Model_6
  
  Morlon_Model_6
  
  treebt<-branching.times(tree)
  
gamStat(treebt)
  
  
  
  #Fits to various models of rate changes
  
  fitdAICrc(treebt, modelset = c("pureBirth", "bd", "DDL", "DDX", "yule2rate", "yule3rate"), ints = 100)->dAIC
  
  dAIC$model->model
  
  dAIC$AIC->aic
  
  as.character(model)->model
  
  cbind(model, aic)->dAIC2
  
  
  
  #Fits a monotonic decay of speciation rate, sensu Rabosky & Lovette (2008) for comparison to LASER models (Yule2, DDL, DDX)above
  
  fit.linear(tree)->monotonic_decay
  
  monotonic_decay
  
  
  
  #Fits various change in sp/ext models from LASER
  
  fitSPVAR(treebt)->SPVAR
  
  SPVAR
  
  
  
  fitEXVAR(treebt)->EXVAR
  
  EXVAR
  
  
  
  fitBOTHVAR(treebt)->BOTHVAR
  
  BOTHVAR
  
  
  
  
  
  
  
  #Model of Strathmann and Slatkin 1983, with hyperbolic decay of high initial speciation rate and constant extinction
  
  fit.ss83 (tree)->hyperbolic_decay
  
  hyperbolic_decay
  
  
  
  rbind(dAIC2, c("Monotonic Decay", monotonic_decay$AIC))->dAIC2
  
  rbind(dAIC2, c("SPVAR", SPVAR$aic))->dAIC2
  
  rbind(dAIC2, c("EXVAR", EXVAR$aic))->dAIC2
  
  rbind(dAIC2, c("BOTHVAR", BOTHVAR$aic))->dAIC2
  
  rbind(dAIC2, c("Hyberbolic Decay", hyperbolic_decay$AIC))->dAIC2
  
  
  
  rbind(dAIC2, c("Morlon Model 1", Morlon_Model_1_CST$aicc))->dAIC2
  
  rbind(dAIC2, c("Morlon Model 2", Morlon_Model_2_EXP$aicc))->dAIC2
  
  rbind(dAIC2, c("Morlon Model 3", Morlon_Model_3_BD$aicc))->dAIC2
  
  rbind(dAIC2, c("Morlon Model 4a", Morlon_Model_4a_BD$aicc))->dAIC2
  
  rbind(dAIC2, c("Morlon Model 4b", Morlon_Model_4b_BD$aicc))->dAIC2
  
  rbind(dAIC2, c("Morlon Model 4c", Morlon_Model_4c_BD$aicc))->dAIC2
  
  rbind(dAIC2, c("Morlon Model 4d", Morlon_Model_4d_BD$aicc))->dAIC2
  
  rbind(dAIC2, c("Morlon Model 5", Morlon_Model_5$aicc))->dAIC2
  
  rbind(dAIC2, c("Morlon Model 6", Morlon_Model_6$aicc))->dAIC2
  
  AICc<-dAIC2
  
  as.numeric(AICc[1,2])->aic
  
  n<-length(tree$tip.label)-2
  
  aicc<-aic+(4/(n-2))
  
  AICc[1,2]<-aicc
  
  
  
  as.numeric(AICc[2,2])->aic
  
  n<-length(tree$tip.label)-2
  
  aicc<-aic+(12/(n-3))
  
  AICc[2,2]<-aicc
  
  
  
  as.numeric(AICc[3,2])->aic
  
  n<-length(tree$tip.label)-2
  
  aicc<-aic+(12/(n-3))
  
  AICc[3,2]<-aicc
  
  
  
  
  
  as.numeric(AICc[4,2])->aic
  
  n<-length(tree$tip.label)-2
  
  aicc<-aic+(12/(n-3))
  
  AICc[4,2]<-aicc
  
  
  
  as.numeric(AICc[5,2])->aic
  
  n<-length(tree$tip.label)-2
  
  aicc<-aic+(24/(n-4))
  
  AICc[5,2]<-aicc
  
  
  
  as.numeric(AICc[6,2])->aic
  
  n<-length(tree$tip.label)-2
  
  aicc<-aic+(12/(n-3))
  
  AICc[6,2]<-aicc
  
  
  
  as.numeric(AICc[7,2])->aic
  
  n<-length(tree$tip.label)-2
  
  aicc<-aic+(24/(n-4))
  
  AICc[7,2]<-aicc
  
  
  
  as.numeric(AICc[8,2])->aic
  
  n<-length(tree$tip.label)-2
  
  aicc<-aic+(24/(n-4))
  
  AICc[8,2]<-aicc
  
  
  
  as.numeric(AICc[9,2])->aic
  
  n<-length(tree$tip.label)-2
  
  aicc<-aic+(40/(n-6))
  
  AICc[9,2]<-aicc
  
  
  
  as.numeric(AICc[10,2])->aic
  
  n<-length(tree$tip.label)-2
  
  aicc<-aic+(12/(n-3))
  
  AICc[10,2]<-aicc
  
  
  
  colnames(AICc)<-c("Model", "AICc")
  
  as.data.frame(AICc)->AICc
  
  as.numeric(as.character(AICc$AICc))->x
  
  AICc[,2]<-x
  
  as.character(AICc[,1])->x
  
  AICc[,1]<-x
  
  AICc[sort.list(AICc[1:11,2]), ]->nonCoal
  
  AICc[12:20,]->Coal
  
  Coal[sort.list(Coal[,2]), ]->CoalSort
  
  rbind(CoalSort, nonCoal)->sortedAIC
  
  #sortedAIC[,2]<-x
  
  #sortedAIC[22:43]<-x
  
  output<-sortedAIC
  
  #return(output)
  
  output[,2]->x
  
  output[,1]->y
  
  akaike.weights(x[1:9])->morlonweights
  
  akaike.weights(x[10:20])->otherweights
  
morlonout<-cbind(y[1:9],x[1:9],morlonweights[[3]], morlonweights[[2]],morlonweights[[1]])
otherout<-cbind(y[10:20],x[10:20],otherweights[[3]],otherweights[[2]],otherweights[[1]])

head<-c("model", "aicc", "weights", "rel.LL", "deltaAIC")
colnames(morlonout)<-colnames(otherout)<-head
morlonout
otherout
out<-rbind(morlonout,otherout)

write.table(out,file=filename,sep="\t")
return(out)
 
}


 
