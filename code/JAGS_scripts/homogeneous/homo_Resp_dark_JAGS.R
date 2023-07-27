### JAGS script for footprint decomposition

model {
  
  # priors
  sigma_fit~dnorm(0.0001,(1/100000))T(0,25) #uninformed
  taufit<- 1/(sigma_fit^2)

  alpha~dnorm(1,(1/1))T(-1,10)
  beta~dnorm(0.08,(1/0.5))T(0,10)


  
  #likelihood
  for (i in 1:length(co2)) {
    # 
    RsT[i]<-alpha*exp(Tsoil[i]*beta)

    # Tower NEE
    TowerNEE[i]<-(RsT[i])
    co2[i]~dnorm(TowerNEE[i],taufit)
    co2_new[i]~dnorm(TowerNEE[i],taufit)
  
    sq.y[i]<-(co2[i]-TowerNEE[i])^2
    sq.sim[i]<-(co2_new[i]-TowerNEE[i])^2
    resid[i]<-(co2[i]-co2_new[i])
    pd[i] <- dnorm(co2[i],TowerNEE[i],taufit)
    log_pd[i] <- log(dnorm(co2[i],TowerNEE[i],taufit))
    
  }
  
  #derived quantities
  Q10<-exp(10*beta)

  
  # predict to new soil temperatures
  for (k in 1:length(T_hat)){
    Rd_hat[k]<-alpha*exp(T_hat[k]*beta)
  }
  
  #Model fit
  sd.y<-sd(co2)
  sd.sim<-sd(co2_new)
  p.sd<-step(sd.sim-sd.y)

  mu.y<-mean(co2)
  mu.sim<-mean(co2_new)
  p.mu<-step(mu.sim-mu.y)

  fit<-sum(sq.y)
  fit.new<-sum(sq.sim)
  p.fit<-step(fit.new-fit)
  
  
}