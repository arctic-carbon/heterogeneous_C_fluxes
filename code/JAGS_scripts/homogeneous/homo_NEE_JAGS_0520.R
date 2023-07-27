### JAGS script for footprint decomposition using the hyperbolic GPP model.

model {
  
  # priors
  sigma_fit~dnorm(0.0001,(1/100000))T(0,25) #uninformed
  taufit<- 1/(sigma_fit^2)
  
  
  # dark respiration priors
  alpha~dnorm(1.200772842,(1/0.2))T(0.95192109,1.43790984)
  beta~dnorm(0.02986636,(1/0.02))T(0.002924673,0.06429345) 

  # GPP less informative priors:
  E0~dnorm(0.013045,(1/1))T(0,0.3)
  pmax~dnorm(10,(1/25))T(0,25)
  
  
  #likelihood
  for (i in 1:length(co2)) {
    # 
    Rs[i]<-alpha*exp(Tsoil[i]*beta)
    gpp[i]<-E0*pmax*Tscale[i]*PAR[i]/(pmax+(E0*PAR[i]))
    nee[i]<-Rs[i]- gpp[i]
 
    # Tower NEE
    TowerNEE[i]<-(nee[i])
    co2[i]~dnorm(TowerNEE[i],taufit)
    co2_new[i]~dnorm(TowerNEE[i],taufit)
  
    sq.y[i]<-(co2[i]-TowerNEE[i])^2
    sq.sim[i]<-(co2_new[i]-TowerNEE[i])^2
    resid[i]<-(co2[i]-co2_new[i])

  }
  
  #derived quantities
  Q10<-exp(10*beta)

  
  # predict to new drivers
  for (k in 1:length(T_hat)){
    Rs_hat[k]<-alpha*exp(T_hat[k]*beta)

    gpp_hat[k]<-E0*pmax*Tscale_hat[k]*PAR_hat[k]/(pmax+(E0*PAR_hat[k]))

    nee_hat[k]<-Rs_hat[k]- gpp_hat[k]

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