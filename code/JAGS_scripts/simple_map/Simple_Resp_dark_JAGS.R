### JAGS script for footprint decomposition using the hyperbolic GPP model.

# for the Simple version, all tundra lumped, all fen lumped, all water lumped, deg p alone.

model {
  
  # priors
  sigma_fit~dnorm(0.0001,(1/100000))T(0,25) #uninformed
  taufit<- 1/(sigma_fit^2)
  
  # tundra priors
  alphaLTS~dnorm(1,(1/1))T(-1,10)
  betaLTS~dnorm(0.08,(1/0.5))T(0,10)

  # fen priors
  alphaFS1~dnorm(1,(1/1))T(-1,10)
  betaFS1~dnorm(0.08,(1/0.5))T(0,10)


  # degraded priors
  alphaDS~dnorm(1,(1/1))T(-1,10)
  betaDS~dnorm(0.08,(1/0.5))T(0,10)

  # water priors
  FWprior~dnorm(0.5,(1/100))T(0,25)
  
  
  #likelihood
  for (i in 1:length(co2)) {
    # tundra
    RsLT[i]<-alphaLTS*exp(Tsoil[i]*betaLTS)
    
    # fen
    RsF1[i]<-alphaFS1*exp(Tsoil[i]*betaFS1)


    # degraded
    RsD[i]<-alphaDS*exp(Tsoil[i]*betaDS)
    
    # water
    Fwater[i]<-FWprior
    
    # Tower NEE
    TowerNEE[i]<-(RsLT[i]*ftp_LT[i]+RsF1[i]*ftp_F1[i]+RsD[i]*ftp_D[i]+Fwater[i]*ftp_W[i])
    co2[i]~dnorm(TowerNEE[i],taufit)
    co2_new[i]~dnorm(TowerNEE[i],taufit)
  
    sq.y[i]<-(co2[i]-TowerNEE[i])^2
    sq.sim[i]<-(co2_new[i]-TowerNEE[i])^2
    resid[i]<-(co2[i]-co2_new[i])

  }
  

  #derived quantities
  Q10LT<-exp(10*betaLTS)
  Q10F1<-exp(10*betaFS1)
  Q10D<-exp(10*betaDS)

  # predict to new soil temperatures
  for (k in 1:length(T_hat)){
    RsLT_hat[k]<-alphaLTS*exp(T_hat[k]*betaLTS)
    RsF1_hat[k]<-alphaFS1*exp(T_hat[k]*betaFS1)
    RsD_hat[k]<-alphaDS*exp(T_hat[k]*betaDS)
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