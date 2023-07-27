### JAGS script for footprint decomposition using the hyperbolic GPP model.

# for the COmplex version, four tundras: green, lichen, edge plateau, deg edge, 
# fen, water is water and water edge, deg p alone.

model {
  
  # priors
  sigma_fit~dnorm(0.0001,(1/100000))T(0,25) #uninformed
  taufit<- 1/(sigma_fit^2)
  
  # tundra priors

  # # tundra priors
  alphaLTS~dnorm(1,(1/1))T(0,10)
  betaLTS~dnorm(0.08,(1/0.5))T(0,1)
  alphaGTS~dnorm(1,(1/1))T(0,10)
  betaGTS~dnorm(0.08,(1/0.5))T(0,1)
  alphaETS~dnorm(1,(1/1))T(0,10)
  betaETS~dnorm(0.08,(1/0.5))T(0,1)

  # fen priors
  alphaFS1~dnorm(1,(1/1))T(0,10)
  betaFS1~dnorm(0.08,(1/0.5))T(0,1)

  # deg edge priors
  alphaEDS~dnorm(1,(1/1))T(0,10)
  betaEDS~dnorm(0.08,(1/0.5))T(0,1)

  # degraded priors
  alphaDS~dnorm(1,(1/1))T(0,10)
  betaDS~dnorm(0.08,(1/0.5))T(0,1)

  # water priors
  FWprior~dnorm(0.5,(1/100))T(0,25)
  
  #likelihood
  for (i in 1:length(co2)) {
    # lichen tundra
    RsLT[i]<-alphaLTS*exp(Tsoil[i]*betaLTS)
    # green tundra
    RsGT[i]<-alphaGTS*exp(Tsoil[i]*betaGTS)
    # edge tundra
    RsET[i]<-alphaETS*exp(Tsoil[i]*betaETS)

    # fen 1
    RsF1[i]<-alphaFS1*exp(Tsoil[i]*betaFS1)

    # degraded
    RsD[i]<-alphaDS*exp(Tsoil[i]*betaDS)
    # degraded edge
    RsED[i]<-alphaEDS*exp(Tsoil[i]*betaEDS)

    # water
    Fwater[i]<-FWprior
    
    # Tower NEE
    TowerNEE[i]<-(RsLT[i]*ftp_LT[i]+RsGT[i]*ftp_GT[i]+RsET[i]*ftp_ET[i]+RsF1[i]*ftp_F1[i]+RsD[i]*ftp_D[i]+Fwater[i]*ftp_W[i]+RsED[i]*ftp_ED[i])
    co2[i]~dnorm(TowerNEE[i],taufit)
    co2_new[i]~dnorm(TowerNEE[i],taufit)
  
    sq.y[i]<-(co2[i]-TowerNEE[i])^2
    sq.sim[i]<-(co2_new[i]-TowerNEE[i])^2
    resid[i]<-(co2[i]-co2_new[i])
  }
  
  #derived quantities
  Q10LT<-exp(10*betaLTS)
  Q10GT<-exp(10*betaGTS)
  Q10ET<-exp(10*betaETS)
  Q10F1<-exp(10*betaFS1)
  Q10D<-exp(10*betaDS)
  Q10ED<-exp(10*betaEDS)
  
  # predict to new soil temperatures
  for (k in 1:length(T_hat)){
    RsLT_hat[k]<-alphaLTS*exp(T_hat[k]*betaLTS)
    RsGT_hat[k]<-alphaGTS*exp(T_hat[k]*betaGTS)
    RsET_hat[k]<-alphaETS*exp(T_hat[k]*betaETS)
    RsF1_hat[k]<-alphaFS1*exp(T_hat[k]*betaFS1)
    RsED_hat[k]<-alphaEDS*exp(T_hat[k]*betaEDS)
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