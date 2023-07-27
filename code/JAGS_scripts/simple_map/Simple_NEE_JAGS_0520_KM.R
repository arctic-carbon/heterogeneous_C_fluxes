### JAGS script for footprint decomposition using the hyperbolic GPP model

# for the Simple version, all tundra lumped, all fen lumped, all water lumped, deg p alone.
# priors for Reco come from Rd fit for 0520 data.

model {
  
  # priors
  sigma_fit~dnorm(0.0001,(1/100000))T(0,25) #uninformed
  taufit<- 1/(sigma_fit^2)
  
  
  # Kormann & Meixner et al. dark respiration priors
  # fen priors
  betaFS1~dnorm(0.085045578,(1/0.05))T(0.003821461,0.339105789)
  alphaFS1~dnorm(0.045399473,(1/0.2))T(-0.683095726,0.995798458)
  # degraded priors
  betaDS~dnorm(0.140782735,(1/0.05))T(0.008164083,0.345750134)
  alphaDS~dnorm(1.139422661,(1/0.2))T(-0.243779134,2.919233687)
  # tundra priors
  betaLTS~dnorm(0.028494103,(1/0.02))T(0.002425661,0.066230167)
  alphaLTS~dnorm(1.321584635,(1/0.2))T(1.035666623,1.604903716)
  
  # water priors
  FWprior~dnorm(1.470009521,(1/1))T(0.073687682,4.782755122)
  
  
  # GPP less informative priors:
  E0LT~dunif(0,0.25)
  pmaxLT~dunif(0,50)
  
  E0F1~dunif(0,0.25)
  pmaxF1~dunif(0,50)
  
  E0DS~dunif(0,0.25)
  pmaxDS~dunif(0,75)
 

  
  #likelihood
  for (i in 1:length(co2)) {
    # tundra
    RsLT[i]<-alphaLTS*exp(Tsoil[i]*betaLTS)
    gppLT[i]<-Tscale[i]*E0LT*pmaxLT*PAR[i]/(pmaxLT+(E0LT*PAR[i]))
    neeLT[i]<-RsLT[i]- gppLT[i]
    
    # fen
    RsF1[i]<-alphaFS1*exp(Tsoil[i]*betaFS1)
    gppF1[i]<-Tscale[i]*E0F1*pmaxF1*PAR[i]/(pmaxF1+(E0F1*PAR[i]))
    neeF1[i]<-RsF1[i]- gppF1[i]
    
    # degraded
    RsD[i]<-alphaDS*exp(Tsoil[i]*betaDS)
    gppD[i]<-Tscale[i]*E0DS*pmaxDS*PAR[i]/(pmaxDS+(E0DS*PAR[i]))
    neeD[i]<-RsD[i]- gppD[i]
    
    # water
    Fwater[i]<-FWprior
    
    # Tower NEE
    TowerNEE[i]<-(neeLT[i]*ftp_LT[i]+neeF1[i]*ftp_F1[i]+neeD[i]*ftp_D[i]+Fwater[i]*ftp_W[i])
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

  # predict to new drivers
  for (k in 1:length(T_hat)){
    RsLT_hat[k]<-alphaLTS*exp(T_hat[k]*betaLTS)
    RsF1_hat[k]<-alphaFS1*exp(T_hat[k]*betaFS1)
    RsD_hat[k]<-alphaDS*exp(T_hat[k]*betaDS)
    
    gppLT_hat[k]<-E0LT*pmaxLT*Tscale_hat[k]*PAR_hat[k]/(pmaxLT+(E0LT*PAR_hat[k]))
    gppF1_hat[k]<-E0F1*pmaxF1*Tscale_hat[k]*PAR_hat[k]/(pmaxF1+(E0F1*PAR_hat[k]))
    gppD_hat[k]<-E0DS*pmaxDS*Tscale_hat[k]*PAR_hat[k]/(pmaxDS+(E0DS*PAR_hat[k]))
    neeLT_hat[k]<-RsLT_hat[k]- gppLT_hat[k]
    neeF1_hat[k]<-RsF1_hat[k]- gppF1_hat[k]
    neeD_hat[k]<-RsD_hat[k]- gppD_hat[k]

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