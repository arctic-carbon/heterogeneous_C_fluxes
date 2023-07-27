### JAGS script for footprint decomposition using the hyperbolic GPP model.

# for the Complex version, four tundras: green, lichen, edge plateau, deg edge, 
# fens, water is water and water edge, deg p alone.
# priors for Reco come from Rd fit for 0620 data.

model {
  
  # priors
  sigma_fit~dnorm(0.0001,(1/100000))T(0,25) #uninformed
  taufit<- 1/(sigma_fit^2)
  
  
  # Kljun et al. dark respiration priors
  # fen priors
  betaFS1~dnorm(0.063235891,(1/0.05))T(0.003791582,0.209601133)
  alphaFS1~dnorm(1.179524496,(1/0.4))T(0.088101774,2.82487014)

  # degraded priors
  betaDS~dnorm(0.101177614,(1/0.05))T(0.00510518,0.277500753)
  alphaDS~dnorm(0.961513627,(1/0.35))T(0.045017444,2.736757938)

  # deg edge priors
  betaEDS~dnorm(0.081743212,(1/0.03))T(0.010310792,0.160916679)
  alphaEDS~dnorm(1.913805445,(1/0.5))T(0.491365967,3.573897918)
  
    # tundra priors
  betaLTS~dnorm(0.048849607,(1/0.02))T(0.002348633,0.196662521)
  alphaLTS~dnorm(0.937868029,(1/0.2))T(0.037204568,2.521195497)
  betaGTS~dnorm(0.052613375,(1/0.02))T(0.002198715,0.230381981)
  alphaGTS~dnorm(0.577327148,(1/0.2))T(0.014221852,2.114433423)
  betaETS~dnorm(0.033021498,(1/0.02))T(0.001710534,0.130559122) #
  alphaETS~dnorm(1.180525406,(1/0.2))T(0.228068651,2.166921493)
  
  # water priors
  FWprior~dnorm(2.759526939,(1/1))T(0.125340807,9.74701211)
  
  # GPP less informative priors:

  E0LT~dunif(0,0.25)
  pmaxLT~dunif(0,25)
  E0GT~dunif(0,0.15)
  pmaxGT~dunif(0,25)
  E0ET~dunif(0,0.25)
  pmaxET~dunif(0,30)
  
  E0F1~dunif(0,0.25)
  pmaxF1~dunif(0,25)

  E0ED~dunif(0,0.25)
  pmaxED~dunif(0,45)
  E0DS~dunif(0,0.25)
  pmaxDS~dunif(0,95)
  

  
  #likelihood
  for (i in 1:length(co2)) {
    # lichen tundra
    RsLT[i]<-alphaLTS*exp(Tsoil[i]*betaLTS)
    gppLT[i]<-E0LT*pmaxLT*Tscale[i]*PAR[i]/(pmaxLT+(E0LT*PAR[i]))
    neeLT[i]<-RsLT[i]- gppLT[i]
    # green tundra
    RsGT[i]<-alphaGTS*exp(Tsoil[i]*betaGTS)
    gppGT[i]<-E0GT*pmaxGT*Tscale[i]*PAR[i]/(pmaxGT+(E0GT*PAR[i]))
    neeGT[i]<-RsGT[i]- gppGT[i]
    # edge tundra
    RsET[i]<-alphaETS*exp(Tsoil[i]*betaETS)
    gppET[i]<-E0ET*pmaxET*Tscale[i]*PAR[i]/(pmaxET+(E0ET*PAR[i]))
    neeET[i]<-RsET[i]- gppET[i]
    
    # fen 1
    RsF1[i]<-alphaFS1*exp(Tsoil[i]*betaFS1)
    gppF1[i]<-E0F1*pmaxF1*Tscale[i]*PAR[i]/(pmaxF1+(E0F1*PAR[i]))
    neeF1[i]<-RsF1[i]- gppF1[i]

    # degraded
    RsD[i]<-alphaDS*exp(Tsoil[i]*betaDS)
    gppD[i]<-E0DS*pmaxDS*Tscale[i]*PAR[i]/(pmaxDS+(E0DS*PAR[i]))
    neeD[i]<-RsD[i]- gppD[i]
    
     # degraded edge
    RsED[i]<-alphaEDS*exp(Tsoil[i]*betaEDS)
    gppED[i]<-E0ED*pmaxED*Tscale[i]*PAR[i]/(pmaxED+(E0ED*PAR[i]))
    neeED[i]<-RsED[i]- gppED[i]
    
    # water
    Fwater[i]<-FWprior
    
    # Tower NEE
    TowerNEE[i]<-(neeLT[i]*ftp_LT[i]+neeGT[i]*ftp_GT[i]+neeET[i]*ftp_ET[i]+neeF1[i]*ftp_F1[i]+neeD[i]*ftp_D[i]+Fwater[i]*ftp_W[i]+neeED[i]*ftp_ED[i])
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
  
  # predict to new drivers
  for (k in 1:length(T_hat)){
    RsLT_hat[k]<-alphaLTS*exp(T_hat[k]*betaLTS)
    RsGT_hat[k]<-alphaGTS*exp(T_hat[k]*betaGTS)
    RsET_hat[k]<-alphaETS*exp(T_hat[k]*betaETS)
    RsF1_hat[k]<-alphaFS1*exp(T_hat[k]*betaFS1)
    RsED_hat[k]<-alphaEDS*exp(T_hat[k]*betaEDS)
    RsD_hat[k]<-alphaDS*exp(T_hat[k]*betaDS)
    
    gppLT_hat[k]<-E0LT*pmaxLT*Tscale_hat[k]*PAR_hat[k]/(pmaxLT+(E0LT*PAR_hat[k]))
    gppGT_hat[k]<-E0GT*pmaxGT*Tscale_hat[k]*PAR_hat[k]/(pmaxGT+(E0GT*PAR_hat[k]))
    gppET_hat[k]<-E0ET*pmaxET*Tscale_hat[k]*PAR_hat[k]/(pmaxET+(E0ET*PAR_hat[k]))
    gppF1_hat[k]<-E0F1*pmaxF1*Tscale_hat[k]*PAR_hat[k]/(pmaxF1+(E0F1*PAR_hat[k]))
    gppED_hat[k]<-E0ED*pmaxED*Tscale_hat[k]*PAR_hat[k]/(pmaxED+(E0ED*PAR_hat[k]))
    gppD_hat[k]<-E0DS*pmaxDS*Tscale_hat[k]*PAR_hat[k]/(pmaxDS+(E0DS*PAR_hat[k]))
    
    neeLT_hat[k]<-RsLT_hat[k]- gppLT_hat[k]
    neeGT_hat[k]<-RsGT_hat[k]- gppGT_hat[k]
    neeET_hat[k]<-RsET_hat[k]- gppET_hat[k]
    neeF1_hat[k]<-RsF1_hat[k]- gppF1_hat[k]
    neeED_hat[k]<-RsED_hat[k]- gppED_hat[k]
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