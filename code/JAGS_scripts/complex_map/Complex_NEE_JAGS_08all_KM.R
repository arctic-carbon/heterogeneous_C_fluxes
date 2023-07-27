### JAGS script for footprint decomposition using the hyperbolic GPP model.

# for the Caomplex version, four tundras: green, lichen, edge plateau, deg edge, 
# fen, water is water and water edge, deg p alone.
# priors for Reco come from Rd fit for 08all data.

model {
  
  # priors
  sigma_fit~dnorm(0.0001,(1/100000))T(0,25) #uninformed
  taufit<- 1/(sigma_fit^2)
  
  
  # Kormann & Meixner et al. dark respiration priors
  # fen priors
  betaFS1~dnorm(0.072035747,(1/0.05))T(0.004593992,0.200733346)
  alphaFS1~dnorm(1.068742209,(1/0.4))T(0.149340011,2.495072375)

  # degraded priors
  betaDS~dnorm(0.059554693,(1/0.05))T(0.002540218,0.260513051)
  alphaDS~dnorm(0.621223363,(1/0.35))T(0.015181369,2.363431537)
  # deg edge priors
  betaEDS~dnorm(0.067782126,(1/0.03))T(0.003119992,0.234228983)
  alphaEDS~dnorm(0.994700438,(1/0.5))T(0.060255369,2.6883232)

    # tundra priors
  betaLTS~dnorm(0.086534313,(1/0.02))T(0.005508464,0.255733424) 
  alphaLTS~dnorm(0.872168761,(1/0.2))T(0.09041995,2.21497457)
  betaGTS~dnorm(0.037893339,(1/0.02))T(0.002190368,0.127582137)
  alphaGTS~dnorm(1.798195275,(1/0.2))T(0.352229753,3.347085782)
  betaETS~dnorm(0.037919738,(1/0.02))T(0.002570578,0.095536941) #
  alphaETS~dnorm(2.233210215,(1/0.2))T(1.062652961,3.522043981)
  
  # water priors
  FWprior~dnorm(3.451312359,(1/1))T(0.182296829,9.785021861)
  
  # GPP less informative priors:

  E0LT~dunif(0,0.25)
  pmaxLT~dunif(0,35)
  E0GT~dunif(0,0.25)
  pmaxGT~dunif(0,35)
  E0ET~dunif(0,0.25)
  pmaxET~dunif(0,35)
  
  E0F1~dunif(0,0.25)
  pmaxF1~dunif(0,50)
  
  E0ED~dunif(0,0.25)
  pmaxED~dunif(0,60)
  # E0DS~dunif(0,0.25)
  # pmaxDS~dunif(0,60)
  

  
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
   # gppD[i]<-E0DS*pmaxDS*Tscale[i]*PAR[i]/(pmaxDS+(E0DS*PAR[i]))
    neeD[i]<-RsD[i]#- gppD[i]
    
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
    #gppD_hat[k]<-E0DS*pmaxDS*Tscale_hat[k]*PAR_hat[k]/(pmaxDS+(E0DS*PAR_hat[k]))
    
    neeLT_hat[k]<-RsLT_hat[k]- gppLT_hat[k]
    neeGT_hat[k]<-RsGT_hat[k]- gppGT_hat[k]
    neeET_hat[k]<-RsET_hat[k]- gppET_hat[k]
    neeF1_hat[k]<-RsF1_hat[k]- gppF1_hat[k]
    neeED_hat[k]<-RsED_hat[k]- gppED_hat[k]
    neeD_hat[k]<-RsD_hat[k]#- gppD_hat[k]
    
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