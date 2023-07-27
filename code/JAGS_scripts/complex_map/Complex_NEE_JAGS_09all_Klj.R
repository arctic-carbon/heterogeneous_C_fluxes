### JAGS script for footprint decomposition using the hyperbolic GPP model.

# for the Complex version, four tundras: green, lichen, edge plateau, deg edge, 
# fens, water is water and water edge, deg p alone.
# priors for Reco come from Rd fit for 09all data.

model {
  
  # priors
  sigma_fit~dnorm(0.0001,(1/100000))T(0,25) #uninformed
  taufit<- 1/(sigma_fit^2)
  
  
  # Kljun et al. dark respiration priors
  # fen priors
  betaFS1~dnorm(0.036279402,(1/0.05))T(0.001666431,0.121615537)
  alphaFS1~dnorm(2.32413041,(1/0.4))T(1.01197915,3.575616329)
 
  # degraded priors
  betaDS~dnorm(0.078626002,(1/0.05))T(0.003168272,0.335911922)
  alphaDS~dnorm(0.809034597,(1/0.35))T(0.022063699,2.604595056)

  # deg edge priors
  betaEDS~dnorm(0.039350764,(1/0.03))T(0.001525128,0.207606361)
  alphaEDS~dnorm(1.295730037,(1/0.5))T(0.150483295,2.675648354)
  
    # tundra priors
  betaLTS~dnorm(0.132711566,(1/0.02))T(0.004319434,0.364419565)
  alphaLTS~dnorm(0.240869873,(1/0.2))T(0.012408701,1.074971568)
  betaGTS~dnorm(0.025067651,(1/0.02))T(0.000900526,0.168049316)
  alphaGTS~dnorm(1.096133606,(1/0.2))T(0.061999263,2.425537449)
  betaETS~dnorm(0.011616745,(1/0.02))T(0.000506695,0.050061508) #
  alphaETS~dnorm(2.133299962,(1/0.2))T(1.303511664,2.820108153)
  
  # water priors
  FWprior~dnorm(3.320262022,(1/1))T(0.394629681,7.304476139)
  
  # GPP less informative priors:

  E0LT~dunif(0,0.75)
  pmaxLT~dunif(0,25)
  E0GT~dunif(0,0.75)
  pmaxGT~dunif(0,25)
  E0ET~dunif(0,0.25)
  pmaxET~dunif(0,25)
  
  E0F1~dunif(0,0.75)
  pmaxF1~dunif(0,25)
 
  E0ED~dunif(0,0.25)
  pmaxED~dunif(0,45)
  # E0DS~dunif(0,0.25)
  # pmaxDS~dunif(0,95)
  

  
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
    # gppD_hat[k]<-E0DS*pmaxDS*Tscale_hat[k]*PAR_hat[k]/(pmaxDS+(E0DS*PAR_hat[k]))
    
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