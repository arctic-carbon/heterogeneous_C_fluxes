### JAGS script for footprint decomposition using constant methane flux in time.

# for the Complex version, four tundras: green, lichen, edge plateau, deg edge, 
# one fens,  water is water and water edge, deg p alone.


model {
  
  # priors
  sigma_fit~dnorm(0.0001,(1/100000))T(0,25) #uninformed
  taufit<- 1/(sigma_fit^2)
  
  
   # fen priors
  ch4FS1~dunif(0,0.5)

  # degraded priors
  ch4DS~dunif(0,0.5)

  # deg edge priors
  ch4EDS~dunif(0,0.5)
  
  
  ch4LTS~dnorm(0,(1/0.00001))
  ch4GTS~dnorm(0,(1/0.00001))
  ch4ETS~dnorm(0,(1/0.00001))
  
  # water priors
  ch4FW~dunif(0,0.5)
  


  
  #likelihood
  for (i in 1:length(ch4)) {
  
    # Tower CH4
    TowerCH4[i]<-(ch4LTS*ftp_LT[i]+ch4GTS*ftp_GT[i]+ch4ETS*ftp_ET[i]+ch4FS1*ftp_F1[i]+ch4DS*ftp_D[i]+ch4FW*ftp_W[i]+ch4EDS*ftp_ED[i])
    ch4[i]~dnorm(TowerCH4[i],taufit)
    ch4_new[i]~dnorm(TowerCH4[i],taufit)
  
    sq.y[i]<-(ch4[i]-TowerCH4[i])^2
    sq.sim[i]<-(ch4_new[i]-TowerCH4[i])^2
    resid[i]<-(ch4[i]-ch4_new[i])

  }
  
  
  #Model fit
  sd.y<-sd(ch4)
  sd.sim<-sd(ch4_new)
  p.sd<-step(sd.sim-sd.y)

  mu.y<-mean(ch4)
  mu.sim<-mean(ch4_new)
  p.mu<-step(mu.sim-mu.y)

  fit<-sum(sq.y)
  fit.new<-sum(sq.sim)
  p.fit<-step(fit.new-fit)

  
}