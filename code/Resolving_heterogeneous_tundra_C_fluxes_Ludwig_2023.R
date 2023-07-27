
# R code used in 'Resolving heterogeneous fluxes from tundra halves the growing season carbon budget'

# Author: Sarah M. Ludwig 1,2

# 1Department of Earth and Environmental Science, Columbia University, New York, NY, United States of America
# 2Lamont-Doherty Earth Observatory, Palisades, NY, United States of America

# Date: July 2023

# This script runs the gap-filling models in MDS and Bayesian analyses.
# It requires artificial gap generating script.
# Each JAGS model is its own additional script, with priors, likelihoods, and model structure
# Code to assess posterior predictive checks and convergence is included..
# Calculated in this code but not exported are RMSE's from artificial gaps.
# Exports include:
  # 1. Summary files of parameters from dark respiration model fits.
  # 2. parameter posterior chains from NEE gap-filling models.
  # 3. predicted landcover NEE posterior chains from NEE gap-filling models.
# landcover and regional carbon budgets were calculated from the predicted landcover NEE posterior chains
# by summarizing over the desired time period and scaling up by the landcover ares provided in this script

# Running the JAGS models can take anywhere from a couple of minutes to a couple of hours, 
# depending on the # of parameters and the length of observations.
# It is unlikely any personal computer could hold all of the MCMC output in memory at once, therefore chains are exported for later analysis.
# It is unlikely any personal computer can hold all JAGS model objects in memory at once, so it is recommended 
# the script be run in chunks with results exported and memory cleared between runs.

# The priors in the JAGS model scripts are uninformative, 
# except in the case of respiration parameters in the NEE models, which are strict.
# For the NEE models, respiration parameter priors have the mean and sd of the posteriors from the dark only respiration model results
# and are further restricted by truncating to within the 95% probability of the posteriors.

#### Libraries ####

library(lubridate)
library(rjags)
library(runjags)
library(parallel)
library(matrixStats)
library(loo)
library(tidyverse)
library(tidyselect)
library(REddyProc)
library(xts)
library(bayestestR)

`%nin%` = Negate(`%in%`)

#### Data set up ####

# 1. scaling region areas by landcover in meter2

LTarea=40573597.9
GTarea=16632266.69
ETarea=24695035.2
EDarea=16394367.7
F1area=15946974.9
F2area=11328203.7
Warea=4211306.3+14920314.6
Darea=7319421.4
totalarea=LTarea+GTarea+ETarea+EDarea+F1area+F2area+Warea+Darea

# 2. set analysis date:

run_version='vdate_time' # set as needed

analysis_date='_0520' # change as needed 
# analysis_date=='_0620'
# analysis_date=='_07all'
# analysis_date=='_08all'
# analysis_date=='_09all'

JAGS_model=paste0('NEE_JAGS',analysis_date)


# 3. load flux timeseries and format properly
setwd("/Users/Ludda/Library/CloudStorage/GoogleDrive-sml2278@columbia.edu/My Drive/YKD/YKD_flux_processed/Footprint_grids/repo")

df_complex=read_csv('eddy_cov_C_fluxes_YKD_complexmap_footprints.csv')
df_simple=read_csv('eddy_cov_C_fluxes_YKD_simplemap_footprints.csv')

ts_end=as.POSIXct('2020-10-01 00:00:00',format="%Y-%m-%d %H:%M",tz='Etc/GMT+9')
ts_start=as.POSIXct('2019-04-30 23:30:00',format="%Y-%m-%d %H:%M",tz='Etc/GMT+9')

df_complex=subset(df_complex,df_complex$datetime<ts_end)
df_complex=subset(df_complex,df_complex$datetime>ts_start)
df_simple=subset(df_simple,df_simple$datetime<ts_end)
df_simple=subset(df_simple,df_simple$datetime>ts_start)

# subset timeseries for appropriate month of analysis
if (analysis_date=='_07all'){
  df_complex=subset(df_complex,df_complex$month%in%c('07'))
} else if (analysis_date=='_08all'){
  df_complex=subset(df_complex,df_complex$month%in%c('08'))
} else if (analysis_date=='_09all'){
  df_complex=subset(df_complex,df_complex$month%in%c('09'))
}else if (analysis_date=='_0520'){
  df_complex=subset(df_complex,df_complex$month%in%c('05'))
} else if (analysis_date=='_0620'){
  df_complex=subset(df_complex,df_complex$month%in%c('06'))
}
if (analysis_date=='_07all'){
  df_simple=subset(df_simple,df_simple$month%in%c('07'))
} else if (analysis_date=='_08all'){
  df_simple=subset(df_simple,df_simple$month%in%c('08'))
} else if (analysis_date=='_09all'){
  df_simple=subset(df_simple,df_simple$month%in%c('09'))
}else if (analysis_date=='_0520'){
  df_simple=subset(df_simple,df_simple$month%in%c('05'))
} else if (analysis_date=='_0620'){
  df_simple=subset(df_simple,df_simple$month%in%c('06'))
}

# 4. use the artificial_gap_creation.R script to create training and testing datasets with random stratified length gaps
# results include: 
  # a. training_complex and testing_complex dataframes used in subsequent Bayesian analysis for the complex map footprint influences
  # b. training_simple and testing_simple dataframes used in subsequent Bayesian analysis for the simple map footprint influences
  # c. observation masks used to properly line up predicted fluxes (_hat) with measured observations in figure creation

# 5. dark data filtered from datasets loaded in step 3.

# for the complex map footprint influences:
df_dark_complex=subset(df_complex,df_complex$PPFD_1_1_1<50)
df_dark_complex=subset(df_dark_complex,df_dark_complex$co2>-3)
df_dark_complex=subset(df_dark_complex,!is.na(df_dark_complex$air_temperature))
df_dark_complex=subset(df_dark_complex,!is.na(df_dark_complex$co2))
df_dark_complex=subset(df_dark_complex,!is.na(df_dark_complex$PPFD_1_1_1))
df_dark_complex=subset(df_dark_complex,!is.na(df_dark_complex$lichen_tundra_Klj))

# for the simple map footprint influences:
df_dark_simple=subset(df_simple,df_simple$PPFD_1_1_1<50)
df_dark_simple=subset(df_dark_simple,df_dark_simple$co2>-3)
df_dark_simple=subset(df_dark_simple,!is.na(df_dark_simple$air_temperature))
df_dark_simple=subset(df_dark_simple,!is.na(df_dark_simple$co2))
df_dark_simple=subset(df_dark_simple,!is.na(df_dark_simple$PPFD_1_1_1))
df_dark_simple=subset(df_dark_simple,!is.na(df_dark_simple$tundra_H))

# 6. load complete timeseries of drivers used for predictions of fluxes, format properly and filter to analysis timeframe. result is 'df_hat' dataframe

df_gapf=read_delim('flux_drivers_YKD_gapfilled.txt',delim="\t")%>%slice(-1)

df_gapf$datetime=as.POSIXct(df_gapf$`Date Time`,format="%Y-%m-%d %H:%M",tz="Etc/GMT+9")
Tair_hat=as.numeric(df_gapf$Tair_f)
PAR_hat=as.numeric(df_gapf$PPFD_1_1_1_f)
Tsoil_hat=as.numeric(df_gapf$SoilT_f)
df_hat1=as.data.frame(cbind(Tair_hat,PAR_hat,Tsoil_hat))
df_hat1$datetime=df_gapf$datetime
df_hat1$month=format(df_hat1$datetime,"%m")
Tmax=40
Tmin=-1.5
Topt=15
df_hat1$Tscale_hat=((df_hat1$Tair_hat-Tmin)*(df_hat1$Tair_hat-Tmax))/((df_hat1$Tair_hat-Tmin)*(df_hat1$Tair_hat-Tmax)-(df_hat1$Tair_hat-Topt)^2)
df_hat1=subset(df_hat1,df_hat1$datetime<ts_end)
df_hat1=subset(df_hat1,df_hat1$datetime>ts_start)

if (analysis_date=='_07all'){
  df_hat=subset(df_hat1,df_hat1$month%in%c('07'))
} else if (analysis_date=='_08all'){
  df_hat=subset(df_hat1,df_hat1$month%in%c('08'))
} else if (analysis_date=='_09all'){
  df_hat=subset(df_hat1,df_hat1$month%in%c('09'))
} else if (analysis_date=='_0520'){
  df_hat=subset(df_hat1,df_hat1$month%in%c('05'))
} else if (analysis_date=='_0620'){
  df_hat=subset(df_hat1,df_hat1$month%in%c('06'))
}


#### MDS NEE gap-filling ####

# Artificial gaps match the ones generated in simple and complex landcover map footprint influence datasets

# load eddy covariance dataset with naming conventions set up for MDS pacakge
EddyData2=read_csv('eddy_cov_C_fluxes_YKD_naming_for_MDS.csv')

#first data in 2020 in Bayse timeseries:
MDS_filter_date=as.POSIXct(data_2020$datetime[1],format="%Y-%m-%d %H:%M",tz="Etc/GMT+9")
EddyData3=subset(EddyData2,EddyData2$DateTime<MDS_filter_date)
additional_timeseries=length(EddyData3$DateTime)

train_MDS=EddyData2
for (i in 1:length(withheld_index)){
  indice=withheld_index[i]+additional_timeseries
  train_MDS$NEE[indice]=NA
}
test_MDS=EddyData2[withheld_index+additional_timeseries,]
test_MDS=subset(test_MDS,!is.na(test_MDS$NEE))

EProc=sEddyProc$new('YKD',train_MDS,c('NEE','Rg','Tair','VPD','PPFD_1_1_1','SoilT'))
EProc$sMDSGapFill('NEE') 

FilledEddyData=EProc$sExportResults()
CombinedData=cbind(train_MDS,FilledEddyData)
# fWriteDataframeToFile(CombinedData,paste0('YKD_MDS_train_gapfill',analysis_date,run_version,'.txt'))

# calculate MDS RMSE against withheld data from artificial gaps for analysis_date:
MDS_fit=left_join(test_MDS,subset(CombinedData,select=c(NEE_f,DateTime)),by='DateTime')
MDS_fit$residuals=(MDS_fit$NEE_f-MDS_fit$NEE)
RMSE_MDS=sqrt(sum(MDS_fit$residuals^2)/length(MDS_fit$residuals));RMSE_MDS;analysis_date 



#### Homogeneous NEE ####

# 1. set up initial values, change as needed
initsNEE=list(list(sigma_fit=1,alpha=0.7,beta=0.028,
                   E0=0.2*0.9,pmax=5.90*0.9),
              list(sigma_fit=2,alpha=1.25,beta=0.04,
                   E0=0.03038,pmax=15),
              list(sigma_fit=2.5,alpha=2,beta=0.055,
                   E0=0.003038*1.1,pmax=25))
initsRdark=list(list(sigma_fit=1,alpha=0.2,beta=0.04*0.9),
                list(sigma_fit=2,alpha=1,beta=0.08),
                list(sigma_fit=2.5,alpha=4,beta=0.18))

# 2. set up data lists for jags model
Tscale=((train_complex$air_temperature-Tmin)*(train_complex$air_temperature-Tmax))/((train_complex$air_temperature-Tmin)*(train_complex$air_temperature-Tmax)-(train_complex$air_temperature-Topt)^2)

data_jags_dark=list(Tsoil=df_dark_simple$air_temperature,
                    co2=df_dark_simple$co2,T_hat=df_hat$Tair_hat)
data_jags_NEE=list(Tsoil=train_complex$air_temperature,PAR=train_complex$PPFD_1_1_1,
                   co2=train_complex$co2,Tscale=Tscale,T_hat=df_hat$Tair_hat,
                   Tscale_hat=df_hat$Tscale_hat,PAR_hat=df_hat$PAR_hat)

# 3. jags run parameters
n.adapt=5000
n.update=5000
n.iter=3000

# 4. dark fit for respiration
set.seed(1)
cl=makeCluster(3)
rjm_dark=run.jags(model="homo_Resp_dark_JAGS.R",data=data_jags_dark,
                  monitor=c('sigma_fit','alpha','beta','Q10',
                            'p.fit','p.mu','p.sd'),
                  method='rjparallel',summarise = FALSE,cl=cl,
                  n.chains=length(initsRdark),inits=initsRdark,burnin=n.update, adapt=n.adapt,sample=n.iter)
stopCluster(cl)

# check model fit params:
# should be > 0.05 and < 0.95
mean(as.matrix(combine.mcmc(as.mcmc.list(rjm_dark,vars="p.fit"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjm_dark,vars="p.mu"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjm_dark,vars="p.sd"))))

# export posterior chains for parameters for dark model
alpha_quant=as.matrix(combine.mcmc(as.mcmc.list(rjm_dark,vars="alpha")))
beta_quant=as.matrix(combine.mcmc(as.mcmc.list(rjm_dark,vars="beta")))
Q10_quant=as.matrix(combine.mcmc(as.mcmc.list(rjm_dark,vars="Q10")))
params_list=list(alpha_quant,beta_quant,Q10_quant)
params_export=matrix(nrow=length(params_list),ncol=5)
for (i in 1:length(params_list)){
  iter=params_list[[i]]
  params_export[i,1]=attributes(iter)$dimnames[[2]]
  param_stats=quantile(iter,probs=c(0.025,0.5,0.975))
  params_export[i,2]=param_stats[2]
  params_export[i,3]=mean(iter)
  params_export[i,4]=param_stats[1]
  params_export[i,5]=param_stats[3]
}
params_export=as.data.frame(params_export)
names(params_export)=c('parameter','median_est','mean_est','lower2_5','upper98_5')
write_csv(params_export,paste0('HOMO_Rd_params',analysis_date,run_version,'.csv'))

# 5.  all light fit for NEE
set.seed(1)
cl=makeCluster(3)
rjm_NEE=run.jags(model=paste0("homo_",JAGS_model,".R"),data=data_jags_NEE,
                 monitor=c('sigma_fit','alpha','beta','pmax','E0',
                           'nee_hat','Q10','Rs_hat','gpp_hat',
                           'p.fit','p.mu','p.sd','log_pd','co2_new'),
                 method='rjparallel',summarise = FALSE,cl=cl,
                 n.chains=length(initsNEE),inits=initsNEE,burnin=n.update, adapt=n.adapt,sample=n.iter)
stopCluster(cl)
# check model fit params:
# should be > 0.05 and < 0.95
mean(as.matrix(combine.mcmc(as.mcmc.list(rjm_NEE,vars="p.fit"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjm_NEE,vars="p.mu"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjm_NEE,vars="p.sd"))))

# compare to withheld data
Homo_resid=t(as.matrix(combine.mcmc(as.mcmc.list(rjm_NEE,vars="nee_hat"))))[withheld_index,]-test_complex$co2
Homo_RMSE=sqrt(colSums((Homo_resid)^2)/nrow(Homo_resid))
hist(Homo_RMSE)

# export posterior chains for predictions for NEE model
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_NEE,vars="nee_hat"))))),paste0('nee_HOMO',analysis_date,run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_NEE,vars="gpp_hat"))))),paste0('gpp_HOMO',analysis_date,run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_NEE,vars="Rs_hat"))))),paste0('Rs_HOMO',analysis_date,run_version,'.csv'))

# export posterior chains for parameters for NEE model
alpha_quant=as.matrix(combine.mcmc(as.mcmc.list(rjm_NEE,vars="alpha")))
beta_quant=as.matrix(combine.mcmc(as.mcmc.list(rjm_NEE,vars="beta")))
Q10_quant=as.matrix(combine.mcmc(as.mcmc.list(rjm_NEE,vars="Q10")))
pmax_quant=as.matrix(combine.mcmc(as.mcmc.list(rjm_NEE,vars="pmax")))
E0_quant=as.matrix(combine.mcmc(as.mcmc.list(rjm_NEE,vars="E0")))
params_list=list(alpha_quant,beta_quant,Q10_quant,pmax_quant,E0_quant)
params_export=matrix(nrow=length(params_list),ncol=5)
for (i in 1:length(params_list)){
  iter=params_list[[i]]
  params_export[i,1]=attributes(iter)$dimnames[[2]]
  param_stats=quantile(iter,probs=c(0.025,0.5,0.975))
  params_export[i,2]=param_stats[2]
  params_export[i,3]=mean(iter)
  params_export[i,4]=param_stats[1]
  params_export[i,5]=param_stats[3]
}
params_export=as.data.frame(params_export)
names(params_export)=c('parameter','median_est','mean_est','lower2_5','upper98_5')
write_csv(params_export,paste0('HOMO_NEE_params',analysis_date,run_version,'.csv'))



#### Complex NEE ####

# 1. set up initial values, change as needed
initsNEE_complex=list(list(sigma_fit=1,alphaDS=1*0.9,betaDS=0.08*0.9,FWprior=2*0.9,E0DS=0.0013045*0.9,pmaxDS=5.90*0.9,
                    alphaLTS=1*0.9,betaLTS=0.04*0.9,E0LT=0.003038*0.9,pmaxLT=5.90*0.9,
                    alphaGTS=1*0.9,betaGTS=0.08*0.9,E0GT=0.003038*0.9,pmaxGT=5.90*0.9,
                    alphaETS=1*0.9,betaETS=0.03,E0ET=0.003038*0.9,pmaxET=5.90*0.9,
                    alphaFS1=1*0.9,betaFS1=0.07*0.9,E0F1=0.0013045*0.9,pmaxF1=5.90*0.9,
                    alphaEDS=1*0.9,betaEDS=0.03*0.9,E0ED=0.0013045*0.9,pmaxED=5.90*0.9),
               list(sigma_fit=2,alphaDS=1,betaDS=0.08,FWprior=2,E0DS=0.013045,pmaxDS=14.33,
                    alphaLTS=1,betaLTS=0.04,E0LT=0.03038,pmaxLT=14.33,
                    alphaGTS=1,betaGTS=0.08,E0GT=0.03038,pmaxGT=14.33,
                    alphaETS=1,betaETS=0.045,E0ET=0.03038,pmaxET=14.33,
                    alphaFS1=1,betaFS1=0.01,E0F1=0.013045,pmaxF1=14.33,
                    alphaEDS=1,betaEDS=0.03,E0ED=0.013045,pmaxED=14.33),
               list(sigma_fit=2.5,alphaDS=1*1.1,betaDS=0.08*1.1,FWprior=2*1.1,E0DS=0.1,pmaxDS=25,
                    alphaLTS=1*1.1,betaLTS=0.04*1.1,E0LT=0.13038*1.1,pmaxLT=25,
                    alphaGTS=1*1.1,betaGTS=0.08*1.1,E0GT=0.13038,pmaxGT=25,
                    alphaETS=1*1.1,betaETS=0.01,E0ET=0.13038*1.1,pmaxET=25,
                    alphaFS1=1*1.1,betaFS1=0.05*1.1,E0F1=0.13045*1.1,pmaxF1=25,
                    alphaEDS=1*1.1,betaEDS=0.03*1.1,E0ED=0.13045*1.1,pmaxED=25))

initsRdark_complex=list(list(sigma_fit=1,alphaDS=0.1*0.9,betaDS=0.02*0.9,FWprior=0.2*0.9,
                      alphaLTS=0.1*0.9,betaLTS=0.02*0.9,
                      alphaGTS=0.1*0.9,betaGTS=0.02*0.9,
                      alphaETS=0.1*0.9,betaETS=0.02*0.9,
                      alphaFS1=0.1*0.9,betaFS1=0.02*0.9,
                      alphaEDS=0.1*0.9,betaEDS=0.02*0.9),
                 list(sigma_fit=2,alphaDS=1,betaDS=0.08,FWprior=2,
                      alphaLTS=1,betaLTS=0.08,
                      alphaGTS=1,betaGTS=0.08,
                      alphaETS=1,betaETS=0.08,
                      alphaFS1=1,betaFS1=0.08,
                      alphaEDS=1,betaEDS=0.08),
                 list(sigma_fit=2.5,alphaDS=5*1.1,betaDS=0.8*1.1,FWprior=12*1.1,
                      alphaLTS=5*1.1,betaLTS=0.8*1.1,
                      alphaGTS=5*1.1,betaGTS=0.8*1.1,
                      alphaETS=5*1.1,betaETS=0.8*1.1,
                      alphaFS1=5*1.1,betaFS1=0.8*1.1,
                      alphaEDS=5*1.1,betaEDS=0.8*1.1))

# 2. jags run parameters
n.adapt=5000
n.update=5000
n.iter=3000

# 3. set up data lists for jags model
Tscale=((train_complex$air_temperature-Tmin)*(train_complex$air_temperature-Tmax))/((train_complex$air_temperature-Tmin)*(train_complex$air_temperature-Tmax)-(train_complex$air_temperature-Topt)^2)

data_complex_ZM=list(Tsoil=train_complex$air_temperature,PAR=train_complex$PPFD_1_1_1,
                   co2=train_complex$co2,ftp_LT=train_complex$lichen_tundra_ZM,ftp_GT=train_complex$green_tundra_ZM,ftp_ET=train_complex$edge_tundra_ZM,
                   ftp_F1=(train_complex$fen1_ZM+train_complex$fen2_ZM),ftp_D=train_complex$degraded_ZM,
                   ftp_W=train_complex$water_ZM,Tscale=Tscale,ftp_ED=train_complex$edge_degraded_ZM,T_hat=df_hat$Tair_hat,
                   Tscale_hat=df_hat$Tscale_hat,PAR_hat=df_hat$PAR_hat)
data_complex_H=list(Tsoil=train_complex$air_temperature,PAR=train_complex$PPFD_1_1_1,
                  co2=train_complex$co2,ftp_LT=train_complex$lichen_tundra_H,ftp_GT=train_complex$green_tundra_H,ftp_ET=train_complex$edge_tundra_H,
                  ftp_F1=(train_complex$fen1_H+train_complex$fen2_H),ftp_D=train_complex$degraded_H,
                  ftp_W=train_complex$water_H,Tscale=Tscale,ftp_ED=train_complex$edge_degraded_H,T_hat=df_hat$Tair_hat,
                  Tscale_hat=df_hat$Tscale_hat,PAR_hat=df_hat$PAR_hat)
data_complex_Klj=list(Tsoil=train_complex$air_temperature,PAR=train_complex$PPFD_1_1_1,
                    co2=train_complex$co2,ftp_LT=train_complex$lichen_tundra_Klj,ftp_GT=train_complex$green_tundra_Klj,ftp_ET=train_complex$edge_tundra_Klj,
                    ftp_F1=(train_complex$fen1_Klj+train_complex$fen2_Klj),ftp_D=train_complex$degraded_Klj,
                    ftp_W=train_complex$water_Klj,Tscale=Tscale,ftp_ED=train_complex$edge_degraded_Klj,T_hat=df_hat$Tair_hat,
                    Tscale_hat=df_hat$Tscale_hat,PAR_hat=df_hat$PAR_hat)

data_complex_dark_ZM=list(Tsoil=df_dark_complex$air_temperature,ftp_GT=df_dark_complex$green_tundra_ZM,ftp_ET=df_dark_complex$edge_tundra_ZM,
                        co2=df_dark_complex$co2,ftp_LT=df_dark_complex$lichen_tundra_ZM,ftp_F1=(df_dark_complex$fen1_ZM+df_dark_complex$fen2_ZM),ftp_D=df_dark_complex$degraded_ZM,
                        ftp_W=df_dark_complex$water_ZM,ftp_ED=df_dark_complex$edge_degraded_ZM,T_hat=df_hat$Tair_hat)
data_complex_dark_H=list(Tsoil=df_dark_complex$air_temperature,ftp_GT=df_dark_complex$green_tundra_H,ftp_ET=df_dark_complex$edge_tundra_H,
                       co2=df_dark_complex$co2,ftp_LT=df_dark_complex$lichen_tundra_H,ftp_F1=(df_dark_complex$fen1_H+df_dark_complex$fen2_H),ftp_D=df_dark_complex$degraded_H,
                       ftp_W=df_dark_complex$water_H,ftp_ED=df_dark_complex$edge_degraded_H,T_hat=df_hat$Tair_hat)
data_complex_dark_Klj=list(Tsoil=df_dark_complex$air_temperature,ftp_GT=df_dark_complex$green_tundra_Klj,ftp_ET=df_dark_complex$edge_tundra_Klj,
                         co2=df_dark_complex$co2,ftp_LT=df_dark_complex$lichen_tundra_Klj,ftp_F1=(df_dark_complex$fen1_Klj+df_dark_complex$fen2_Klj),ftp_D=df_dark_complex$degraded_Klj,
                         ftp_W=df_dark_complex$water_Klj,ftp_ED=df_dark_complex$edge_degraded_Klj,T_hat=df_hat$Tair_hat)


# 4. dark fit for respiration
set.seed(1)
cl=makeCluster(3)
rjmd_complex_ZM=run.jags(model="Complex_Resp_dark_JAGS.R",data=data_complex_dark_ZM,
                  monitor=c('sigma_fit','alphaDS','betaDS','Q10D',
                            'alphaGTS','betaGTS','alphaETS','betaETS',
                            'Q10GT','Q10ET',
                            'alphaLTS','betaLTS','alphaFS1','Q10LT','FWprior',
                            'betaFS1','Q10F1','Q10ED',
                            'p.fit','p.mu','p.sd',
                            'betaEDS','alphaEDS'),
                  method='rjparallel',summarise = FALSE,cl=cl,
                  n.chains=length(initsRdark_complex),inits=initsRdark_complex,burnin=n.update, adapt=n.adapt,sample=n.iter)
stopCluster(cl)

set.seed(1)
cl=makeCluster(3)
rjmd_complex_H=run.jags(model="Complex_Resp_dark_JAGS.R",data=data_complex_dark_H,
                 monitor=c('sigma_fit','alphaDS','betaDS','Q10D',
                           'alphaGTS','betaGTS','alphaETS','betaETS',
                           'Q10GT','Q10ET',
                           'alphaLTS','betaLTS','alphaFS1','Q10LT','FWprior',
                           'betaFS1','Q10F1','Q10ED',
                           'p.fit','p.mu','p.sd',
                           'betaEDS','alphaEDS'),
                 method='rjparallel',summarise = FALSE,cl=cl,
                 n.chains=length(initsRdark_complex),inits=initsRdark_complex,burnin=n.update, adapt=n.adapt,sample=n.iter)
stopCluster(cl)

set.seed(1)
cl=makeCluster(3)
rjmd_complex_Klj=run.jags(model="Complex_Resp_dark_JAGS.R",data=data_complex_dark_Klj,
                   monitor=c('sigma_fit','alphaDS','betaDS','Q10D',
                             'alphaGTS','betaGTS','alphaETS','betaETS',
                             'Q10GT','Q10ET',
                             'alphaLTS','betaLTS','alphaFS1','Q10LT','FWprior',
                             'betaFS1','Q10F1','Q10ED',
                             'p.fit','p.mu','p.sd',
                             'betaEDS','alphaEDS'),
                   method='rjparallel',summarise = FALSE,cl=cl,
                   n.chains=length(initsRdark_complex),inits=initsRdark_complex,burnin=n.update, adapt=n.adapt,sample=n.iter)
stopCluster(cl)

# check model fit params:
# should be > 0.05 and < 0.95
mean(as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_ZM,vars="p.fit"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_ZM,vars="p.mu"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_ZM,vars="p.sd"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_H,vars="p.fit"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_H,vars="p.mu"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_H,vars="p.sd"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_Klj,vars="p.fit"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_Klj,vars="p.mu"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_Klj,vars="p.sd"))))

# Convergence tests:
# should be close to 1
flux_paramsd_complex=c('sigma_fit','alphaDS','betaDS',
                'alphaGTS','betaGTS','alphaETS','betaETS',
                'alphaLTS','betaLTS','alphaFS1','FWprior',
                'betaFS1','betaEDS','alphaEDS')
# plot(rjmd_complex_Klj,vars=flux_paramsd_complex)
# plot(rjmd_complex_ZM,vars=flux_paramsd_complex)
# plot(rjmd_complex_H,vars=flux_paramsd_complex)

gelman_listd_complex_ZM=list()
for (i in 1:length(flux_paramsd_complex)){
  param=flux_paramsd_complex[i]
  result=gelman.diag(as.mcmc.list(rjmd_complex_ZM,vars=param))
  gelman_listd_complex_ZM[i]=result$psrf[1]
};gelman_listd_complex_ZM
gelman_listd_complex_H=list()
for (i in 1:length(flux_paramsd_complex)){
  param=flux_paramsd_complex[i]
  result=gelman.diag(as.mcmc.list(rjmd_complex_H,vars=param))
  gelman_listd_complex_H[i]=result$psrf[1]
};gelman_listd_complex_H
gelman_listd_complex_Klj=list()
for (i in 1:length(flux_paramsd_complex)){
  param=flux_paramsd_complex[i]
  result=gelman.diag(as.mcmc.list(rjmd_complex_Klj,vars=param))
  gelman_listd_complex_Klj[i]=result$psrf[1]
};gelman_listd_complex_Klj


# export posterior chains for parameters for dark model

# Parameters for ZM model
alphaDS_complex_ZM=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_ZM,vars="alphaDS")))
betaDS_complex_ZM=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_ZM,vars="betaDS")))
alphaLTS_complex_ZM=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_ZM,vars="alphaLTS")))
betaLTS_complex_ZM=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_ZM,vars="betaLTS")))
alphaGTS_complex_ZM=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_ZM,vars="alphaGTS")))
betaGTS_complex_ZM=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_ZM,vars="betaGTS")))
alphaETS_complex_ZM=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_ZM,vars="alphaETS")))
betaETS_complex_ZM=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_ZM,vars="betaETS")))
alphaFS1_complex_ZM=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_ZM,vars="alphaFS1")))
betaFS1_complex_ZM=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_ZM,vars="betaFS1")))
Q10LT_complex_ZM=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_ZM,vars="Q10LT")))
Q10GT_complex_ZM=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_ZM,vars="Q10GT")))
Q10ET_complex_ZM=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_ZM,vars="Q10ET")))
Q10F1_complex_ZM=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_ZM,vars="Q10F1")))
Q10D_complex_ZM=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_ZM,vars="Q10D")))
Q10ED_complex_ZM=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_ZM,vars="Q10ED")))
betaEDS_complex_ZM=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_ZM,vars="betaEDS")))
alphaEDS_complex_ZM=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_ZM,vars="alphaEDS")))
FWprior_complex_ZM=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_ZM,vars="FWprior")))

# Parameters for Hsieh model
alphaDS_complex_H=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_H,vars="alphaDS")))
betaDS_complex_H=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_H,vars="betaDS")))
alphaLTS_complex_H=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_H,vars="alphaLTS")))
betaLTS_complex_H=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_H,vars="betaLTS")))
alphaGTS_complex_H=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_H,vars="alphaGTS")))
betaGTS_complex_H=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_H,vars="betaGTS")))
alphaETS_complex_H=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_H,vars="alphaETS")))
betaETS_complex_H=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_H,vars="betaETS")))
alphaFS1_complex_H=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_H,vars="alphaFS1")))
betaFS1_complex_H=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_H,vars="betaFS1")))
Q10LT_complex_H=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_H,vars="Q10LT")))
Q10GT_complex_H=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_H,vars="Q10GT")))
Q10ET_complex_H=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_H,vars="Q10ET")))
Q10F1_complex_H=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_H,vars="Q10F1")))
Q10D_complex_H=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_H,vars="Q10D")))
Q10ED_complex_H=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_H,vars="Q10ED")))
betaEDS_complex_H=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_H,vars="betaEDS")))
alphaEDS_complex_H=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_H,vars="alphaEDS")))
FWprior_complex_H=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_H,vars="FWprior")))

# Parameters for Kljun model
alphaDS_complex_Klj=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_Klj,vars="alphaDS")))
betaDS_complex_Klj=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_Klj,vars="betaDS")))
alphaLTS_complex_Klj=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_Klj,vars="alphaLTS")))
betaLTS_complex_Klj=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_Klj,vars="betaLTS")))
alphaGTS_complex_Klj=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_Klj,vars="alphaGTS")))
betaGTS_complex_Klj=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_Klj,vars="betaGTS")))
alphaETS_complex_Klj=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_Klj,vars="alphaETS")))
betaETS_complex_Klj=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_Klj,vars="betaETS")))
alphaFS1_complex_Klj=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_Klj,vars="alphaFS1")))
betaFS1_complex_Klj=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_Klj,vars="betaFS1")))
Q10LT_complex_Klj=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_Klj,vars="Q10LT")))
Q10GT_complex_Klj=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_Klj,vars="Q10GT")))
Q10ET_complex_Klj=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_Klj,vars="Q10ET")))
Q10F1_complex_Klj=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_Klj,vars="Q10F1")))
Q10D_complex_Klj=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_Klj,vars="Q10D")))
Q10ED_complex_Klj=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_Klj,vars="Q10ED")))
betaEDS_complex_Klj=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_Klj,vars="betaEDS")))
alphaEDS_complex_Klj=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_Klj,vars="alphaEDS")))
FWprior_complex_Klj=as.matrix(combine.mcmc(as.mcmc.list(rjmd_complex_Klj,vars="FWprior")))

# Hsieh Respiration
params_list=list(Q10ED_complex_H,Q10D_complex_H,Q10F1_complex_H,Q10LT_complex_H,Q10GT_complex_H,Q10ET_complex_H,
                 betaFS1_complex_H,alphaFS1_complex_H,betaDS_complex_H,alphaDS_complex_H,
                 betaEDS_complex_H,alphaEDS_complex_H,betaLTS_complex_H,alphaLTS_complex_H,betaGTS_complex_H,alphaGTS_complex_H,
                 betaETS_complex_H,alphaETS_complex_H,FWprior_complex_H)
params_export=matrix(nrow=length(params_list),ncol=5)
for (i in 1:length(params_list)){
  iter=params_list[[i]]
  params_export[i,1]=attributes(iter)$dimnames[[2]]
  param_stats=quantile(iter,probs=c(0.025,0.5,0.975))
  params_export[i,2]=param_stats[2]
  params_export[i,3]=mean(iter)
  params_export[i,4]=param_stats[1]
  params_export[i,5]=param_stats[3]
}
params_export=as.data.frame(params_export)
names(params_export)=c('parameter','median_est','mean_est','lower2_5','upper97_5')
write_csv(params_export,paste0('Hsieh_Rd_params',analysis_date,'_complex_',run_version,'.csv'))

# KM Respiration
params_list=list(Q10ED_complex_ZM,Q10D_complex_ZM,Q10F1_complex_ZM,Q10LT_complex_ZM,Q10GT_complex_ZM,Q10ET_complex_ZM,
                 betaFS1_complex_ZM,alphaFS1_complex_ZM,betaDS_complex_ZM,alphaDS_complex_ZM,
                 betaEDS_complex_ZM,alphaEDS_complex_ZM,betaLTS_complex_ZM,alphaLTS_complex_ZM,betaGTS_complex_ZM,alphaGTS_complex_ZM,
                 betaETS_complex_ZM,alphaETS_complex_ZM,FWprior_complex_ZM)
params_export=matrix(nrow=length(params_list),ncol=5)
for (i in 1:length(params_list)){
  iter=params_list[[i]]
  params_export[i,1]=attributes(iter)$dimnames[[2]]
  param_stats=quantile(iter,probs=c(0.025,0.5,0.975))
  params_export[i,2]=param_stats[2]
  params_export[i,3]=mean(iter)
  params_export[i,4]=param_stats[1]
  params_export[i,5]=param_stats[3]
}
params_export=as.data.frame(params_export)
names(params_export)=c('parameter','median_est','mean_est','lower2_5','upper97_5')
write_csv(params_export,paste0('KormannMeixner_Rd_params',analysis_date,'_complex_',run_version,'.csv'))

# KLjun Respiration
params_list=list(Q10ED_complex_Klj,Q10D_complex_Klj,Q10F1_complex_Klj,Q10LT_complex_Klj,Q10GT_complex_Klj,Q10ET_complex_Klj,
                 betaFS1_complex_Klj,alphaFS1_complex_Klj,betaDS_complex_Klj,alphaDS_complex_Klj,
                 betaEDS_complex_Klj,alphaEDS_complex_Klj,betaLTS_complex_Klj,alphaLTS_complex_Klj,betaGTS_complex_Klj,alphaGTS_complex_Klj,
                 betaETS_complex_Klj,alphaETS_complex_Klj,FWprior_complex_Klj)
params_export=matrix(nrow=length(params_list),ncol=5)
for (i in 1:length(params_list)){
  iter=params_list[[i]]
  params_export[i,1]=attributes(iter)$dimnames[[2]]
  param_stats=quantile(iter,probs=c(0.025,0.5,0.975))
  params_export[i,2]=param_stats[2]
  params_export[i,3]=mean(iter)
  params_export[i,4]=param_stats[1]
  params_export[i,5]=param_stats[3]
}
params_export=as.data.frame(params_export)
names(params_export)=c('parameter','median_est','mean_est','lower2_5','upper97_5')
write_csv(params_export,paste0('Kljun_Rd_params',analysis_date,'_complex_',run_version,'.csv'))


# 5.  all light fit for NEE

full_var_list=c('sigma_fit','Q10D','Q10GT','Q10ET', 'Q10LT','FWprior',
                'Q10F1','Q10ED','p.fit','p.mu','p.sd',
                'pmaxLT','pmaxGT','pmaxET','pmaxF1','pmaxED',
                'E0LT','E0GT','E0ET','E0F1','E0ED',
                'alphaFS2','betaFS2','E0DS','pmaxDS',
                'alphaGTS','betaGTS','alphaETS','betaETS',
                'alphaDS','betaDS', 'alphaFS1','betaFS1','alphaLTS','betaLTS',
                'gppLT_hat','gppGT_hat','gppET_hat','gppF1_hat','gppED_hat',
                'neeLT_hat','neeGT_hat','neeET_hat','neeF1_hat','neeED_hat',
                'neeD_hat','gppD_hat',
                'betaEDS','alphaEDS','RsLT_hat','RsF1_hat',
                'RsED_hat','RsD_hat','RsGT_hat','RsET_hat')

set.seed(1)
cl=makeCluster(3)
rjm_complex_ZM=run.jags(model=paste0("Complex_",JAGS_model,"_KM.R"),data=data_complex_ZM,
                 monitor=full_var_list,
                 method='rjparallel',summarise = FALSE,cl=cl,
                 n.chains=length(initsNEE_complex),inits=initsNEE_complex,burnin=n.update, adapt=n.adapt,sample=n.iter)
stopCluster(cl)

set.seed(1)
cl=makeCluster(3)
rjm_complex_H=run.jags(model=paste0("Complex_",JAGS_model,"_H.R"),data=data_complex_H,
                monitor=full_var_list,
                method='rjparallel',summarise = FALSE,cl=cl,
                n.chains=length(initsNEE_complex),inits=initsNEE_complex,burnin=n.update, adapt=n.adapt,sample=n.iter)
stopCluster(cl)

set.seed(1)
cl=makeCluster(3)
rjm_complex_Klj=run.jags(model=paste0("Complex_",JAGS_model,"_Klj.R"),data=data_complex_Klj,
                  monitor=full_var_list,
                  method='rjparallel',summarise = FALSE,cl=cl,
                  n.chains=length(initsNEE_complex),inits=initsNEE_complex,burnin=n.update, adapt=n.adapt,sample=n.iter)
stopCluster(cl)

# check model fit params:
# should be > 0.05 and < 0.95
mean(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="p.fit"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="p.mu"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="p.sd"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="p.fit"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="p.mu"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="p.sd"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="p.fit"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="p.mu"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="p.sd"))))

# Convergence tests
flux_params=c('FWprior',
              'sigma_fit',
              'alphaGTS','betaGTS','alphaETS','betaETS',
              'alphaDS','betaDS', 'alphaFS1','betaFS1','alphaLTS','betaLTS',
              'pmaxLT','pmaxGT','pmaxF1','pmaxED',
              'E0LT','E0GT','E0F1','E0ED',
              'betaEDS','alphaEDS','E0ET','pmaxET','E0DS','pmaxDS')
# comment out E0DS and pmaxDS for 08all and 09all runs

# plot(rjm_complex_ZM,vars=flux_params)
# plot(rjm_complex_H,vars=flux_params)
# plot(rjm_complex_Klj,vars=flux_params)

gelman_complex_list_ZM=list()
for (i in 1:length(flux_params)){
  param=flux_params[i]
  result=gelman.diag(as.mcmc.list(rjm_complex_ZM,vars=param))
  gelman_complex_list_ZM[i]=result$psrf[1]
};gelman_complex_list_ZM

gelman_complex_list_H=list()
for (i in 1:length(flux_params)){
  param=flux_params[i]
  result=gelman.diag(as.mcmc.list(rjm_complex_H,vars=param))
  gelman_complex_list_H[i]=result$psrf[1]
};gelman_complex_list_H
gelman_complex_list_Klj=list()
for (i in 1:length(flux_params)){
  param=flux_params[i]
  result=gelman.diag(as.mcmc.list(rjm_complex_Klj,vars=param))
  gelman_complex_list_Klj[i]=result$psrf[1]
};gelman_complex_list_Klj

# compare to withheld data

# Kljun 
pred_neeLT_complex_Klj=t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="neeLT_hat")))[,withheld_index])*test_complex$lichen_tundra_Klj
pred_neeGT_complex_Klj=t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="neeGT_hat")))[,withheld_index])*test_complex$green_tundra_Klj
pred_neeET_complex_Klj=t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="neeET_hat")))[,withheld_index])*test_complex$edge_tundra_Klj
pred_neeED_complex_Klj=t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="neeED_hat")))[,withheld_index])*test_complex$edge_degraded_Klj
pred_neeF1_complex_Klj=t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="neeF1_hat")))[,withheld_index])*(test_complex$fen1_Klj+test_complex$fen2_Klj)
pred_neeD_complex_Klj=t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="neeD_hat")))[,withheld_index])*test_complex$degraded_Klj
pred_water_complex_Klj=matrix(nrow=nrow(pred_neeLT_complex_Klj),ncol=ncol(pred_neeLT_complex_Klj),data=1)
for (i in 1:nrow(pred_water_complex_Klj)){
  pred_water_complex_Klj[i,]=as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="FWprior")))*test_complex$water_Klj[i]
}
pred_NEE_complex_Klj=pred_neeLT_complex_Klj+pred_neeGT_complex_Klj+pred_neeET_complex_Klj+pred_neeED_complex_Klj+pred_neeF1_complex_Klj+pred_neeD_complex_Klj+pred_water_complex_Klj
residuals_NEE_complex_Klj=pred_NEE_complex_Klj-test_complex$co2 
residuals_NEE_complex_Klj=na.omit(residuals_NEE_complex_Klj)
RMSE_NEE_complex_Klj=sqrt(colSums((residuals_NEE_complex_Klj)^2)/nrow(residuals_NEE_complex_Klj))
hist(RMSE_NEE_complex_Klj)

# Kormann & Meixner 
pred_neeLT_complex_ZM=t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="neeLT_hat")))[,withheld_index])*test_complex$lichen_tundra_ZM
pred_neeGT_complex_ZM=t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="neeGT_hat")))[,withheld_index])*test_complex$green_tundra_ZM
pred_neeET_complex_ZM=t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="neeET_hat")))[,withheld_index])*test_complex$edge_tundra_ZM
pred_neeED_complex_ZM=t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="neeED_hat")))[,withheld_index])*test_complex$edge_degraded_ZM
pred_neeF1_complex_ZM=t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="neeF1_hat")))[,withheld_index])*(test_complex$fen1_ZM+test_complex$fen2_ZM)
pred_neeD_complex_ZM=t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="neeD_hat")))[,withheld_index])*test_complex$degraded_ZM
pred_water_complex_ZM=matrix(nrow=nrow(pred_neeLT_complex_ZM),ncol=ncol(pred_neeLT_complex_ZM),data=1)
for (i in 1:nrow(pred_water_complex_ZM)){
  pred_water_complex_ZM[i,]=as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="FWprior")))*test_complex$water_ZM[i]
}
pred_NEE_complex_ZM=pred_neeLT_complex_ZM+pred_neeGT_complex_ZM+pred_neeET_complex_ZM+pred_neeED_complex_ZM+pred_neeF1_complex_ZM+pred_neeD_complex_ZM+pred_water_complex_ZM
residuals_NEE_complex_ZM=pred_NEE_complex_ZM-test_complex$co2
residuals_NEE_complex_ZM=na.omit(residuals_NEE_complex_ZM)
RMSE_NEE_complex_ZM=sqrt(colSums((residuals_NEE_complex_ZM)^2)/nrow(residuals_NEE_complex_ZM))
hist(RMSE_NEE_complex_ZM)

# Hsieh 
pred_neeLT_complex_H=t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="neeLT_hat")))[,withheld_index])*test_complex$lichen_tundra_H
pred_neeGT_complex_H=t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="neeGT_hat")))[,withheld_index])*test_complex$green_tundra_H
pred_neeET_complex_H=t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="neeET_hat")))[,withheld_index])*test_complex$edge_tundra_H
pred_neeED_complex_H=t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="neeED_hat")))[,withheld_index])*test_complex$edge_degraded_H
pred_neeF1_complex_H=t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="neeF1_hat")))[,withheld_index])*(test_complex$fen1_H+test_complex$fen2_H)
pred_neeD_complex_H=t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="neeD_hat")))[,withheld_index])*test_complex$degraded_H
pred_water_complex_H=matrix(nrow=nrow(pred_neeLT_complex_H),ncol=ncol(pred_neeLT_complex_H),data=1)
for (i in 1:nrow(pred_water_complex_H)){
  pred_water_complex_H[i,]=as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="FWprior")))*test_complex$water_H[i]
}
pred_NEE_complex_H=pred_neeLT_complex_H+pred_neeGT_complex_H+pred_neeET_complex_H+pred_neeED_complex_H+pred_neeF1_complex_H+pred_neeD_complex_H+pred_water_complex_H
residuals_NEE_complex_H=pred_NEE_complex_H-test_complex$co2
residuals_NEE_complex_H=na.omit(residuals_NEE_complex_H)
RMSE_NEE_complex_H=sqrt(colSums((residuals_NEE_complex_H)^2)/nrow(residuals_NEE_complex_H))
hist(RMSE_NEE_complex_H)

# export complex landcover nee prediction posterior chains

# Kormann and Meixner
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="neeLT_hat"))))),paste0('neeLT_ZM',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="neeGT_hat"))))),paste0('neeGT_ZM',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="neeET_hat"))))),paste0('neeET_ZM',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="neeED_hat"))))),paste0('neeED_ZM',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="neeF1_hat"))))),paste0('neeF1_ZM',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="neeD_hat"))))),paste0('neeD_ZM',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="FWprior"))))),paste0('Fwater_ZM',analysis_date,'_complex_',run_version,'.csv'))

# Hsieh
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="neeLT_hat"))))),paste0('neeLT_H',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="neeGT_hat"))))),paste0('neeGT_H',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="neeET_hat"))))),paste0('neeET_H',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="neeED_hat"))))),paste0('neeED_H',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="neeF1_hat"))))),paste0('neeF1_H',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="neeD_hat"))))),paste0('neeD_H',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="FWprior"))))),paste0('Fwater_H',analysis_date,'_complex_',run_version,'.csv'))

# Kljun
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="neeLT_hat"))))),paste0('neeLT_Klj',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="neeGT_hat"))))),paste0('neeGT_Klj',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="neeET_hat"))))),paste0('neeET_Klj',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="neeED_hat"))))),paste0('neeED_Klj',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="neeF1_hat"))))),paste0('neeF1_Klj',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="neeD_hat"))))),paste0('neeD_Klj',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="FWprior"))))),paste0('Fwater_Klj',analysis_date,'_complex_',run_version,'.csv'))

# export complex parameter posterior chains

# Kormann and Maixner
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="alphaDS"))))),paste0('alphaDS_ZM',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="betaDS"))))),paste0('betaDS_ZM',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="alphaLTS"))))),paste0('alphaLTS_ZM',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="betaLTS"))))),paste0('betaLTS_ZM',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="alphaGTS"))))),paste0('alphaGTS_ZM',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="betaGTS"))))),paste0('betaGTS_ZM',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="alphaETS"))))),paste0('alphaETS_ZM',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="betaETS"))))),paste0('betaETS_ZM',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="alphaFS1"))))),paste0('alphaFS1_ZM',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="betaFS1"))))),paste0('betaFS1_ZM',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="Q10LT"))))),paste0('Q10LT_ZM',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="Q10GT"))))),paste0('Q10GT_ZM',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="Q10ET"))))),paste0('Q10ET_ZM',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="Q10ED"))))),paste0('Q10ED_ZM',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="Q10F1"))))),paste0('Q10F1_ZM',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="Q10D"))))),paste0('Q10D_ZM',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="betaEDS"))))),paste0('betaEDS_ZM',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="alphaEDS"))))),paste0('alphaEDS_ZM',analysis_date,'_complex_',run_version,'.csv'))

write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="pmaxLT"))))),paste0('pmaxLT_ZM',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="pmaxGT"))))),paste0('pmaxGT_ZM',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="pmaxET"))))),paste0('pmaxET_ZM',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="pmaxED"))))),paste0('pmaxED_ZM',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="pmaxF1"))))),paste0('pmaxF1_ZM',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="pmaxDS"))))),paste0('pmaxDS_ZM',analysis_date,'_complex_',run_version,'.csv'))

write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="E0LT"))))),paste0('E0LT_ZM',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="E0GT"))))),paste0('E0GT_ZM',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="E0ET"))))),paste0('E0ET_ZM',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="E0ED"))))),paste0('E0ED_ZM',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="E0F1"))))),paste0('E0F1_ZM',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_ZM,vars="E0DS"))))),paste0('E0DS_ZM',analysis_date,'_complex_',run_version,'.csv'))

# Hsieh
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="alphaDS"))))),paste0('alphaDS_H',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="betaDS"))))),paste0('betaDS_H',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="alphaLTS"))))),paste0('alphaLTS_H',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="betaLTS"))))),paste0('betaLTS_H',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="alphaGTS"))))),paste0('alphaGTS_H',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="betaGTS"))))),paste0('betaGTS_H',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="alphaETS"))))),paste0('alphaETS_H',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="betaETS"))))),paste0('betaETS_H',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="alphaFS1"))))),paste0('alphaFS1_H',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="betaFS1"))))),paste0('betaFS1_H',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="Q10LT"))))),paste0('Q10LT_H',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="Q10GT"))))),paste0('Q10GT_H',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="Q10ET"))))),paste0('Q10ET_H',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="Q10ED"))))),paste0('Q10ED_H',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="Q10F1"))))),paste0('Q10F1_H',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="Q10D"))))),paste0('Q10D_H',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="betaEDS"))))),paste0('betaEDS_H',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="alphaEDS"))))),paste0('alphaEDS_H',analysis_date,'_complex_',run_version,'.csv'))

write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="pmaxLT"))))),paste0('pmaxLT_H',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="pmaxGT"))))),paste0('pmaxGT_H',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="pmaxET"))))),paste0('pmaxET_H',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="pmaxED"))))),paste0('pmaxED_H',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="pmaxF1"))))),paste0('pmaxF1_H',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="pmaxDS"))))),paste0('pmaxDS_H',analysis_date,'_complex_',run_version,'.csv'))

write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="E0LT"))))),paste0('E0LT_H',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="E0GT"))))),paste0('E0GT_H',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="E0ET"))))),paste0('E0ET_H',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="E0ED"))))),paste0('E0ED_H',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="E0F1"))))),paste0('E0F1_H',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_H,vars="E0DS"))))),paste0('E0DS_H',analysis_date,'_complex_',run_version,'.csv'))

# Kljun
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="alphaDS"))))),paste0('alphaDS_Klj',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="betaDS"))))),paste0('betaDS_Klj',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="alphaLTS"))))),paste0('alphaLTS_Klj',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="betaLTS"))))),paste0('betaLTS_Klj',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="alphaGTS"))))),paste0('alphaGTS_Klj',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="betaGTS"))))),paste0('betaGTS_Klj',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="alphaETS"))))),paste0('alphaETS_Klj',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="betaETS"))))),paste0('betaETS_Klj',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="alphaFS1"))))),paste0('alphaFS1_Klj',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="betaFS1"))))),paste0('betaFS1_Klj',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="Q10LT"))))),paste0('Q10LT_Klj',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="Q10GT"))))),paste0('Q10GT_Klj',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="Q10ET"))))),paste0('Q10ET_Klj',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="Q10ED"))))),paste0('Q10ED_Klj',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="Q10F1"))))),paste0('Q10F1_Klj',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="Q10D"))))),paste0('Q10D_Klj',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="betaEDS"))))),paste0('betaEDS_Klj',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="alphaEDS"))))),paste0('alphaEDS_Klj',analysis_date,'_complex_',run_version,'.csv'))

write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="pmaxLT"))))),paste0('pmaxLT_Klj',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="pmaxGT"))))),paste0('pmaxGT_Klj',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="pmaxET"))))),paste0('pmaxET_Klj',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="pmaxED"))))),paste0('pmaxED_Klj',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="pmaxF1"))))),paste0('pmaxF1_Klj',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="pmaxDS"))))),paste0('pmaxDS_Klj',analysis_date,'_complex_',run_version,'.csv'))

write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="E0LT"))))),paste0('E0LT_Klj',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="E0GT"))))),paste0('E0GT_Klj',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="E0ET"))))),paste0('E0ET_Klj',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="E0ED"))))),paste0('E0ED_Klj',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="E0F1"))))),paste0('E0F1_Klj',analysis_date,'_complex_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_complex_Klj,vars="E0DS"))))),paste0('E0DS_Klj',analysis_date,'_complex_',run_version,'.csv'))





#### Simple NEE ####

# 1. set up initial values, change as needed
initsNEE_simple=list(list(sigma_fit=1,alphaDS=1*0.9,betaDS=0.08*0.85,FWprior=0.5*0.9,E0DS=0.005,pmaxDS=5.90,
                          alphaLTS=1*0.9,betaLTS=0.03*0.85,E0LT=0.005,pmaxF1=5.90*0.9,
                          alphaFS1=1*0.9,betaFS1=0.1*0.85,E0F1=0.005,pmaxLT=5.90*0.9),
                     list(sigma_fit=2,alphaDS=1,betaDS=0.08,FWprior=2,E0DS=0.01,pmaxDS=15,
                          alphaLTS=1,betaLTS=0.03,E0LT=0.03038,pmaxF1=14.33,
                          alphaFS1=1,betaFS1=0.09,E0F1=0.013045,pmaxLT=14.33),
                     list(sigma_fit=2.5,alphaDS=1*1.1,betaDS=0.08*2,FWprior=5*1.1,E0DS=0.13045,pmaxDS=25,
                          alphaLTS=1*1.1,betaLTS=0.03*1.3,E0LT=0.13045,pmaxF1=25*1.1,
                          alphaFS1=1*1.1,betaFS1=0.08*1.2,E0F1=0.13045,pmaxLT=25*1.1))

initsRdark_simple=list(list(sigma_fit=1,alphaDS=0.1*0.9,betaDS=0.01*0.9,FWprior=0.2*0.9,
                            alphaLTS=0.1*0.9,betaLTS=0.01*0.9,
                            alphaFS1=0.1*0.9,betaFS1=0.01*0.9),
                       list(sigma_fit=2,alphaDS=1,betaDS=0.04,FWprior=2,
                            alphaLTS=1,betaLTS=0.04,
                            alphaFS1=1,betaFS1=0.04),
                       list(sigma_fit=2.5,alphaDS=5*1.1,betaDS=0.18*1.1,FWprior=12*1.1,
                            alphaLTS=5*1.1,betaLTS=0.18*1.1,
                            alphaFS1=5*1.1,betaFS1=0.18*1.1))

# 2. set up data lists for jags model
Tscale=((train_simple$air_temperature-Tmin)*(train_simple$air_temperature-Tmax))/((train_simple$air_temperature-Tmin)*(train_simple$air_temperature-Tmax)-(train_simple$air_temperature-Topt)^2)
data_simple_ZM=list(Tsoil=train_simple$air_temperature,PAR=train_simple$PPFD_1_1_1,
                    co2=train_simple$co2,ftp_LT=train_simple$tundra_ZM,
                    ftp_F1=train_simple$fen1_ZM,ftp_D=train_simple$degraded_ZM,
                    ftp_W=train_simple$water_ZM,Tscale=Tscale,T_hat=df_hat$Tair_hat,
                    Tscale_hat=df_hat$Tscale_hat,PAR_hat=df_hat$PAR_hat)
data_simple_H=list(Tsoil=train_simple$air_temperature,PAR=train_simple$PPFD_1_1_1,
                   co2=train_simple$co2,ftp_LT=train_simple$tundra_H,
                   ftp_F1=train_simple$fen1_H,ftp_D=train_simple$degraded_H,
                   ftp_W=train_simple$water_H,Tscale=Tscale,T_hat=df_hat$Tair_hat,
                   Tscale_hat=df_hat$Tscale_hat,PAR_hat=df_hat$PAR_hat)
data_simple_Klj=list(Tsoil=train_simple$air_temperature,PAR=train_simple$PPFD_1_1_1,
                     co2=train_simple$co2,ftp_LT=train_simple$tundra_Klj,
                     ftp_F1=train_simple$fen1_Klj,ftp_D=train_simple$degraded_Klj,
                     ftp_W=train_simple$water_Klj,Tscale=Tscale,T_hat=df_hat$Tair_hat,
                     Tscale_hat=df_hat$Tscale_hat,PAR_hat=df_hat$PAR_hat)

data_simple_dark_ZM=list(Tsoil=df_dark_simple$air_temperature,
                         co2=df_dark_simple$co2,ftp_LT=df_dark_simple$tundra_ZM,ftp_F1=df_dark_simple$fen1_ZM,ftp_D=df_dark_simple$degraded_ZM,
                         ftp_W=df_dark_simple$water_ZM,T_hat=df_hat$Tair_hat)
data_simple_dark_H=list(Tsoil=df_dark_simple$air_temperature,
                        co2=df_dark_simple$co2,ftp_LT=df_dark_simple$tundra_H,ftp_F1=df_dark_simple$fen1_H,ftp_D=df_dark_simple$degraded_H,
                        ftp_W=df_dark_simple$water_H,T_hat=df_hat$Tair_hat)
data_simple_dark_Klj=list(Tsoil=df_dark_simple$air_temperature,
                          co2=df_dark_simple$co2,ftp_LT=df_dark_simple$tundra_Klj,ftp_F1=df_dark_simple$fen1_Klj,ftp_D=df_dark_simple$degraded_Klj,
                          ftp_W=df_dark_simple$water_Klj,T_hat=df_hat$Tair_hat)

# 3. jags run parameters
n.adapt=5000
n.update=5000
n.iter=3000

# 4. dark fit for respiration
set.seed(1)
cl=makeCluster(3)
rjmd_simple_ZM=run.jags(model="Simple_Resp_dark_JAGS.R",data=data_simple_dark_ZM,
                        monitor=c('sigma_fit','alphaDS','betaDS','Q10D',
                                  'alphaLTS','betaLTS','alphaFS1','Q10LT','FWprior',
                                  'betaFS1','Q10F1','Fwater',
                                  'p.fit','p.mu','p.sd',
                                  'RsLT_hat','RsF1_hat', 'RsD_hat'),
                        method='rjparallel',summarise = FALSE,cl=cl,
                        n.chains=length(initsRdark_simple),inits=initsRdark_simple,burnin=n.update, adapt=n.adapt,sample=n.iter)
stopCluster(cl)

set.seed(1)
cl=makeCluster(3)
rjmd_simple_H=run.jags(model="Simple_Resp_dark_JAGS.R",data=data_simple_dark_H,
                       monitor=c('sigma_fit','alphaDS','betaDS','Q10D',
                                 'alphaLTS','betaLTS','alphaFS1','Q10LT','FWprior',
                                 'betaFS1','Q10F1','Fwater',
                                 'p.fit','p.mu','p.sd',
                                 'RsLT_hat','RsF1_hat', 'RsD_hat'),
                       method='rjparallel',summarise = FALSE,cl=cl,
                       n.chains=length(initsRdark_simple),inits=initsRdark_simple,burnin=n.update, adapt=n.adapt,sample=n.iter)
stopCluster(cl)

set.seed(1)
cl=makeCluster(3)
rjmd_simple_Klj=run.jags(model="Simple_Resp_dark_JAGS.R",data=data_simple_dark_Klj,
                         monitor=c('sigma_fit','alphaDS','betaDS','Q10D',
                                   'alphaLTS','betaLTS','alphaFS1','Q10LT','FWprior',
                                   'betaFS1','Q10F1','Fwater',
                                   'p.fit','p.mu','p.sd',
                                   'RsLT_hat','RsF1_hat', 'RsD_hat'),
                         method='rjparallel',summarise = FALSE,cl=cl,
                         n.chains=length(initsRdark_simple),inits=initsRdark_simple,burnin=n.update, adapt=n.adapt,sample=n.iter)
stopCluster(cl)

# check model fit params:
# should be > 0.05 and < 0.95
mean(as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_ZM,vars="p.fit"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_ZM,vars="p.mu"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_ZM,vars="p.sd"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_H,vars="p.fit"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_H,vars="p.mu"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_H,vars="p.sd"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_Klj,vars="p.fit"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_Klj,vars="p.mu"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_Klj,vars="p.sd"))))

#Convergence tests:
flux_paramsd_simple=c('sigma_fit','alphaDS','betaDS',
                      'alphaLTS','betaLTS','alphaFS1','FWprior',
                      'betaFS1')

# plot(rjmd_simple_ZM,vars=flux_paramsd_simple)
# plot(rjmd_simple_H,vars=flux_paramsd_simple)
# plot(rjmd_simple_Klj,vars=flux_paramsd_simple)

gelman_simple_listd_ZM=list()
for (i in 1:length(flux_paramsd_simple)){
  param=flux_paramsd_simple[i]
  result=gelman.diag(as.mcmc.list(rjmd_simple_ZM,vars=param))
  gelman_simple_listd_ZM[i]=result$psrf[1]
};gelman_simple_listd_ZM
gelman_simple_listd_H=list()
for (i in 1:length(flux_paramsd_simple)){
  param=flux_paramsd_simple[i]
  result=gelman.diag(as.mcmc.list(rjmd_simple_H,vars=param))
  gelman_simple_listd_H[i]=result$psrf[1]
};gelman_simple_listd_H
gelman_simple_listd_Klj=list()
for (i in 1:length(flux_paramsd_simple)){
  param=flux_paramsd_simple[i]
  result=gelman.diag(as.mcmc.list(rjmd_simple_Klj,vars=param))
  gelman_simple_listd_Klj[i]=result$psrf[1]
};gelman_simple_listd_Klj

# export posterior chains for parameters for dark model

# Parameters for ZM model
alphaDS_simple_ZM=as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_ZM,vars="alphaDS")))
betaDS_simple_ZM=as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_ZM,vars="betaDS")))
alphaLTS_simple_ZM=as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_ZM,vars="alphaLTS")))
betaLTS_simple_ZM=as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_ZM,vars="betaLTS")))
alphaFS1_simple_ZM=as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_ZM,vars="alphaFS1")))
betaFS1_simple_ZM=as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_ZM,vars="betaFS1")))
Q10LT_simple_ZM=as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_ZM,vars="Q10LT")))
Q10F1_simple_ZM=as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_ZM,vars="Q10F1")))
Q10D_simple_ZM=as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_ZM,vars="Q10D")))
FWprior_simple_ZM=as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_ZM,vars="FWprior")))

# Parameters for Hsieh model
alphaDS_simple_H=as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_H,vars="alphaDS")))
betaDS_simple_H=as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_H,vars="betaDS")))
alphaLTS_simple_H=as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_H,vars="alphaLTS")))
betaLTS_simple_H=as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_H,vars="betaLTS")))
alphaFS1_simple_H=as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_H,vars="alphaFS1")))
betaFS1_simple_H=as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_H,vars="betaFS1")))
Q10LT_simple_H=as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_H,vars="Q10LT")))
Q10F1_simple_H=as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_H,vars="Q10F1")))
Q10D_simple_H=as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_H,vars="Q10D")))
FWprior_simple_H=as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_H,vars="FWprior")))

# Parameters for Kljun model
alphaDS_simple_Klj=as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_Klj,vars="alphaDS")))
betaDS_simple_Klj=as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_Klj,vars="betaDS")))
alphaLTS_simple_Klj=as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_Klj,vars="alphaLTS")))
betaLTS_simple_Klj=as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_Klj,vars="betaLTS")))
alphaFS1_simple_Klj=as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_Klj,vars="alphaFS1")))
betaFS1_simple_Klj=as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_Klj,vars="betaFS1")))
Q10LT_simple_Klj=as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_Klj,vars="Q10LT")))
Q10F1_simple_Klj=as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_Klj,vars="Q10F1")))
Q10D_simple_Klj=as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_Klj,vars="Q10D")))
FWprior_simple_Klj=as.matrix(combine.mcmc(as.mcmc.list(rjmd_simple_Klj,vars="FWprior")))

# Hsieh Respiration
params_list=list(Q10D_simple_H,Q10F1_simple_H,Q10LT_simple_H,
                 betaFS1_simple_H,alphaFS1_simple_H,betaDS_simple_H,alphaDS_simple_H,
                 betaLTS_simple_H,alphaLTS_simple_H,FWprior_simple_H)
params_export=matrix(nrow=length(params_list),ncol=5)
for (i in 1:length(params_list)){
  iter=params_list[[i]]
  params_export[i,1]=attributes(iter)$dimnames[[2]]
  param_stats=quantile(iter,probs=c(0.025,0.5,0.975))
  params_export[i,2]=param_stats[2]
  params_export[i,3]=mean(iter)
  params_export[i,4]=param_stats[1]
  params_export[i,5]=param_stats[3]
}
params_export=as.data.frame(params_export)
names(params_export)=c('parameter','median_est','mean_est','lower2_5','upper97_5')
write_csv(params_export,paste0('Hsieh_Rd_params',analysis_date,'_simple_',run_version,'.csv'))

# KM Respiration
params_list=list(Q10D_simple_ZM,Q10F1_simple_ZM,Q10LT_simple_ZM,
                 betaFS1_simple_ZM,alphaFS1_simple_ZM,betaDS_simple_ZM,alphaDS_simple_ZM,
                 betaLTS_simple_ZM,alphaLTS_simple_ZM,FWprior_simple_ZM)
params_export=matrix(nrow=length(params_list),ncol=5)
for (i in 1:length(params_list)){
  iter=params_list[[i]]
  params_export[i,1]=attributes(iter)$dimnames[[2]]
  param_stats=quantile(iter,probs=c(0.025,0.5,0.975))
  params_export[i,2]=param_stats[2]
  params_export[i,3]=mean(iter)
  params_export[i,4]=param_stats[1]
  params_export[i,5]=param_stats[3]
}
params_export=as.data.frame(params_export)
names(params_export)=c('parameter','median_est','mean_est','lower2_5','upper97_5')
write_csv(params_export,paste0('KormannMeixner_Rd_params',analysis_date,'_simple_',run_version,'.csv'))

# KLjun Respiration
params_list=list(Q10D_simple_Klj,Q10F1_simple_Klj,Q10LT_simple_Klj,
                 betaFS1_simple_Klj,alphaFS1_simple_Klj,betaDS_simple_Klj,alphaDS_simple_Klj,
                 betaLTS_simple_Klj,alphaLTS_simple_Klj,FWprior_simple_Klj)

params_export=matrix(nrow=length(params_list),ncol=5)
for (i in 1:length(params_list)){
  iter=params_list[[i]]
  params_export[i,1]=attributes(iter)$dimnames[[2]]
  param_stats=quantile(iter,probs=c(0.025,0.5,0.975))
  params_export[i,2]=param_stats[2]
  params_export[i,3]=mean(iter)
  params_export[i,4]=param_stats[1]
  params_export[i,5]=param_stats[3]
}
params_export=as.data.frame(params_export)
names(params_export)=c('parameter','median_est','mean_est','lower2_5','upper97_5')
write_csv(params_export,paste0('Kljun_Rd_params',analysis_date,'_simple_',run_version,'.csv'))


# 5.  all light fit for NEE

full_var_list_simple=c('sigma_fit','Q10D', 'Q10LT',
                       'Q10F1','p.fit','p.mu','p.sd',
                       'pmaxLT','pmaxF1','FWprior',
                       'E0LT','E0F1','E0DS','pmaxDS',
                       'alphaDS','betaDS', 'alphaFS1','betaFS1','alphaLTS','betaLTS',
                       'gppLT_hat','gppF1_hat','gppD_hat',
                       'neeLT_hat','neeF1_hat','RsD_hat',
                       'RsLT_hat','RsF1_hat','neeD_hat')

set.seed(1)
cl=makeCluster(3)
rjm_simple_ZM=run.jags(model=paste0("Simple_",JAGS_model,"_KM.R"),data=data_simple_ZM,
                       monitor=full_var_list4,
                       method='rjparallel',summarise = FALSE,cl=cl,
                       n.chains=length(initsNEE_simple),inits=initsNEE_simple,burnin=n.update, adapt=n.adapt,sample=n.iter)
stopCluster(cl)

set.seed(1)
cl=makeCluster(3)
rjm_simple_H=run.jags(model=paste0("Simple_",JAGS_model,"_H.R"),data=data_simple_H,
                      monitor=full_var_list4,
                      method='rjparallel',summarise = FALSE,cl=cl,
                      n.chains=length(initsNEE_simple),inits=initsNEE_simple,burnin=n.update, adapt=n.adapt,sample=n.iter)
stopCluster(cl)

set.seed(1)
cl=makeCluster(3)
rjm_simple_Klj=run.jags(model=paste0("Simple_",JAGS_model,"_Klj.R"),data=data_simple_Klj,
                        monitor=full_var_list4,
                        method='rjparallel',summarise = FALSE,cl=cl,
                        n.chains=length(initsNEE_simple),inits=initsNEE_simple,burnin=n.update, adapt=n.adapt,sample=n.iter)
stopCluster(cl)

# check model fit params:
# should be > 0.05 and < 0.95
mean(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_ZM,vars="p.fit"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_ZM,vars="p.mu"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_ZM,vars="p.sd"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_H,vars="p.fit"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_H,vars="p.mu"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_H,vars="p.sd"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_Klj,vars="p.fit"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_Klj,vars="p.mu"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_Klj,vars="p.sd"))))

# Convergence tests
flux_params_simple=c( 'sigma_fit',
                      'FWprior',#'E0DS',
                      'alphaDS','betaDS', 'alphaFS1','betaFS1','alphaLTS','betaLTS',
                      'pmaxLT','pmaxF1',#'pmaxDS',
                      'E0LT','E0F1')
# uncomment pmaxDS and E0DS for 0520, 0620, 07all runs

# plot(rjm_simple_ZM,vars=flux_params_simple)
# plot(rjm_simple_H,vars=flux_params_simple)
# plot(rjm_simple_Klj,vars=flux_params_simple)

gelman_simple_list_ZM=list()
for (i in 1:length(flux_params_simple)){
  param=flux_params_simple[i]
  result=gelman.diag(as.mcmc.list(rjm_simple_ZM,vars=param))
  gelman_simple_list_ZM[i]=result$psrf[1]
};gelman_simple_list_ZM

gelman_simple_list_H=list()
for (i in 1:length(flux_params_simple)){
  param=flux_params_simple[i]
  result=gelman.diag(as.mcmc.list(rjm_simple_H,vars=param))
  gelman_simple_list_H[i]=result$psrf[1]
};gelman_simple_list_H
gelman_simple_list_Klj=list()
for (i in 1:length(flux_params_simple)){
  param=flux_params_simple[i]
  result=gelman.diag(as.mcmc.list(rjm_simple_Klj,vars=param))
  gelman_simple_list_Klj[i]=result$psrf[1]
};gelman_simple_list_Klj

# compare to withheld data

# Kljun 
pred_neeLT_simple_Klj=t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_Klj,vars="neeLT_hat")))[,withheld_index])*test_simple$tundra_Klj
pred_neeF1_simple_Klj=t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_Klj,vars="neeF1_hat")))[,withheld_index])*(test_simple$fen1_Klj)
pred_neeD_simple_Klj=t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_Klj,vars="neeD_hat")))[,withheld_index])*test_simple$degraded_Klj
pred_water_simple_Klj=matrix(nrow=nrow(pred_neeLT_simple_Klj),ncol=ncol(pred_neeLT_simple_Klj),data=1)
for (i in 1:nrow(pred_water_simple_Klj)){
  pred_water_simple_Klj[i,]=as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_Klj,vars="FWprior")))*test_simple$water_Klj[i]
}
pred_NEE_simple_Klj=pred_neeLT_simple_Klj+pred_neeF1_simple_Klj+pred_neeD_simple_Klj+pred_water_simple_Klj
residuals_NEE_simple_Klj=pred_NEE_simple_Klj-test_simple$co2
residuals_NEE_simple_Klj=na.omit(residuals_NEE_simple_Klj)
RMSE_NEE_simple_Klj=sqrt(colSums((residuals_NEE_simple_Klj)^2)/nrow(residuals_NEE_simple_Klj))
hist(RMSE_NEE_simple_Klj)

# Kormann & Meixner 
pred_neeLT_simple_ZM=t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_ZM,vars="neeLT_hat")))[,withheld_index])*test_simple$tundra_ZM
pred_neeF1_simple_ZM=t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_ZM,vars="neeF1_hat")))[,withheld_index])*(test_simple$fen1_ZM)
pred_neeD_simple_ZM=t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_ZM,vars="neeD_hat")))[,withheld_index])*test_simple$degraded_ZM
pred_water_simple_ZM=matrix(nrow=nrow(pred_neeD_simple_ZM),ncol=ncol(pred_neeD_simple_ZM),data=1)
for (i in 1:nrow(pred_water_simple_ZM)){
  pred_water_simple_ZM[i,]=as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_ZM,vars="FWprior")))*test_simple$water_ZM[i]
}
pred_NEE_simple_ZM=pred_neeLT_simple_ZM+pred_neeF1_simple_ZM+pred_neeD_simple_ZM+pred_water_simple_ZM
residuals_NEE_simple_ZM=pred_NEE_simple_ZM-test_simple$co2
residuals_NEE_simple_ZM=na.omit(residuals_NEE_simple_ZM)
RMSE_NEE_simple_ZM=sqrt(colSums((residuals_NEE_simple_ZM)^2)/nrow(residuals_NEE_simple_ZM))
hist(RMSE_NEE_simple_ZM)

# Hsieh 
pred_neeLT_simple_H=t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_H,vars="neeLT_hat")))[,withheld_index])*test_simple$tundra_H
pred_neeF1_simple_H=t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_H,vars="neeF1_hat")))[,withheld_index])*(test_simple$fen1_H)
pred_neeD_simple_H=t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_H,vars="neeD_hat")))[,withheld_index])*test_simple$degraded_H
pred_water_simple_H=matrix(nrow=nrow(pred_neeD_simple_H),ncol=ncol(pred_neeD_simple_H),data=1)
for (i in 1:nrow(pred_water_simple_H)){
  pred_water_simple_H[i,]=as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_H,vars="FWprior")))*test_simple$water_H[i]
}
pred_NEE_simple_H=pred_neeLT_simple_H+pred_neeF1_simple_H+pred_neeD_simple_H+pred_water_simple_H
residuals_NEE_simple_H=pred_NEE_simple_H-test_simple$co2
residuals_NEE_simple_H=na.omit(residuals_NEE_simple_H)
RMSE_NEE_simple_H=sqrt(colSums((residuals_NEE_simple_H)^2)/nrow(residuals_NEE_simple_H))
hist(RMSE_NEE_simple_H)


# export posterior chains for landcover nee predictions

# Kormann & Meixner
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_ZM,vars="neeLT_hat"))))),paste0('neeLT_ZM',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_ZM,vars="neeF1_hat"))))),paste0('neeF1_ZM',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_ZM,vars="neeD_hat"))))),paste0('neeD_ZM',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_ZM,vars="FWprior"))))),paste0('Fwater_ZM',analysis_date,'_simple_',run_version,'.csv'))

# Hsieh
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_H,vars="neeLT_hat"))))),paste0('neeLT_H',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_H,vars="neeF1_hat"))))),paste0('neeF1_H',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_H,vars="neeD_hat"))))),paste0('neeD_H',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_H,vars="FWprior"))))),paste0('Fwater_H',analysis_date,'_simple_',run_version,'.csv'))

# Kljun
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_Klj,vars="neeLT_hat"))))),paste0('neeLT_Klj',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_Klj,vars="neeF1_hat"))))),paste0('neeF1_Klj',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_Klj,vars="neeD_hat"))))),paste0('neeD_Klj',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_Klj,vars="FWprior"))))),paste0('Fwater_Klj',analysis_date,'_simple_',run_version,'.csv'))

# export posterior chains for parameters for NEE model

# Kormann & Meixner
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_ZM,vars="alphaDS"))))),paste0('alphaDS_ZM',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_ZM,vars="betaDS"))))),paste0('betaDS_ZM',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_ZM,vars="alphaLTS"))))),paste0('alphaLTS_ZM',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_ZM,vars="betaLTS"))))),paste0('betaLTS_ZM',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_ZM,vars="alphaFS1"))))),paste0('alphaFS1_ZM',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_ZM,vars="betaFS1"))))),paste0('betaFS1_ZM',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_ZM,vars="Q10LT"))))),paste0('Q10LT_ZM',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_ZM,vars="Q10F1"))))),paste0('Q10F1_ZM',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_ZM,vars="Q10D"))))),paste0('Q10D_ZM',analysis_date,'_simple_',run_version,'.csv'))

write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_ZM,vars="pmaxLT"))))),paste0('pmaxLT_ZM',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_ZM,vars="pmaxF1"))))),paste0('pmaxF1_ZM',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_ZM,vars="pmaxDS"))))),paste0('pmaxDS_ZM',analysis_date,'_simple_',run_version,'.csv'))

write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_ZM,vars="E0LT"))))),paste0('E0LT_ZM',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_ZM,vars="E0F1"))))),paste0('E0F1_ZM',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_ZM,vars="E0DS"))))),paste0('E0DS_ZM',analysis_date,'_simple_',run_version,'.csv'))

# Hsieh
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_H,vars="alphaDS"))))),paste0('alphaDS_H',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_H,vars="betaDS"))))),paste0('betaDS_H',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_H,vars="alphaLTS"))))),paste0('alphaLTS_H',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_H,vars="betaLTS"))))),paste0('betaLTS_H',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_H,vars="alphaFS1"))))),paste0('alphaFS1_H',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_H,vars="betaFS1"))))),paste0('betaFS1_H',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_H,vars="Q10LT"))))),paste0('Q10LT_H',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_H,vars="Q10F1"))))),paste0('Q10F1_H',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_H,vars="Q10D"))))),paste0('Q10D_H',analysis_date,'_simple_',run_version,'.csv'))

write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_H,vars="pmaxLT"))))),paste0('pmaxLT_H',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_H,vars="pmaxF1"))))),paste0('pmaxF1_H',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_H,vars="pmaxDS"))))),paste0('pmaxDS_H',analysis_date,'_simple_',run_version,'.csv'))

write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_H,vars="E0LT"))))),paste0('E0LT_H',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_H,vars="E0F1"))))),paste0('E0F1_H',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_H,vars="E0DS"))))),paste0('E0DS_H',analysis_date,'_simple_',run_version,'.csv'))

# Kljun
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_Klj,vars="alphaDS"))))),paste0('alphaDS_Klj',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_Klj,vars="betaDS"))))),paste0('betaDS_Klj',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_Klj,vars="alphaLTS"))))),paste0('alphaLTS_Klj',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_Klj,vars="betaLTS"))))),paste0('betaLTS_Klj',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_Klj,vars="alphaFS1"))))),paste0('alphaFS1_Klj',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_Klj,vars="betaFS1"))))),paste0('betaFS1_Klj',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_Klj,vars="Q10LT"))))),paste0('Q10LT_Klj',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_Klj,vars="Q10F1"))))),paste0('Q10F1_Klj',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_Klj,vars="Q10D"))))),paste0('Q10D_Klj',analysis_date,'_simple_',run_version,'.csv'))

write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_Klj,vars="pmaxLT"))))),paste0('pmaxLT_Klj',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_Klj,vars="pmaxF1"))))),paste0('pmaxF1_Klj',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_Klj,vars="pmaxDS"))))),paste0('pmaxDS_Klj',analysis_date,'_simple_',run_version,'.csv'))

write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_Klj,vars="E0LT"))))),paste0('E0LT_Klj',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_Klj,vars="E0F1"))))),paste0('E0F1_Klj',analysis_date,'_simple_',run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_simple_Klj,vars="E0DS"))))),paste0('E0DS_Klj',analysis_date,'_simple_',run_version,'.csv'))


#### CH4  ####

# 1. set up initial values, change as needed
initsCH4=list(list(sigma_fit=0.01,ch4DS=0.01*0.9,
                     ch4LTS=0.001*0.9,
                     ch4GTS=0.001*0.9,
                     ch4ETS=0.001*0.9,
                     ch4FS1=0.01*0.9,
                     ch4EDS=0.01*0.9,ch4FW=0.01*0.9),
                 list(sigma_fit=0.01,ch4DS=0.01,
                      ch4LTS=0.001,
                      ch4GTS=0.001,
                      ch4ETS=0.001,
                      ch4FS1=0.01,
                      ch4EDS=0.01,ch4FW=0.01),
                 list(sigma_fit=0.01,ch4DS=0.01*2,
                      ch4LTS=0.001*2,
                      ch4GTS=0.001*2,
                      ch4ETS=0.001*2,
                      ch4FS1=0.01*2,
                      ch4EDS=0.01*2,ch4FW=0.01*2))


# 2. set up data lists for jags model
train_bays_ch4=subset(train_complex,!is.na(train_complex$ch4))

data_ch4_jags_ZM=list(ch4=train_bays_ch4$ch4,ftp_LT=train_bays_ch4$lichen_tundra_ZM,ftp_GT=train_bays_ch4$green_tundra_ZM,ftp_ET=train_bays_ch4$edge_tundra_ZM,
                   ftp_F1=(train_bays_ch4$fen1_ZM+train_bays_ch4$fen2_ZM),ftp_D=train_bays_ch4$degraded_ZM,
                   ftp_W=train_bays_ch4$water_ZM,ftp_ED=train_bays_ch4$edge_degraded_ZM)
data_ch4_jags_H=list(ch4=train_bays_ch4$ch4,ftp_LT=train_bays_ch4$lichen_tundra_H,ftp_GT=train_bays_ch4$green_tundra_H,ftp_ET=train_bays_ch4$edge_tundra_H,
                  ftp_F1=(train_bays_ch4$fen1_H+train_bays_ch4$fen2_H),ftp_D=train_bays_ch4$degraded_H,
                  ftp_W=train_bays_ch4$water_H,ftp_ED=train_bays_ch4$edge_degraded_H)
data_ch4_jags_Klj=list(ch4=train_bays_ch4$ch4,ftp_LT=train_bays_ch4$lichen_tundra_Klj,ftp_GT=train_bays_ch4$green_tundra_Klj,ftp_ET=train_bays_ch4$edge_tundra_Klj,
                    ftp_F1=(train_bays_ch4$fen1_Klj+train_bays_ch4$fen2_Klj),ftp_D=train_bays_ch4$degraded_Klj,
                    ftp_W=train_bays_ch4$water_Klj,ftp_ED=train_bays_ch4$edge_degraded_Klj)

# 3. jags run parameters
n.adapt=5000
n.update=5000
n.iter=3000

cat_ch4_params=c('sigma_fit','ch4DS','ch4FW',
                 'ch4LTS','ch4GTS','ch4ETS',
                  'ch4FS1','ch4EDS')

# 4. run jags models for ch4
set.seed(1)
cl=makeCluster(3)
rjm_ch4_ZM=run.jags(model="CH4_JAGS.R",data=data_ch4_jags_ZM,
                  monitor=c('sigma_fit','ch4DS','ch4FW',
                            'ch4LTS','ch4FS1','ch4EDS',
                            'ch4GTS','ch4ETS',
                            'p.fit','p.mu','p.sd'),
                  method='rjparallel',summarise = FALSE,cl=cl,
                  n.chains=length(initsCH4),inits=initsCH4,burnin=n.update, adapt=n.adapt,sample=n.iter)
stopCluster(cl)

set.seed(1)
cl=makeCluster(3)
rjm_ch4_H=run.jags(model="CH4_JAGS.R",data=data_ch4_jags_H,
                     monitor=c('sigma_fit','ch4DS','ch4FW',
                               'ch4LTS','ch4FS1','ch4EDS',
                               'ch4GTS','ch4ETS',
                               'p.fit','p.mu','p.sd'),
                     method='rjparallel',summarise = FALSE,cl=cl,
                     n.chains=length(initsCH4),inits=initsCH4,burnin=n.update, adapt=n.adapt,sample=n.iter)
stopCluster(cl)

set.seed(1)
cl=makeCluster(3)
rjm_ch4_Klj=run.jags(model="CH4_JAGS.R",data=data_ch4_jags_Klj,
                     monitor=c('sigma_fit','ch4DS','ch4FW',
                               'ch4LTS','ch4FS1','ch4EDS',
                               'ch4GTS','ch4ETS',
                               'p.fit','p.mu','p.sd'),
                     method='rjparallel',summarise = FALSE,cl=cl,
                     n.chains=length(initsCH4),inits=initsCH4,burnin=n.update, adapt=n.adapt,sample=n.iter)
stopCluster(cl)

# check model fit params:
# should be > 0.05 and < 0.95
mean(as.matrix(combine.mcmc(as.mcmc.list(rjm_ch4_ZM,vars="p.fit"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjm_ch4_ZM,vars="p.mu"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjm_ch4_ZM,vars="p.sd"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjm_ch4_H,vars="p.fit"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjm_ch4_H,vars="p.mu"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjm_ch4_H,vars="p.sd"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjm_ch4_Klj,vars="p.fit"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjm_ch4_Klj,vars="p.mu"))))
mean(as.matrix(combine.mcmc(as.mcmc.list(rjm_ch4_Klj,vars="p.sd"))))

# Convergence tests:

gelman_list_ZM_ch4=list()
for (i in 1:length(cat_ch4_params)){
  param=cat_ch4_params[i]
  result=gelman.diag(as.mcmc.list(rjm_ch4_ZM,vars=param))
  gelman_list_ZM_ch4[i]=result$psrf[1]
};gelman_list_ZM_ch4

gelman_list_H_ch4=list()
for (i in 1:length(cat_ch4_params)){
  param=cat_ch4_params[i]
  result=gelman.diag(as.mcmc.list(rjm_ch4_H,vars=param))
  gelman_list_H_ch4[i]=result$psrf[1]
};gelman_list_H_ch4

gelman_list_Klj_ch4=list()
for (i in 1:length(cat_ch4_params)){
  param=cat_ch4_params[i]
  result=gelman.diag(as.mcmc.list(rjm_ch4_Klj,vars=param))
  gelman_list_Klj_ch4[i]=result$psrf[1]
};gelman_list_Klj_ch4


# compare to withheld data

# Klun
pred_water_ch4_Klj=t(replicate(n=length(test_complex$ch4),as.vector(combine.mcmc(as.mcmc.list(rjm_ch4_Klj,vars="ch4FW")))))*test_complex$water_Klj
pred_fen_ch4_Klj=t(replicate(n=length(test_complex$ch4),as.vector(combine.mcmc(as.mcmc.list(rjm_ch4_Klj,vars="ch4FS1")))))*(test_complex$fen1_Klj+test_complex$fen2_Klj)
pred_deg_ch4_Klj=t(replicate(n=length(test_complex$ch4),as.vector(combine.mcmc(as.mcmc.list(rjm_ch4_Klj,vars="ch4DS")))))*test_complex$degraded_Klj
pred_degE_ch4_Klj=t(replicate(n=length(test_complex$ch4),as.vector(combine.mcmc(as.mcmc.list(rjm_ch4_Klj,vars="ch4EDS")))))*test_complex$edge_degraded_Klj
pred_LT_ch4_Klj=t(replicate(n=length(test_complex$ch4),as.vector(combine.mcmc(as.mcmc.list(rjm_ch4_Klj,vars="ch4LTS")))))*test_complex$lichen_tundra_Klj
pred_ET_ch4_Klj=t(replicate(n=length(test_complex$ch4),as.vector(combine.mcmc(as.mcmc.list(rjm_ch4_Klj,vars="ch4ETS")))))*test_complex$edge_tundra_Klj
pred_GT_ch4_Klj=t(replicate(n=length(test_complex$ch4),as.vector(combine.mcmc(as.mcmc.list(rjm_ch4_Klj,vars="ch4GTS")))))*test_complex$green_tundra_Klj

pred_ch4_Klj=pred_water_ch4_Klj+pred_fen_ch4_Klj+pred_deg_ch4_Klj+pred_degE_ch4_Klj+pred_LT_ch4_Klj+pred_ET_ch4_Klj+pred_GT_ch4_Klj
residuals_ch4_Klj=pred_ch4_Klj-test_complex$ch4
residuals_ch4_Klj=na.omit(residuals_ch4_Klj)
RMSE_ch4_Klj=sqrt(colSums((residuals_ch4_Klj)^2)/nrow(residuals_ch4_Klj))
hist(RMSE_ch4_Klj)

# Hsieh
pred_water_ch4_H=t(replicate(n=length(test_complex$ch4),as.vector(combine.mcmc(as.mcmc.list(rjm_ch4_H,vars="ch4FW")))))*test_complex$water_H
pred_fen_ch4_H=t(replicate(n=length(test_complex$ch4),as.vector(combine.mcmc(as.mcmc.list(rjm_ch4_H,vars="ch4FS1")))))*(test_complex$fen1_H+test_complex$fen2_H)
pred_deg_ch4_H=t(replicate(n=length(test_complex$ch4),as.vector(combine.mcmc(as.mcmc.list(rjm_ch4_H,vars="ch4DS")))))*test_complex$degraded_H
pred_degE_ch4_H=t(replicate(n=length(test_complex$ch4),as.vector(combine.mcmc(as.mcmc.list(rjm_ch4_H,vars="ch4EDS")))))*test_complex$edge_degraded_H
pred_LT_ch4_H=t(replicate(n=length(test_complex$ch4),as.vector(combine.mcmc(as.mcmc.list(rjm_ch4_H,vars="ch4LTS")))))*test_complex$lichen_tundra_H
pred_ET_ch4_H=t(replicate(n=length(test_complex$ch4),as.vector(combine.mcmc(as.mcmc.list(rjm_ch4_H,vars="ch4ETS")))))*test_complex$edge_tundra_H
pred_GT_ch4_H=t(replicate(n=length(test_complex$ch4),as.vector(combine.mcmc(as.mcmc.list(rjm_ch4_H,vars="ch4GTS")))))*test_complex$green_tundra_H

pred_ch4_H=pred_water_ch4_H+pred_fen_ch4_H+pred_deg_ch4_H+pred_degE_ch4_H+pred_LT_ch4_H+pred_ET_ch4_H+pred_GT_ch4_H
residuals_ch4_H=pred_ch4_H-test_complex$ch4
residuals_ch4_H=na.omit(residuals_ch4_H)
RMSE_ch4_H=sqrt(colSums((residuals_ch4_H)^2)/nrow(residuals_ch4_H))
hist(RMSE_ch4_H)

# Kormann & Meixner
pred_water_ch4_ZM=t(replicate(n=length(test_complex$ch4),as.vector(combine.mcmc(as.mcmc.list(rjm_ch4_ZM,vars="ch4FW")))))*test_complex$water_ZM
pred_fen_ch4_ZM=t(replicate(n=length(test_complex$ch4),as.vector(combine.mcmc(as.mcmc.list(rjm_ch4_ZM,vars="ch4FS1")))))*(test_complex$fen1_ZM+test_complex$fen2_ZM)
pred_deg_ch4_ZM=t(replicate(n=length(test_complex$ch4),as.vector(combine.mcmc(as.mcmc.list(rjm_ch4_ZM,vars="ch4DS")))))*test_complex$degraded_ZM
pred_degE_ch4_ZM=t(replicate(n=length(test_complex$ch4),as.vector(combine.mcmc(as.mcmc.list(rjm_ch4_ZM,vars="ch4EDS")))))*test_complex$edge_degraded_ZM
pred_LT_ch4_ZM=t(replicate(n=length(test_complex$ch4),as.vector(combine.mcmc(as.mcmc.list(rjm_ch4_ZM,vars="ch4LTS")))))*test_complex$lichen_tundra_ZM
pred_ET_ch4_ZM=t(replicate(n=length(test_complex$ch4),as.vector(combine.mcmc(as.mcmc.list(rjm_ch4_ZM,vars="ch4ETS")))))*test_complex$edge_tundra_ZM
pred_GT_ch4_ZM=t(replicate(n=length(test_complex$ch4),as.vector(combine.mcmc(as.mcmc.list(rjm_ch4_ZM,vars="ch4GTS")))))*test_complex$green_tundra_ZM

pred_ch4_ZM=pred_water_ch4_ZM+pred_fen_ch4_ZM+pred_deg_ch4_ZM+pred_degE_ch4_ZM+pred_LT_ch4_ZM+pred_ET_ch4_ZM+pred_GT_ch4_ZM
residuals_ch4_ZM=pred_ch4_ZM-test_complex$ch4
residuals_ch4_ZM=na.omit(residuals_ch4_ZM)
RMSE_ch4_ZM=sqrt(colSums((residuals_ch4_ZM)^2)/nrow(residuals_ch4_ZM))
hist(RMSE_ch4_ZM)

# export posterior chains for landcover ch4 predictions

# Kormann and Meixner
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_ch4_ZM,vars="ch4FW"))))),paste0('ch4FW_ZM',analysis_date,run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_ch4_ZM,vars="ch4DS"))))),paste0('ch4DS_ZM',analysis_date,run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_ch4_ZM,vars="ch4EDS"))))),paste0('ch4EDS_ZM',analysis_date,run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_ch4_ZM,vars="ch4FS1"))))),paste0('ch4FS1_ZM',analysis_date,run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_ch4_ZM,vars="ch4LTS"))))),paste0('ch4LTS_ZM',analysis_date,run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_ch4_ZM,vars="ch4ETS"))))),paste0('ch4ETS_ZM',analysis_date,run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_ch4_ZM,vars="ch4GTS"))))),paste0('ch4GTS_ZM',analysis_date,run_version,'.csv'))

# Hsieh
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_ch4_H,vars="ch4FW"))))),paste0('ch4FW_H',analysis_date,run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_ch4_H,vars="ch4DS"))))),paste0('ch4DS_H',analysis_date,run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_ch4_H,vars="ch4EDS"))))),paste0('ch4EDS_H',analysis_date,run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_ch4_H,vars="ch4FS1"))))),paste0('ch4FS1_H',analysis_date,run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_ch4_H,vars="ch4LTS"))))),paste0('ch4LTS_H',analysis_date,run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_ch4_H,vars="ch4ETS"))))),paste0('ch4ETS_H',analysis_date,run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_ch4_H,vars="ch4GTS"))))),paste0('ch4GTS_H',analysis_date,run_version,'.csv'))

# Kljun
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_ch4_Klj,vars="ch4FW"))))),paste0('ch4FW_Klj',analysis_date,run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_ch4_Klj,vars="ch4DS"))))),paste0('ch4DS_Klj',analysis_date,run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_ch4_Klj,vars="ch4EDS"))))),paste0('ch4EDS_Klj',analysis_date,run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_ch4_Klj,vars="ch4FS1"))))),paste0('ch4FS1_Klj',analysis_date,run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_ch4_Klj,vars="ch4LTS"))))),paste0('ch4LTS_Klj',analysis_date,run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_ch4_Klj,vars="ch4ETS"))))),paste0('ch4ETS_Klj',analysis_date,run_version,'.csv'))
write_csv(data.frame(t(as.matrix(combine.mcmc(as.mcmc.list(rjm_ch4_Klj,vars="ch4GTS"))))),paste0('ch4GTS_Klj',analysis_date,run_version,'.csv'))





