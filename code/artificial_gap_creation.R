
# Code for generating artificial gaps in carbon fluxes

# author: Sarah M. Ludwig

# to be used in conjunction with Resolving_heterogeneous_tundra_C_fluxes_Ludwig_2023.R

#### artificial gap parameters: ####

# these I changed to create different random gaps as desired
seed_num_long=1
seed_num_med=2
seed_num_small=3

# these are the same for every analysis
long_gap=10
med_gap=4

#### for the complex map footprint influence dataset: ####

data_2020=subset(df_complex,df_complex$datetime>as.POSIXct('2020-04-01 00:00:00',format="%Y-%m-%d %H:%M",tz='Etc/GMT+9'))
data_2019=subset(df_complex,df_complex$datetime<as.POSIXct('2020-04-01 00:00:00',format="%Y-%m-%d %H:%M",tz='Etc/GMT+9'))
# need to create  new 'gaps'
all_obs_int=which(!is.na(data_2020$co2))
set.seed(seed_num_long) 
# find positions of 1% of obs for random long gaps
withhold_int_long=sample.int(length(all_obs_int),floor(length(data_2020$co2)*0.01))
start_index_long=all_obs_int[withhold_int_long]
missing_vector=vector(mode='numeric',length=length(start_index_long))
withheld_index_long=start_index_long
for (i in 1:long_gap){
  missing_vector=start_index_long+i
  withheld_index_long=c(withheld_index_long,missing_vector)
}
train_index_long=all_obs_int[all_obs_int%nin%withheld_index_long]

# find positions of 1% of remaining obs for random medium gaps
set.seed(seed_num_med)
withhold_int_med=sample.int(length(train_index_long),floor(length(data_2020$co2)*0.01))
start_index_med=train_index_long[withhold_int_med]
missing_vector=vector(mode='numeric',length=length(start_index_med))
withheld_index_med=start_index_med
for (i in 1:med_gap){
  missing_vector=start_index_med+i
  withheld_index_med=c(withheld_index_med,missing_vector)
}
train_index_long_med=train_index_long[train_index_long%nin%withheld_index_med]

# find positions of 5% of remaining obs for random small gaps
set.seed(seed_num_small)
withhold_int_small=sample.int(length(train_index_long_med),floor(length(data_2020$co2)*0.05))
withheld_index_small=train_index_long_med[withhold_int_small]
train_index=train_index_long_med[train_index_long_med%nin%withheld_index_small]

train_complex=rbind(data_2019,df_complex[train_index+length(data_2019$co2),])
train_complex=subset(train_complex,!is.na(train_complex$co2))
train_complex=subset(train_complex,!is.na(train_complex$air_temperature))
train_complex=subset(train_complex,!is.na(train_complex$PPFD_1_1_1))
train_complex=subset(train_complex,!is.na(train_complex$lichen_tundra_Klj))

withheld_index_all=c(withheld_index_small,withheld_index_med,withheld_index_long)
withheld_index_all=unique(withheld_index_all)
test_complex=df_complex[withheld_index_all+length(data_2019$co2),]
withheld_index=withheld_index_all[which(!is.na(test_complex$co2))]
test_complex=subset(test_complex,!is.na(test_complex$co2))

obs_mask_df=subset(df_complex,df_complex$datetime>as.POSIXct('2020-04-01 00:00:00',format="%Y-%m-%d %H:%M",tz='Etc/GMT+9'))
for (i in 1:length(withheld_index)){
  indice=withheld_index[i]
  obs_mask_df$co2[indice]=NA
}
obs_mask=is.na(obs_mask_df$co2)

#### for the simple map footprint influence dataset: ####

data_2020=subset(df_simple,df_simple$datetime>as.POSIXct('2020-04-01 00:00:00',format="%Y-%m-%d %H:%M",tz='Etc/GMT+9'))
data_2019=subset(df_simple,df_simple$datetime<as.POSIXct('2020-04-01 00:00:00',format="%Y-%m-%d %H:%M",tz='Etc/GMT+9'))
# need to create  new 'gaps'
all_obs_int=which(!is.na(data_2020$co2))
set.seed(seed_num_long) #512 good seed
# find positions of 1% of obs for random long gaps
withhold_int_long=sample.int(length(all_obs_int),floor(length(data_2020$co2)*0.01))
start_index_long=all_obs_int[withhold_int_long]
missing_vector=vector(mode='numeric',length=length(start_index_long))
withheld_index_long=start_index_long
for (i in 1:long_gap){
  missing_vector=start_index_long+i
  withheld_index_long=c(withheld_index_long,missing_vector)
}
train_index_long=all_obs_int[all_obs_int%nin%withheld_index_long]
# find positions of 1% of remaining obs for random medium gaps
set.seed(seed_num_med)
withhold_int_med=sample.int(length(train_index_long),floor(length(data_2020$co2)*0.01))
start_index_med=train_index_long[withhold_int_med]
missing_vector=vector(mode='numeric',length=length(start_index_med))
withheld_index_med=start_index_med
for (i in 1:med_gap){
  missing_vector=start_index_med+i
  withheld_index_med=c(withheld_index_med,missing_vector)
}
train_index_long_med=train_index_long[train_index_long%nin%withheld_index_med]

# find positions of 5% of remaining obs for random small gaps
set.seed(seed_num_small)
withhold_int_small=sample.int(length(train_index_long_med),floor(length(data_2020$co2)*0.05))
withheld_index_small=train_index_long_med[withhold_int_small]
train_index=train_index_long_med[train_index_long_med%nin%withheld_index_small]


train_simple=rbind(data_2019,df_simple[train_index+length(data_2019$co2),])
train_simple=subset(train_simple,!is.na(train_simple$co2))
no_temp=which(is.na(train_simple$air_temperature))
removed_temp=train_simple[no_temp,]
train_simple=subset(train_simple,!is.na(train_simple$air_temperature))
train_simple=subset(train_simple,!is.na(train_simple$PPFD_1_1_1))
train_simple=subset(train_simple,!is.na(train_simple$tundra_H))

withheld_index_all=c(withheld_index_small,withheld_index_med,withheld_index_long)
withheld_index_all_unq=unique(withheld_index_all)
test_simple=df_simple[withheld_index_all_unq+length(data_2019$co2),]
withheld_index=withheld_index_all_unq[which(!is.na(test_simple$co2))]
test_simple=subset(test_simple,!is.na(test_simple$co2))

obs_mask_df=subset(df_simple,df_simple$datetime>as.POSIXct('2020-04-01 00:00:00',format="%Y-%m-%d %H:%M",tz='Etc/GMT+9'))
temp_mask=which(obs_mask_df$datetime%in%removed_temp$datetime)

for (i in 1:length(c(withheld_index,temp_mask))){
  indice=c(withheld_index,temp_mask)[i]
  obs_mask_df$co2[indice]=NA
}
obs_mask=is.na(obs_mask_df$co2)

