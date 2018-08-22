rm(list=ls())


library(rethinking)
library(here)
setwd(here())

source('scripts/scaling_function.R')

## read in predicted recovery years, rename 'base' to rates
load('data/recovery_trajectory.Rdata'); rates<-base

## load reef-scale predictors (not uploaded - can be provided by request to the authors)
load(file='data/recovery_predictors_clean.Rdata')

## match in reef predictors to recovery year dataframe
rates$Herb_biomass<-pred$Herb_biomass[match(rates$location, pred$Region)]
rates$Depth<-pred$Depth[match(rates$location, pred$Region)]
rates$Init_complex<-pred$Init_complex[match(rates$location, pred$Region)]
rates$Init_totalcoral<-pred$Init_totalcoral[match(rates$location, pred$Region)]
rates$Nutrients_CN_ratio<-pred$Nutrients_CN_ratio[match(rates$location, pred$Region)]
rates$Wave_exposure_joules<-pred$Wave_exposure_joules[match(rates$location, pred$Region)]
rates$Sea_urchin_m2_2008<-pred$Sea_urchin_m2_2008[match(rates$location, pred$Region)]
rates$Coral_Juv_dens<-pred$Coral_Juv_dens[match(rates$location, pred$Region)]
rates$Manage<-pred$Manage[match(rates$location, pred$Region)]

## strip dots for rethinking model
colnames(rates)[colnames(rates)=='recovery.year']<-'recoveryyear'

## change recovery year to real time
rates$recoveryyear<-rates$recoveryyear+6


pdf(file='figures/recovery_diagnostic_jackknife.pdf', height=12, width=12)

## Loop through dataset with one data point removed each iteration
# fit linear model predicting recovery year
# save model and posteriors

for(i in 1:12){

	scaled<-scaler(rates[-i,], ID = c('location', 'recoveryyear'))

	### rates ~ predictors
		m <- map2stan(
			alist(

				## response distribution
			    recoveryyear ~ dnorm( mu , sigma) ,
			    
			    ## model structure
			    mu <- ratea + rateb*Herb_biomass +
					  rated*Depth + 
					  ratec*Coral_Juv_dens +
			    	  ratee*Init_complex +
			    	  ratef*Init_totalcoral + 
			    	  rateg*Wave_exposure_joules + 
			    	  rateh*Nutrients_CN_ratio +
			    	  ratej*Manage,

			    ## fixed priors for exp. covariates
			    c(rateb, ratec, rated, ratee, ratef, rateg, rateh, ratej) ~ dnorm(0, 2),
			    
			    ## priors from mean year value
			    c(ratea) ~ dnorm(17, 5),

			    ## error priors
			    c(sigma) ~ dcauchy(0 , 2 )

		), data=scaled, warmup=1500, iter=3000, chains=1)

	## extract posteriors
	post<-extract.samples(m, n = 1000)

	save(m, scaled, post, file=paste('data/jacknife', i, '.Rdata'))

	## inspect param estimates
	plot(precis(m, depth=2), main=paste('jacknife_', i))

}

dev.off()




