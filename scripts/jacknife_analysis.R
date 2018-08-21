#!/bin/env Rscript
# rm(list=ls())
# setwd('/Users/robins64/Documents/git_repos/seychelles-benthos')

library(rethinking)
library(dplyr)
library(stringr)
library(tidyr)
source('plot-cor-functions.R')
source('scaling_function.R')

## read in predicted recovery years
load('data/bayes_predicted_recovery_time.Rdata'); rates<-base
rates$baseline.minus<-NULL; rates$baseline.plus<-NULL


## add predictors
load(file='data/recovery_predictors_clean.Rdata')
pred$Region<-gsub('_', ' ', pred$Region)

rates$Herb_biomass<-pred$Herb_biomass[match(rates$location, pred$Region)]
rates$Depth<-pred$Depth[match(rates$location, pred$Region)]
rates$Init_complex<-pred$Init_complex[match(rates$location, pred$Region)]
rates$Init_totalcoral<-pred$Init_totalcoral[match(rates$location, pred$Region)]
rates$Nutrients_CN_ratio<-pred$Nutrients_CN_ratio[match(rates$location, pred$Region)]
rates$Wave_exposure_joules<-pred$Wave_exposure_joules[match(rates$location, pred$Region)]
rates$Sea_urchin_m2_2008<-pred$Sea_urchin_m2_2008[match(rates$location, pred$Region)]
rates$Coral_Juv_dens<-pred$Coral_Juv_dens[match(rates$location, pred$Region)]
rates$Manage<-pred$Manage[match(rates$location, pred$Region)]

rates<-data.frame(rates)

colnames(rates)[colnames(rates)=='recovery.year']<-'recoveryyear'

## change recovery year to real time
rates$recoveryyear<-rates$recoveryyear+6

## setup post list for saving posterior distributions
post.list<-numeric()

pdf(file='results/jacknife/recovery_diagnostic_jackknife.pdf', height=12, width=12)

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

	save(m, scaled, post, file=paste('results/jacknife/jacknife', i, '.Rdata'))

	## inspect param estimates
	plot(precis(m, depth=2), main=paste('jacknife_', i))

	## add posteriors to master list
	post.list[[paste0("jack", i)]] <- post

}

dev.off()


save(post.list, file='results/jacknife/jacknife_posts.Rdata')

