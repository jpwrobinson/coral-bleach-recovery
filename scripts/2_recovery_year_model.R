#!/bin/env Rscript


library(rethinking)
library(here)

setwd(here())

source('scripts/scaling_function.R')


## ------- ------- ------- ------- ------- ------- ------- ##
              ### Fit model to estimate recovery year ###
## ------- ------- ------- ------- ------- ------- ------- ##

## load predictors (can be provided by request to the authors)
load(file='data/recovery_predictors_clean.Rdata') ## data frame is 'rates'


## change recovery year to real time
rates$recoveryyear<-rates$recoveryyear+6

## scale exp. covariates to mean = 0 
scaled<-scaler(rates, ID = c('location', 'recoveryyear'))

### rec.year ~ predictors; linear model

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
	    	  rateh*N_percent +
	    	  ratej*Manage,

	    ## fixed priors for exp. covariates
	    c(rateb, ratec, rated, ratee, ratef, rateg, rateh, ratej) ~ dnorm(0, 2),
	    
	    ## priors from mean year value
	    c(ratea) ~ dnorm(17, 5),

	    ## error priors
	    c(sigma) ~ dcauchy(0 , 2 )

), data=scaled, warmup=1500, iter=7000, chains = 1)


rec.year.m<-m

save(rec.year.m, file='data/recovery_year_model.Rdata')

