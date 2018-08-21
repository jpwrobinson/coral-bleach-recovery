#!/bin/env Rscript

library(rethinking)
library(tidyverse)
library(here)

setwd(here())
## function to scale and center continuous predictors (and dummy variables for categorical)
source('scaling_function.R')

## load UVC benthic data (can be provided by request to the authors)
load('data/SC_site.Rdata') ## dataframe is 'focal'

## ------- ------- ------- ------- ------- ------- ------- ##
              ### Fit models to total coral ###
## ------- ------- ------- ------- ------- ------- ------- ##

## scale exp params
scaled<-scaler(focal, ID = c('state', 'location', 'cover', 'year', 'availsubstrate'))
scaled$cover<-round(scaled$cover, 0)

## change 0s to 1% cover for distribution fitting (this is 2 replicates in 2005)
scaled$cover[scaled$cover==0]<-1
## log-transform cover and available substrate
scaled$logcover<-log(scaled$cover)
scaled$logavailsubstrate<-log(scaled$availsubstrate)

## change year to recovery years
scaled$year<-scaled$year-2004


### Compare fits from different growth models
pdf(file='figures/logistic_models_predictions.pdf', height=7, width = 14)

## 1) Simple 3-param logistic (Osborne et al. 2011), from SSlogis
## intrinsic growth rate and asymptote by reef, estimated from observations with no other constraints

initVals.logistic <- getInitial(logcover ~ SSlogis(year, Asym, xmid, scal), data = scaled);initVals.logistic

m1 <- map2stan(
	alist(

		## response distribution
	    cover ~ dgampois( mu , scale ) ,
	    
	    ## model structure
	    log(mu) <- Asym/(1 + exp((xmid-year)/rate)), ## 3 param logistic, SSlogis
	    Asym <- Asymfix + Asymr[location],
		rate <- ratefix + rater[location],


	    ## random priors
	    c(Asymr, rater)[location] ~  dmvnorm2(0, rateError, Rho), 

	    ## fixed priors from model without exp. params
	    c(Asymfix) ~ dnorm(3.6, 1),
	    c(xmid) ~ dnorm(-0.9, 1),
	    c(ratefix) ~ dnorm(6, 1),

	    ## error priors
	    c(scale, rateError) ~ dcauchy(0 , 2 ),
	    Rho ~ dlkjcorr(4)

), data=scaled, warmup=1500, iter=3000, chains=3,constraints=list(scale="lower=0"))

precis(m1, 2)

## 2) simple 3-param logistic growth with constrained K (Osborne et al. 2017)
## K is constrained to available hard substrate area (rock + coral + rubble)

initVals.logistic <- getInitial(logcover ~ SSlogis(year, logavailsubstrate, xmid, scal), data = scaled);initVals.logistic


m2 <- map2stan(
	alist(

		## response distribution
	    cover ~ dgampois( mu , scale ) ,
	    
	    ## model structure
		log(mu) <- logavailsubstrate/(1 + exp((xmid-year)/rate)), 
		rate <- ratefix + rater[location],

	    ## random priors
	   	rater[location] ~ dnorm(0, sigmar), 

	    ## fixed priors from model without exp. params
	    c(xmid) ~ dnorm(-0.9, 1),
	    c(ratefix) ~ dnorm(6, 1),

	    ## error priors
	    c(scale, sigmar) ~ dcauchy(0 , 2 )

), data=scaled,  warmup=1500, iter=3000, chains=3, constraints=list(scale="lower=0"))

precis(m2, depth=2)

## 3) 4-param logistic growth 
## K is unconstrained, and for 2 asymptotes

m3 <- map2stan(
	alist(

		## response distribution
	    cover ~ dgampois( mu , scale ) ,
	    
	    ## model structure
	    log(mu) <-  Asym+(Asym2-Asym)/(1+exp((xmid-year)/rate)), ## 4 param logistic
	    Asym <- Asymfix + Asymr[location],
	    Asym2 <- Asymfix2 + Asymr2[location],
		rate <- ratefix + rater[location],



	    ## random priors
	    c(Asymr, Asymr2, rater)[location] ~  dmvnorm2(0, rateError, Rho), 

	    ## fixed priors from model without exp. params
	    c(Asymfix) ~ dnorm(3.6, 1),
	    c(Asymfix2) ~ dnorm(3.6, 1),
	    c(xmid) ~ dnorm(-0.9, 1),
	    c(ratefix) ~ dnorm(6, 1),

	    ## error priors
	    c(scale, rateError) ~ dcauchy(0 , 2 ),
	    Rho ~ dlkjcorr(4)

), data=scaled, warmup=1500, iter=3000, chains=3,constraints=list(scale="lower=0"))

precis(m3, depth = 2)

## 4) 4-param logistic growth 
## Asym2 is constrained to available substrate, and for 2 asymptotes

m4 <- map2stan(
	alist(

		## response distribution
	    cover ~ dgampois( mu , scale ) ,
	    
	    ## model structure
	    log(mu) <-  Asym+(Asym2-Asym)/(1+exp((xmid-year)/rate)), ## 4 param logistic
	    Asym <- Asymfix + Asymr[location],
	    Asym2 <- logavailsubstrate,
		rate <- ratefix + rater[location],



	    ## random priors
	    c(Asymr, rater)[location] ~  dmvnorm2(0, rateError, Rho), 

	    ## fixed priors from model without exp. params
	    c(Asymfix) ~ dnorm(3.6, 1),
	    # c(Asymfix2) ~ dnorm(3.6, 1),
	    c(xmid) ~ dnorm(-0.9, 1),
	    c(ratefix) ~ dnorm(6, 1),

	    ## error priors
	    c(scale, rateError) ~ dcauchy(0 , 2 ),
	    Rho ~ dlkjcorr(4)

), data=scaled, warmup=1500, iter=3000, chains=3,constraints=list(scale="lower=0"))

precis(m4, depth = 2)

### Simulate predicted cover for each model, compare fits
year.vec<-c(1:10)
pred<-expand.grid(location = unique(scaled$location), year = year.vec)
pred$logavailsubstrate<-scaled$logavailsubstrate[match(pred$location, scaled$location)]
pred$availsubstrate<-scaled$availsubstrate[match(pred$location, scaled$location)]

## simulate post predictions

## Model 1 
means<-link(m1, data=pred)
pred.cover<-apply(means$mu, 2, mean)
pred.cover.PI <- apply( means$mu , 2 , HPDI , prob=0.97 )
##. convert into plotable form
site.preds1<-pred
site.preds1$pcover<-pred.cover
site.preds1$lower<-pred.cover.PI[1,]
site.preds1$upper<-pred.cover.PI[2,]

ggplot(site.preds1, aes(year, pcover)) + 
			geom_point(data = scaled, aes(year, cover)) +
			geom_line() + geom_ribbon(aes(ymin=lower, ymax=upper), fill=alpha('grey', 0.5)) +
			facet_wrap(~location) + labs(title='Model 1 - simple logistic')

## Model 2
means<-link(m2, data=pred)
pred.cover<-apply(means$mu, 2, mean)
pred.cover.PI <- apply( means$mu , 2 , HPDI , prob=0.97 )
##. convert into plotable form
site.preds2<-pred
site.preds2$pcover<-pred.cover
site.preds2$lower<-pred.cover.PI[1,]
site.preds2$upper<-pred.cover.PI[2,]

ggplot(site.preds2, aes(year, pcover)) + 
			geom_point(data = scaled, aes(year, cover)) +
			geom_line() + geom_ribbon(aes(ymin=lower, ymax=upper), fill=alpha('grey', 0.5)) +
			facet_wrap(~location) + labs(title='Model 2 - simple logistic contrained K')

## Model 3
means<-link(m3, data=pred)
pred.cover<-apply(means$mu, 2, mean)
pred.cover.PI <- apply( means$mu , 2 , HPDI , prob=0.97 )
##. convert into plotable form
site.preds3<-pred
site.preds3$pcover<-pred.cover
site.preds3$lower<-pred.cover.PI[1,]
site.preds3$upper<-pred.cover.PI[2,]

ggplot(site.preds3, aes(year, pcover)) + 
			geom_point(data = scaled, aes(year, cover)) +
			geom_line() + geom_ribbon(aes(ymin=lower, ymax=upper), fill=alpha('grey', 0.5)) +
			facet_wrap(~location) + labs(title='Model 3 - 4 param logistic')

## Model 4 
means<-link(m4, data=pred)
pred.cover<-apply(means$mu, 2, mean)
pred.cover.PI <- apply( means$mu , 2 , HPDI , prob=0.97 )
##. convert into plotable form
site.preds4<-pred
site.preds4$pcover<-pred.cover
site.preds4$lower<-pred.cover.PI[1,]
site.preds4$upper<-pred.cover.PI[2,]

ggplot(site.preds4, aes(year, pcover)) + 
			geom_point(data = scaled, aes(year, cover)) +
			geom_line() + geom_ribbon(aes(ymin=lower, ymax=upper), fill=alpha('grey', 0.6)) +
			facet_wrap(~location) + labs(title='Model 4 - 4 param logistic, constrained K')



waic<-compare(m1, m2, m3, m4)
waic<-waic@output

write.csv(data.frame(waic), 'results/waic_logistic_compare.csv')


## save top model
rec.trajectory.m<-m4
save(rec.trajectory.m, file='results/recovery_year_model.Rdata')


dev.off()
