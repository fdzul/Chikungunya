######################################################################
# ruizMoreno.R                                                       #
# Author:  Dario Ghersi                                              #
# Version: 20140918                                                  #
# Goal:    implementation of the ODE system in:                      #
#          Ruiz-Moreno et al, PLoS Neglected Tropical Diseases 2012  #
######################################################################

require(deSolve)

returnTempDepParams <- function(temp) {
    ## return the temperature dependent parameters

    ## mosquito egg mortality (1/survival, capped at 1)
    if (temp <= 5) {
        mG <- 1
    }
    else if (temp > 5 && temp <= 15) {
        mG <- min(1 / (9.4 - 1.46 * temp + 0.092 * temp^2), 1)
    }
    else if (temp > 15 && temp <= 20) {
        mG <- min(1 / (-167.9 + 11.74 * temp), 1)
    }
    else if (temp > 20 && temp <= 39.54) {
        mG <- min(1 / (-5.08 + 7.23 * temp - 0.18 * temp^2), 1)
    }
    else {
        mG <- 1
    }
    
    ## proportion of eggs entering diapause state
    if (temp < 0) {
        zG <- 0
    }
    else if (temp >=0 && temp <= 10) {
        zG <- -0.1 * temp + 1
    }
    else { # temp > 10
        zG <- 1
    }

    ## proportion of eggs leaving diapause state
    if (temp < 0) {
        vD <- 0
    }
    else if (temp >= 0 && temp <= 10) {
        vD = 0.1 * temp
    }
    else { # temp > 10
        vD <- 1
    }

    ## 1 / egg maturation period
    xG <- 1 / (13.18 - 0.845 * temp + 0.016 * temp^2)

    ## mortality of larvae (1 / survival of larvae, capped at 1)
    if (temp < 10) {
        mL <- 1
    }
    else if (temp >= 10 && temp <= 38.5) {
        mL <- min(1 / (-136.08 + 17.8 * temp -0.37 * temp^2), 1)
    }
    else { # temp > 38.5
        mL <- 1
    }

    ## 1 / larvae maturation
    xL <- 1 / (117.58 - 7.6 * temp + 0.1321528 * temp^2)

    ## adult mortality (1 / adult survival, capped at 1)
    if (temp < -5) {
        mAs <- 1
    }
    else if (temp >= -5 && temp < 15) {
        mAs <- 1 / (5.6 + 1.12 * temp)
    }
    else if (temp > 15 && temp <= 38.5) {
        mAs <- 1 / (40.27 -1.38 * temp + 0.01 * temp^2)
    }
    else { # temp > 38.5
        mAs <- 1
    }
    mAe <- mAi <- mAs # the mortality of the adults is the same
                      # regardless of their infective state

    ## 1 / extrinsic incubation period
    if (temp < 10) {
        sAe <- 1 / 4
    }
    else if (temp >= 10 && temp <= 32) {
        sAe <- 1 / (113 / 22 - (2.4 / 22) * temp)
    }
    else { # temp > 32
        sAe <- 1 / (1.5)
    }

    ## put all the parameters in a list
    params <- list(mG=mG, zG=zG, vD=vD, xG=xG, mL=mL, xL=xL, mAs=mAs,
                   mAe=mAe, mAi=mAi, sAe=sAe)
    
    return(params)
}

######################################################################

returnParams <- function() {
    ## set up the temperature-independent parameters
    
    ## adult mosquito population
    Am<-10000000
    
    ## number of mosquito eggs per mosquito per day
    f <- 4.4
    
    ## symptomatic/asymptomaic
    c <- 0.75
   
    ## biting rate
    b <- 1/3

    ## infective period
    inf <- 1/14

    ## probability of human to mosquito transition
    Tmh <- 0.977

    ## probability of mosquito to human transmission
    Thm <- 0.726
 
    ## density dependent factor for reproduction (median)
    Pm <- 0.0000002 
    
    ## natural mortality rate of human adults 
    mu <- 7.9/1000 
    
    ## disease induced mortality rate of human adults
    dmu <- 1/1000

    ##density function
    rH <- 0.000002



    ## put all the parameters in a list
    params <- list(f=f, c=c, Pm=Pm, mu=mu,b=b,Thm=Thm,inf=inf,dmu=dmu,Am=Am,Tmh=Tmh,rH=rH)
    
    return(params)
}

returnInits <- function() {
    ## number of eggs
	G <- 1000
	
	## number of eggs undergoing diapause
	D<- 1000

	## number of larvae
	L <- 1000

	## number of susceptible adults
	As<-1000
	
	## number of exposed adults
	Ae <- 4000

	#number of infected adults
	Ai <- 10
	
    ## put all the parameters in a list
    params <- list(G=G,D=D,L=L,As=As,Ae=Ae,Ai=Ai,rH=rH)
    
    return(params)
}
