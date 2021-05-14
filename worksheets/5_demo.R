#' ---
#' title: "Class 6 live demo: Transition estimation"
#' author: "The team"
#' date: "last updated: `r Sys.Date()`"
#' output: html_document
#' ---
#' 
#' 
## ----setup, include=FALSE, echo=FALSE---------------------------------------------------------------------
options(htmltools.dir.version = FALSE)
knitr::opts_chunk$set(comment = "", echo = TRUE, warning = FALSE, message = FALSE)
library(tidyverse)
theme_set(theme_light(base_size = 14))
update_geom_defaults("point", list(size = 2)) 
library(here)
library(nimble)

#' 
#' ## Introduction
#' 
#' In this demo, we illustrate how to fit multistate models to capture-recapture data. We first consider states as sites and estimate movements. Next we treat states as breeding states for studying life-history trade-offs. 
#' 
## ---------------------------------------------------------------------------------------------------------
library(tidyverse)
library(nimble)
library(MCMCvis)

#' 
#' 
#' ## Multisite model with 2 sites
#' 
#' We are going to analyze real capture-recapture data on > 20,000 Canada geese (Branta canadensis) marked and resighted in the east coast of the US in the mid-Atlantic (New York, Pennsylvania, New Jersey), Chesapeake (Delaware, Maryland, Virginia), and Carolinas (North and South Carolina). The data were kindly provided by Jay Hestbeck. We analyse a subset of 500 geese from the original dataset, we'll get back to the whole dataset in the last live demo. 
#' 
#' Read in the data. 
## ---------------------------------------------------------------------------------------------------------
geese <- read_csv("geese.csv", col_names = TRUE)
y <- as.matrix(geese)
head(y)

#' 
#' For simplicity, we start by considering only 2 sites, dropping Carolinas. 
## ---------------------------------------------------------------------------------------------------------
y2 <- y
y2[y2==3] <- 0 # act as if there was no detection in site 3 Carolinas
mask <- apply(y2, 1, sum)
y2 <- y2[mask!=0,] # remove rows w/ 0s only
head(y2)

#' 
#' Get occasion of first capture.
## ---------------------------------------------------------------------------------------------------------
get.first <- function(x) min(which(x != 0))
first <- apply(y2, 1, get.first)
first

#' 
#' Write the model. 
## ---------------------------------------------------------------------------------------------------------
multisite <- nimbleCode({
  
  # -------------------------------------------------
  # Parameters:
  # phiA: survival probability site A
  # phiB: survival probability site B
  # psiAB: movement probability from site A to site B
  # psiBA: movement probability from site B to site A
  # pA: recapture probability site A
  # pB: recapture probability site B
  # piA: prop of being in site A at initial capture
  # -------------------------------------------------
  # States (z):
  # 1 alive at A
  # 2 alive at B
  # 3 dead
  # Observations (y):  
  # 1 not seen
  # 2 seen at site A 
  # 3 seen at site B
  # -------------------------------------------------
  
  # priors
  phiA ~ dunif(0, 1)
  phiB ~ dunif(0, 1)
  psiAB ~ dunif(0, 1)
  psiBA ~ dunif(0, 1)
  pA ~ dunif(0, 1)
  pB ~ dunif(0, 1)
  piA ~ dunif(0, 1)
  
  # probabilities of state z(t+1) given z(t)
  gamma[1,1] <- phiA * (1 - psiAB)
  gamma[1,2] <- phiA * psiAB
  gamma[1,3] <- 1 - phiA
  gamma[2,1] <- phiB * psiBA
  gamma[2,2] <- phiB * (1 - psiBA)
  gamma[2,3] <- 1 - phiB
  gamma[3,1] <- 0
  gamma[3,2] <- 0
  gamma[3,3] <- 1
  
  delta[1] <- piA          # Pr(alive in A t = 1)
  delta[2] <- 1 - piA      # Pr(alive in B t = 1)
  delta[3] <- 0            # Pr(dead t = 1) = 0
  
  # probabilities of y(t) given z(t)
  omega[1,1] <- 1 - pA     # Pr(alive A t -> non-detected t)
  omega[1,2] <- pA         # Pr(alive A t -> detected A t)
  omega[1,3] <- 0          # Pr(alive A t -> detected B t)
  omega[2,1] <- 1 - pB     # Pr(alive B t -> non-detected t)
  omega[2,2] <- 0          # Pr(alive B t -> detected A t)
  omega[2,3] <- pB         # Pr(alive B t -> detected B t)
  omega[3,1] <- 1          # Pr(dead t -> non-detected t)
  omega[3,2] <- 0          # Pr(dead t -> detected A t)
  omega[3,3] <- 0          # Pr(dead t -> detected B t)
  
  # likelihood 
  for (i in 1:N){
    # latent state at first capture
    z[i,first[i]] ~ dcat(delta[1:3])
    for (t in (first[i]+1):K){
      # draw z(t) given z(t-1)
      z[i,t] ~ dcat(gamma[z[i,t-1],1:3])
      # draw y(t) given z(t)
      y[i,t] ~ dcat(omega[z[i,t],1:3])
    }
  }
})

#' 
#' Data in a list. Remember to add 1. 
## ---------------------------------------------------------------------------------------------------------
my.data <- list(y = y2 + 1)

#' 
#' Constants in a list. 
## ---------------------------------------------------------------------------------------------------------
my.constants <- list(first = first, 
                     K = ncol(y2), 
                     N = nrow(y2))

#' 
#' Initial values. 
## ---------------------------------------------------------------------------------------------------------
zinits <- y2 # say states are observations, detections in A or B are taken as alive in same sites 
zinits[zinits==0] <- sample(c(1,2), sum(zinits==0), replace = TRUE) # non-detections become alive in site A or B
initial.values <- function(){list(phiA = runif(1, 0, 1), 
                                  phiB = runif(1, 0, 1), 
                                  psiAB = runif(1, 0, 1), 
                                  psiBA = runif(1, 0, 1), 
                                  pA = runif(1, 0, 1), 
                                  pB = runif(1, 0, 1), 
                                  piA = runif(1, 0, 1),
                                  z = zinits)}

#' 
#' Parameters to monitor.
## ---------------------------------------------------------------------------------------------------------
parameters.to.save <- c("phiA", "phiB","psiAB", "psiBA", "pA", "pB", "piA")
parameters.to.save

#' 
#' MCMC settings.
## ---------------------------------------------------------------------------------------------------------
n.iter <- 10000
n.burnin <- 5000
n.chains <- 2

#' 
#' Run nimble.
## ---------------------------------------------------------------------------------------------------------
mcmc.multisite <- nimbleMCMC(code = multisite, 
                             constants = my.constants,
                             data = my.data,              
                             inits = initial.values,
                             monitors = parameters.to.save,
                             niter = n.iter,
                             nburnin = n.burnin, 
                             nchains = n.chains)

#' 
#' Let's inspect the results. First the numerical summaries. 
## ---------------------------------------------------------------------------------------------------------
MCMCsummary(mcmc.multisite, round = 2)

#' 
#' Second, a caterpillar plot of the parameter estimates. 
## ---- echo = TRUE, fig.width = 7.5, fig.asp = 0.618, dev = "svg"------------------------------------------
MCMCplot(mcmc.multisite)

#' 
#' Actually, the initial state of a goose is known exactly. It is alive at site of initial capture. Therefore, we don't need to try and estimate the probability of initial states. Let's rewrite the model. 
## ---------------------------------------------------------------------------------------------------------
multisite <- nimbleCode({
  
  # -------------------------------------------------
  # Parameters:
  # phiA: survival probability site A
  # phiB: survival probability site B
  # psiAB: movement probability from site A to site B
  # psiBA: movement probability from site B to site A
  # pA: recapture probability site A
  # pB: recapture probability site B
  # -------------------------------------------------
  # States (z):
  # 1 alive at A
  # 2 alive at B
  # 3 dead
  # Observations (y):  
  # 1 not seen
  # 2 seen at A 
  # 3 seen at B
  # -------------------------------------------------
  
  # priors
  phiA ~ dunif(0, 1)
  phiB ~ dunif(0, 1)
  psiAB ~ dunif(0, 1)
  psiBA ~ dunif(0, 1)
  pA ~ dunif(0, 1)
  pB ~ dunif(0, 1)

  # probabilities of state z(t+1) given z(t)
  gamma[1,1] <- phiA * (1 - psiAB)
  gamma[1,2] <- phiA * psiAB
  gamma[1,3] <- 1 - phiA
  gamma[2,1] <- phiB * psiBA
  gamma[2,2] <- phiB * (1 - psiBA)
  gamma[2,3] <- 1 - phiB
  gamma[3,1] <- 0
  gamma[3,2] <- 0
  gamma[3,3] <- 1
  
  # probabilities of y(t) given z(t)
  omega[1,1] <- 1 - pA     # Pr(alive A t -> non-detected t)
  omega[1,2] <- pA         # Pr(alive A t -> detected A t)
  omega[1,3] <- 0          # Pr(alive A t -> detected B t)
  omega[2,1] <- 1 - pB     # Pr(alive B t -> non-detected t)
  omega[2,2] <- 0          # Pr(alive B t -> detected A t)
  omega[2,3] <- pB         # Pr(alive B t -> detected B t)
  omega[3,1] <- 1          # Pr(dead t -> non-detected t)
  omega[3,2] <- 0          # Pr(dead t -> detected A t)
  omega[3,3] <- 0          # Pr(dead t -> detected B t)
  
  # likelihood 
  for (i in 1:N){
    # latent state at first capture
    z[i,first[i]] <- y[i,first[i]] - 1 # if seen at site A (y = 2), state is alive in A (y - 1 = z = 1) with prob = 1
    for (t in (first[i]+1):K){         # if seen at site B (y = 3), state is alive in A (y - 1 = z = 2) with prob = 1
      # draw (t) given z(t-1)
      z[i,t] ~ dcat(gamma[z[i,t-1],1:3])
      # draw y(t) given z(t)
      y[i,t] ~ dcat(omega[z[i,t],1:3])
    }
  }
})

#' 
#' Initial values without $p_A$. 
## ---------------------------------------------------------------------------------------------------------
zinits <- y2
zinits[zinits==0] <- sample(c(1,2), sum(zinits==0), replace = TRUE)
initial.values <- function(){list(phiA = runif(1, 0, 1), 
                                  phiB = runif(1, 0, 1), 
                                  psiAB = runif(1, 0, 1), 
                                  psiBA = runif(1, 0, 1), 
                                  pA = runif(1, 0, 1), 
                                  pB = runif(1, 0, 1), 
                                  z = zinits)}

#' 
#' Parameters to monitor. 
## ---------------------------------------------------------------------------------------------------------
parameters.to.save <- c("phiA", "phiB","psiAB", "psiBA", "pA", "pB")
parameters.to.save

#' 
#' Run nimble.
## ---------------------------------------------------------------------------------------------------------
mcmc.multisite <- nimbleMCMC(code = multisite, 
                             constants = my.constants,
                             data = my.data,              
                             inits = initial.values,
                             monitors = parameters.to.save,
                             niter = n.iter,
                             nburnin = n.burnin, 
                             nchains = n.chains)

#' 
#' Have a look to the results. Note that convergence is better. 
## ---------------------------------------------------------------------------------------------------------
MCMCsummary(mcmc.multisite, round = 2)

#' 
#' ## Multisite model with 3 sites
#' 
#' Now we proceed with the analysis of the 3 sites. Two methods exist to force the movement probabilities from a site to sum to 1 and to be between 0 and 1 at the same time. 
#' 
#' ### Multinomial logit link
#' 
#' Get the date of first capture. 
## ---------------------------------------------------------------------------------------------------------
get.first <- function(x) min(which(x != 0))
first <- apply(y, 1, get.first)

#' 
#' Write the model. 
## ---------------------------------------------------------------------------------------------------------
multisite <- nimbleCode({
  
  # -------------------------------------------------
  # Parameters:
  # phiA: survival probability site A
  # phiB: survival probability site B
  # phiC: survival probability site B
  # psiAA: movement probability from site A to site A (reference)
  # psiAB = psiA[1]: movement probability from site A to site B
  # psiAC = psiA[2]: movement probability from site A to site C 
  # psiBA = psiB[1]: movement probability from site B to site A
  # psiBB: movement probability from site B to site B (reference)
  # psiBC = psiB[2]: movement probability from site B to site C
  # psiCA = psiC[1]: movement probability from site C to site A
  # psiCB = psiC[2]: movement probability from site C to site B
  # psiCC: movement probability from site C to site C (reference)
  # pA: recapture probability site A
  # pB: recapture probability site B
  # pC: recapture probability site C
  # -------------------------------------------------
  # States (z):
  # 1 alive at A
  # 2 alive at B
  # 2 alive at C
  # 3 dead
  # Observations (y):  
  # 1 not seen
  # 2 seen at A 
  # 3 seen at B
  # 3 seen at C
  # -------------------------------------------------
  
  # Priors
  phiA ~ dunif(0, 1)
  phiB ~ dunif(0, 1)
  phiC ~ dunif(0, 1)
  pA ~ dunif(0, 1)
  pB ~ dunif(0, 1)
  pC ~ dunif(0, 1)
  # transitions: multinomial logit
  # normal priors on logit of all but one transition probs
  for (i in 1:2){
    lpsiA[i] ~ dnorm(0, sd = 1000)
    lpsiB[i] ~ dnorm(0, sd = 1000)
    lpsiC[i] ~ dnorm(0, sd = 1000)
  }
  # constrain the transitions such that their sum is < 1
  for (i in 1:2){
    psiA[i] <- exp(lpsiA[i]) / (1 + exp(lpsiA[1]) + exp(lpsiA[2]))
    psiB[i] <- exp(lpsiB[i]) / (1 + exp(lpsiB[1]) + exp(lpsiB[2]))
    psiC[i] <- exp(lpsiC[i]) / (1 + exp(lpsiC[1]) + exp(lpsiC[2]))
  }
  # last transition probability
  psiA[3] <- 1 - psiA[1] - psiA[2]
  psiB[3] <- 1 - psiB[1] - psiB[2]
  psiC[3] <- 1 - psiC[1] - psiC[2]

  # probabilities of state z(t+1) given z(t)
  gamma[1,1] <- phiA * psiA[1]
  gamma[1,2] <- phiA * psiA[2]
  gamma[1,3] <- phiA * psiA[3]
  gamma[1,4] <- 1 - phiA
  gamma[2,1] <- phiB * psiB[1]
  gamma[2,2] <- phiB * psiB[2]
  gamma[2,3] <- phiB * psiB[3]
  gamma[2,4] <- 1 - phiB
  gamma[3,1] <- phiC * psiC[1]
  gamma[3,2] <- phiC * psiC[2]
  gamma[3,3] <- phiC * psiC[3]
  gamma[3,4] <- 1 - phiC
  gamma[4,1] <- 0
  gamma[4,2] <- 0
  gamma[4,3] <- 0
  gamma[4,4] <- 1
  
  # probabilities of y(t) given z(t)
  omega[1,1] <- 1 - pA     # Pr(alive A t -> non-detected t)
  omega[1,2] <- pA         # Pr(alive A t -> detected A t)
  omega[1,3] <- 0          # Pr(alive A t -> detected B t)
  omega[1,4] <- 0          # Pr(alive A t -> detected C t)
  omega[2,1] <- 1 - pB     # Pr(alive B t -> non-detected t)
  omega[2,2] <- 0          # Pr(alive B t -> detected A t)
  omega[2,3] <- pB         # Pr(alive B t -> detected B t)
  omega[2,4] <- 0          # Pr(alive B t -> detected C t)
  omega[3,1] <- 1 - pC     # Pr(alive C t -> non-detected t)
  omega[3,2] <- 0          # Pr(alive C t -> detected A t)
  omega[3,3] <- 0          # Pr(alive C t -> detected B t)
  omega[3,4] <- pC         # Pr(alive C t -> detected C t)
  omega[4,1] <- 1          # Pr(dead t -> non-detected t)
  omega[4,2] <- 0          # Pr(dead t -> detected A t)
  omega[4,3] <- 0          # Pr(dead t -> detected B t)
  omega[4,4] <- 0          # Pr(dead t -> detected C t)
  
  # likelihood 
  for (i in 1:N){
    # latent state at first capture
    z[i,first[i]] <- y[i,first[i]] - 1
    for (t in (first[i]+1):K){
      # z(t) given z(t-1)
      z[i,t] ~ dcat(gamma[z[i,t-1],1:4])
      # y(t) given z(t)
      y[i,t] ~ dcat(omega[z[i,t],1:4])
    }
  }
})

#' 
#' List of data. 
## ---------------------------------------------------------------------------------------------------------
my.data <- list(y = y + 1)

#' 
#' List of constants.
## ---------------------------------------------------------------------------------------------------------
my.constants <- list(first = first, 
                     K = ncol(y), 
                     N = nrow(y))

#' 
#' Initial values. 
## ---------------------------------------------------------------------------------------------------------
zinits <- y
zinits[zinits==0] <- sample(c(1,2,3), sum(zinits==0), replace = TRUE)
initial.values <- function(){list(phiA = runif(1, 0, 1), 
                                  phiB = runif(1, 0, 1), 
                                  phiC = runif(1, 0, 1), 
                                  lpsiA = rnorm(2, 0, 10),
                                  lpsiB = rnorm(2, 0, 10),
                                  lpsiC = rnorm(2, 0, 10),
                                  pA = runif(1, 0, 1), 
                                  pB = runif(1, 0, 1), 
                                  pC = runif(1, 0, 1),
                                  z = zinits)}

#' 
#' Parameters to monitor.
## ---------------------------------------------------------------------------------------------------------
parameters.to.save <- c("phiA", "phiB", "phiC", "psiA", "psiB", "psiC","pA", "pB", "pC")
parameters.to.save

#' 
#' Run nimble. 
## ---------------------------------------------------------------------------------------------------------
mcmc.multisite <- nimbleMCMC(code = multisite, 
                             constants = my.constants,
                             data = my.data,              
                             inits = initial.values,
                             monitors = parameters.to.save,
                             niter = n.iter,
                             nburnin = n.burnin, 
                             nchains = n.chains)

#' 
#' Inspect the results. 
## ---------------------------------------------------------------------------------------------------------
MCMCsummary(mcmc.multisite, round = 2)

#' 
#' Caterpillar plot. 
## ---- echo = TRUE, fig.width = 7.5, fig.asp = 0.618, dev = "svg"------------------------------------------
MCMCplot(mcmc.multisite)

#' 
#' 
#' ### Dirichlet prior
#' 
#' Another method is to use Dirichlet priors for the movement parameters. Let's write the model. 
## ---------------------------------------------------------------------------------------------------------
multisite <- nimbleCode({
  
  # -------------------------------------------------
  # Parameters:
  # phiA: survival probability site A
  # phiB: survival probability site B
  # phiC: survival probability site B
  # psiAA = psiA[1]: movement probability from site A to site A (reference)
  # psiAB = psiA[2]: movement probability from site A to site B
  # psiAC = psiA[3]: movement probability from site A to site C 
  # psiBA = psiB[1]: movement probability from site B to site A
  # psiBB = psiB[2]: movement probability from site B to site B (reference)
  # psiBC = psiB[3]: movement probability from site B to site C
  # psiCA = psiC[1]: movement probability from site C to site A
  # psiCB = psiC[2]: movement probability from site C to site B
  # psiCC = psiC[3]: movement probability from site C to site C (reference)
  # pA: recapture probability site A
  # pB: recapture probability site B
  # pC: recapture probability site C
  # -------------------------------------------------
  # States (z):
  # 1 alive at A
  # 2 alive at B
  # 2 alive at C
  # 3 dead
  # Observations (y):  
  # 1 not seen
  # 2 seen at A 
  # 3 seen at B
  # 3 seen at C
  # -------------------------------------------------
  
  # Priors
  phiA ~ dunif(0, 1)
  phiB ~ dunif(0, 1)
  phiC ~ dunif(0, 1)
  pA ~ dunif(0, 1)
  pB ~ dunif(0, 1)
  pC ~ dunif(0, 1)
  # transitions: Dirichlet priors
  psiA[1:3] ~ ddirch(alpha[1:3])
  psiB[1:3] ~ ddirch(alpha[1:3])
  psiC[1:3] ~ ddirch(alpha[1:3])

  # probabilities of state z(t+1) given z(t)
  gamma[1,1] <- phiA * psiA[1]
  gamma[1,2] <- phiA * psiA[2]
  gamma[1,3] <- phiA * psiA[3]
  gamma[1,4] <- 1 - phiA
  gamma[2,1] <- phiB * psiB[1]
  gamma[2,2] <- phiB * psiB[2]
  gamma[2,3] <- phiB * psiB[3]
  gamma[2,4] <- 1 - phiB
  gamma[3,1] <- phiC * psiC[1]
  gamma[3,2] <- phiC * psiC[2]
  gamma[3,3] <- phiC * psiC[3]
  gamma[3,4] <- 1 - phiC
  gamma[4,1] <- 0
  gamma[4,2] <- 0
  gamma[4,3] <- 0
  gamma[4,4] <- 1
  
  # probabilities of y(t) given z(t)
  omega[1,1] <- 1 - pA     # Pr(alive A t -> non-detected t)
  omega[1,2] <- pA         # Pr(alive A t -> detected A t)
  omega[1,3] <- 0          # Pr(alive A t -> detected B t)
  omega[1,4] <- 0          # Pr(alive A t -> detected C t)
  omega[2,1] <- 1 - pB     # Pr(alive B t -> non-detected t)
  omega[2,2] <- 0          # Pr(alive B t -> detected A t)
  omega[2,3] <- pB         # Pr(alive B t -> detected B t)
  omega[2,4] <- 0          # Pr(alive B t -> detected C t)
  omega[3,1] <- 1 - pC     # Pr(alive C t -> non-detected t)
  omega[3,2] <- 0          # Pr(alive C t -> detected A t)
  omega[3,3] <- 0          # Pr(alive C t -> detected B t)
  omega[3,4] <- pC         # Pr(alive C t -> detected C t)
  omega[4,1] <- 1          # Pr(dead t -> non-detected t)
  omega[4,2] <- 0          # Pr(dead t -> detected A t)
  omega[4,3] <- 0          # Pr(dead t -> detected B t)
  omega[4,4] <- 0          # Pr(dead t -> detected C t)
  
  # likelihood 
  for (i in 1:N){
    # latent state at first capture
    z[i,first[i]] <- y[i,first[i]] - 1
    for (t in (first[i]+1):K){
      # z(t) given z(t-1)
      z[i,t] ~ dcat(gamma[z[i,t-1],1:4])
      # y(t) given z(t)
      y[i,t] ~ dcat(omega[z[i,t],1:4])
    }
  }
})

#' 
#' Data in a list. 
## ---------------------------------------------------------------------------------------------------------
my.data <- list(y = y + 1)

#' 
#' Constants in a list. 
## ---------------------------------------------------------------------------------------------------------
my.constants <- list(first = first, 
                     K = ncol(y), 
                     N = nrow(y),
                     alpha = c(1, 1, 1))

#' 
#' Initial values. 
## ---------------------------------------------------------------------------------------------------------
zinits <- y
zinits[zinits==0] <- sample(c(1,2,3), sum(zinits==0), replace = TRUE)
initial.values <- function(){list(phiA = runif(1, 0, 1), 
                                  phiB = runif(1, 0, 1), 
                                  phiC = runif(1, 0, 1), 
                                  psiA = rep(1/3, 3),
                                  psiB = rep(1/3, 3),
                                  psiC = rep(1/3, 3),
                                  pA = runif(1, 0, 1), 
                                  pB = runif(1, 0, 1), 
                                  pC = runif(1, 0, 1),
                                  z = zinits)}

#' 
#' Run nimble.
## ---------------------------------------------------------------------------------------------------------
mcmc.multisite <- nimbleMCMC(code = multisite, 
                             constants = my.constants,
                             data = my.data,              
                             inits = initial.values,
                             monitors = parameters.to.save,
                             niter = n.iter,
                             nburnin = n.burnin, 
                             nchains = n.chains)

#' 
#' Inspect results. Compare to the estimates we obtained with the multinomial logit link method. 
## ---------------------------------------------------------------------------------------------------------
MCMCsummary(mcmc.multisite, round = 2)

#' 
#' Caterpillar plots. 
## ---- echo = TRUE, fig.width = 7.5, fig.asp = 0.618, dev = "svg"------------------------------------------
MCMCplot(mcmc.multisite)

#' 
#' 
#' ## Multistate model
#' 
#' Now we illustrate how one can replace site by breeding state and address another question in evolution ecology: life-history trade-offs. We use data on Soory shearwaters that were kindly provided by David Fletcher. Birds are seen as breeders or non-breeders. 
## ---------------------------------------------------------------------------------------------------------
titis <- read_csv2("titis.csv", col_names = FALSE)
y <-  as.matrix(titis)
head(y)

#' 
#' The code is very similar to that for the model above with two sites, i.e. the structure is the same, but the transition probabilities have different interpretations. 
## ---------------------------------------------------------------------------------------------------------
multistate <- nimbleCode({
  
  # -------------------------------------------------
  # Parameters:
  # phiB: survival probability state B
  # phiNB: survival probability state NB
  # psiBNB: transition probability from B to NB
  # psiNBB: transition probability from NB to B
  # pB: recapture probability B
  # pNB: recapture probability NB
  # -------------------------------------------------
  # States (z):
  # 1 alive B
  # 2 alive NB
  # 3 dead
  # Observations (y):  
  # 1 not seen
  # 2 seen as B 
  # 3 seen as NB
  # -------------------------------------------------
  
  # priors
  phiB ~ dunif(0, 1)
  phiNB ~ dunif(0, 1)
  psiBNB ~ dunif(0, 1)
  psiNBB ~ dunif(0, 1)
  pB ~ dunif(0, 1)
  pNB ~ dunif(0, 1)
  
  # probabilities of state z(t+1) given z(t)
  gamma[1,1] <- phiB * (1 - psiBNB)
  gamma[1,2] <- phiB * psiBNB
  gamma[1,3] <- 1 - phiB
  gamma[2,1] <- phiNB * psiNBB
  gamma[2,2] <- phiNB * (1 - psiNBB)
  gamma[2,3] <- 1 - phiNB
  gamma[3,1] <- 0
  gamma[3,2] <- 0
  gamma[3,3] <- 1
  
  # probabilities of y(t) given z(t)
  omega[1,1] <- 1 - pB    # Pr(alive B t -> non-detected t)
  omega[1,2] <- pB        # Pr(alive B t -> detected B t)
  omega[1,3] <- 0         # Pr(alive B t -> detected NB t)
  omega[2,1] <- 1 - pNB   # Pr(alive NB t -> non-detected t)
  omega[2,2] <- 0         # Pr(alive NB t -> detected B t)
  omega[2,3] <- pNB       # Pr(alive NB t -> detected NB t)
  omega[3,1] <- 1         # Pr(dead t -> non-detected t)
  omega[3,2] <- 0         # Pr(dead t -> detected N t)
  omega[3,3] <- 0         # Pr(dead t -> detected NB t)
  
  # likelihood 
  for (i in 1:N){
    # latent state at first capture
    z[i,first[i]] <- y[i,first[i]] - 1
    for (t in (first[i]+1):K){
      # z(t) given z(t-1)
      z[i,t] ~ dcat(gamma[z[i,t-1],1:3])
      # y(t) given z(t)
      y[i,t] ~ dcat(omega[z[i,t],1:3])
    }
  }
})

#' 
#' Data in a list. 
## ---------------------------------------------------------------------------------------------------------
my.data <- list(y = y + 1)

#' 
#' Get data of first capture. 
## ---------------------------------------------------------------------------------------------------------
get.first <- function(x) min(which(x != 0))
first <- apply(y, 1, get.first)

#' 
#' Constants in a list. 
## ---------------------------------------------------------------------------------------------------------
my.constants <- list(first = first, 
                     K = ncol(y), 
                     N = nrow(y))

#' 
#' Initial values. Same as before, we initialize the latent states with the detections and non-detections. Then non-detections are replaced by a 1 or a 2 for breeder on non-breeder. 
## ---------------------------------------------------------------------------------------------------------
zinits <- y
zinits[zinits == 0] <- sample(c(1,2), sum(zinits == 0), replace = TRUE)
initial.values <- function(){list(phiNB = runif(1, 0, 1), 
                                  phiB = runif(1, 0, 1), 
                                  psiNBB = runif(1, 0, 1), 
                                  psiBNB = runif(1, 0, 1), 
                                  pNB = runif(1, 0, 1), 
                                  pB = runif(1, 0, 1), 
                                  z = zinits)}

#' 
#' Parameters to be monitored. 
## ---------------------------------------------------------------------------------------------------------
parameters.to.save <- c("phiNB", "phiB","psiNBB", "psiBNB", "pNB", "pB")
parameters.to.save

#' 
#' MCMC details. 
## ---------------------------------------------------------------------------------------------------------
n.iter <- 5000
n.burnin <- 2500
n.chains <- 2

#' 
#' Run nimble.
## ---------------------------------------------------------------------------------------------------------
mcmc.multistate <- nimbleMCMC(code = multistate, 
                             constants = my.constants,
                             data = my.data,              
                             inits = initial.values,
                             monitors = parameters.to.save,
                             niter = n.iter,
                             nburnin = n.burnin, 
                             nchains = n.chains)

#' 
#' Let's inspect the results. 
## ---------------------------------------------------------------------------------------------------------
MCMCsummary(mcmc.multistate, round = 2)

#' 
#' Caterpillar plot of the parameter estimates. 
## ---- echo = TRUE, fig.width = 7.5, fig.asp = 0.618, dev = "svg"------------------------------------------
MCMCplot(mcmc.multistate)

#' 
#' Non-breeder individuals seem to have a survival higher than breeder individuals, suggesting a trade-off between reproduction and survival. Let's compare graphically the survival of breeder and non-breeder individuals. First we gather the values generated for $\phi_B$ and $\phi_{NB}$ for the two chains. 
## ---------------------------------------------------------------------------------------------------------
phiB <- c(mcmc.multistate$chain1[,"phiB"], mcmc.multistate$chain2[,"phiB"])
phiNB <- c(mcmc.multistate$chain1[,"phiNB"], mcmc.multistate$chain2[,"phiNB"])
df <- data.frame(param = c(rep("phiB", length(phiB)), rep("phiNB", length(phiB))), 
                 value = c(phiB, phiNB))
head(df)

#' 
#' Then, we plot a histogram of the two posterior distributions. 
## ---- echo = TRUE, fig.width = 7.5, fig.asp = 0.618, dev = "svg"------------------------------------------
df %>%
  ggplot(aes(x = value, fill = param)) +
  geom_density(color = "white", alpha = 0.6, position = 'identity') +
  scale_fill_manual(values = c("#69b3a2", "#404080")) +
  labs(fill = "", x = "survival")

#' 
#' To formally test this difference, one could fit a model with no state effect on survival and compare the two models with WAIC. 
#' 
#' <!-- knitr::purl(here::here("worksheets","5_demo.Rmd"), documentation = 2) -->
