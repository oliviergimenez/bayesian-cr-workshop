#' ---
#' title: "Class 7 live demo: Uncertainty in state assignment"
#' author: "The team"
#' date: "last updated: `r Sys.Date()`"
#' output: html_document
#' ---
#' 
#' 
## ----setup, include=FALSE, echo=FALSE------------------------------------------------------
options(htmltools.dir.version = FALSE)
knitr::opts_chunk$set(comment = "", message = TRUE, warning = TRUE)
library(tidyverse)
theme_set(theme_light(base_size = 14))
update_geom_defaults("point", list(size = 2)) 
library(here)
library(nimble)

#' 
#' ## Introduction
#' 
#' In this demo, we illustrate how to fit models with several states in which the assignment of individuals to states is difficult to make with certainty. We start by breeding states then continue with disease states and end by illustrating how to incorporate individual heterogeneity in capture-recapture models. 
#' 
## ------------------------------------------------------------------------------------------
library(tidyverse)
library(nimble)
library(MCMCvis)

#' 
#' ## Breeding state uncertainty
#' 
#' Let's get back to the analysis of the Sooty shearwater data. For some individuals, we have some uncertainty in assigning individuals in the breeding or non-breeding state (code 3 in the data). 
## ------------------------------------------------------------------------------------------
titis <- read_csv2("titis_with_uncertainty.csv", col_names = FALSE)
y <- as.matrix(titis)
head(y)

#' 
#' Let's write down the model code. 
## ------------------------------------------------------------------------------------------
multievent <- nimbleCode({
  
  # -------------------------------------------------
  # Parameters:
  # phiB: survival probability state B
  # phiNB: survival probability state NB
  # psiBNB: transition probability from B to NB
  # psiNBB: transition probability from NB to B
  # pB: recapture probability B
  # pNB: recapture probability NB
  # piB prob. of being in initial state breeder
  # betaNB prob to ascertain the breeding status of an individual encountered as non-breeder
  # betaB prob to ascertain the breeding status of an individual encountered as breeder
  # -------------------------------------------------
  # States (z):
  # 1 alive B
  # 2 alive NB
  # 3 dead
  # Observations (y):  
  # 1 = non-detected
  # 2 = seen and ascertained as breeder
  # 3 = seen and ascertained as non-breeder
  # 4 = not ascertained
  # -------------------------------------------------
  
  # priors
  phiB ~ dunif(0, 1)
  phiNB ~ dunif(0, 1)
  psiBNB ~ dunif(0, 1)
  psiNBB ~ dunif(0, 1)
  pB ~ dunif(0, 1)
  pNB ~ dunif(0, 1)
  piB ~ dunif(0, 1)
  betaNB ~ dunif(0, 1)
  betaB ~ dunif(0, 1)
  
  # vector of initial stats probs
  delta[1] <- piB      # prob. of being in initial state B
  delta[2] <- 1 - piB  # prob. of being in initial state NB
  delta[3] <- 0        # prob. of being in initial state dead
  
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
  omega[1,1] <- 1 - pB             # Pr(alive B t -> non-detected t)
  omega[1,2] <- pB * betaB         # Pr(alive B t -> detected B t)
  omega[1,3] <- 0                  # Pr(alive B t -> detected NB t)
  omega[1,4] <- pB * (1 - betaB)   # Pr(alive B t -> detected U t)
  omega[2,1] <- 1 - pNB            # Pr(alive NB t -> non-detected t)
  omega[2,2] <- 0                  # Pr(alive NB t -> detected B t)
  omega[2,3] <- pNB * betaNB       # Pr(alive NB t -> detected NB t)
  omega[2,4] <- pNB * (1 - betaNB) # Pr(alive NB t -> detected U t)
  omega[3,1] <- 1                  # Pr(dead t -> non-detected t)
  omega[3,2] <- 0                  # Pr(dead t -> detected N t)
  omega[3,3] <- 0                  # Pr(dead t -> detected NB t)
  omega[3,4] <- 0                  # Pr(dead t -> detected U t)
  
  omega.init[1,1] <- 0                  # Pr(alive B t = 1 -> non-detected t = 1)
  omega.init[1,2] <- betaB              # Pr(alive B t = 1 -> detected B t = 1)
  omega.init[1,3] <- 0                  # Pr(alive B t = 1 -> detected NB t = 1)
  omega.init[1,4] <- 1 - betaB          # Pr(alive B t = 1 -> detected U t = 1)
  omega.init[2,1] <- 0                  # Pr(alive NB t = 1 -> non-detected t = 1)
  omega.init[2,2] <- 0                  # Pr(alive NB t = 1 -> detected B t = 1)
  omega.init[2,3] <- betaNB             # Pr(alive NB t = 1 -> detected NB t = 1)
  omega.init[2,4] <- 1 - betaNB         # Pr(alive NB t = 1 -> detected U t = 1)
  omega.init[3,1] <- 1                  # Pr(dead t = 1 -> non-detected t = 1)
  omega.init[3,2] <- 0                  # Pr(dead t = 1 -> detected N t = 1)
  omega.init[3,3] <- 0                  # Pr(dead t = 1 -> detected NB t = 1)
  omega.init[3,4] <- 0                  # Pr(dead t = 1 -> detected U t = 1)
  
  # likelihood 
  for (i in 1:N){
    # latent state at first capture
    z[i,first[i]] ~ dcat(delta[1:3])
    y[i,first[i]] ~ dcat(omega.init[z[i,first[i]],1:4])
    for (t in (first[i]+1):K){
      # z(t) given z(t-1)
      z[i,t] ~ dcat(gamma[z[i,t-1],1:3])
      # y(t) given z(t)
      y[i,t] ~ dcat(omega[z[i,t],1:4])
    }
  }
})

#' 
#' Get the data of first capture. 
## ------------------------------------------------------------------------------------------
get.first <- function(x) min(which(x != 0))
first <- apply(y, 1, get.first)

#' 
#' Constants in a list. 
## ------------------------------------------------------------------------------------------
my.constants <- list(first = first, 
                     K = ncol(y), 
                     N = nrow(y))

#' 
#' Data in a list. 
## ------------------------------------------------------------------------------------------
my.data <- list(y = y + 1)

#' 
#' Initial values. 
## ------------------------------------------------------------------------------------------
zinit <- y
zinit[zinit==3] <- sample(c(1,2), sum(zinit==3), replace = TRUE)
for (i in 1:nrow(y)) {
  for (j in 1:ncol(y)) {
    if (j > first[i] & y[i,j]==0) {zinit[i,j] <- which(rmultinom(1, 1, c(1/2,1/2))==1)}
    if (j < first[i]) {zinit[i,j] <- 0}
  }
}
zinit <- as.matrix(zinit)
initial.values <- function(){list(phiNB = runif(1, 0, 1), 
                                  phiB = runif(1, 0, 1), 
                                  psiNBB = runif(1, 0, 1), 
                                  psiBNB = runif(1, 0, 1), 
                                  pNB = runif(1, 0, 1), 
                                  pB = runif(1, 0, 1),
                                  piB = runif(1, 0, 1),
                                  betaB = runif(1, 0, 1),
                                  betaNB = runif(1, 0, 1),
                                  z = zinit)}

#' 
#' Parameters to be monitored. 
## ------------------------------------------------------------------------------------------
parameters.to.save <- c("phiB", 
                        "phiNB", 
                        "psiNBB", 
                        "psiBNB", 
                        "piB", 
                        "pB", 
                        "pNB", 
                        "betaNB", 
                        "betaB")

#' 
#' MCMC details. 
## ------------------------------------------------------------------------------------------
n.iter <- 5000
n.burnin <- 2500
n.chains <- 2

#' 
#' Run nimble. 
## ------------------------------------------------------------------------------------------
mcmc.multievent <- nimbleMCMC(code = multievent, 
                              constants = my.constants,
                              data = my.data,              
                              inits = initial.values,
                              monitors = parameters.to.save,
                              niter = n.iter,
                              nburnin = n.burnin, 
                              nchains = n.chains)

#' 
#' Let's have a look to the results. Breeders were mostly assigned as uncertain, while non-breeders were positively assigned.
## ------------------------------------------------------------------------------------------
MCMCsummary(mcmc.multievent, round = 2)

#' 
#' Caterpillar plot of the parameter estimates. 
## ---- echo = TRUE, fig.width = 7.5, fig.asp = 0.618, dev = "svg"---------------------------
MCMCplot(mcmc.multievent)

#' 
#' 
#' ## Disease state uncertainty
#' 
#' We now turn to the analysis of data from a study of the dynamics of conjunctivitis in the house finch (Carpodacus mexicanus Müller). A challenge in modeling disease dynamics involves the estimation of recovery and infection rates, which can be complicated when disease state is uncertain (seen at distance in the house finch case study). In the data we have 0 for a non-detection, 1 for a detection of a healthy individual, 2 for the detection of a sick individual and 3 for an individual for which we cannot say whether it is sane or ill. 
## ------------------------------------------------------------------------------------------
y <- as.matrix(read.table("hofi_2000.txt"))
head(y)

#' 
#' Write the model. The model is similar to the model of the previous section. 
## ------------------------------------------------------------------------------------------
hmm.disease <- nimbleCode({
  
  # priors
  phiH ~ dunif(0, 1)   # prior survival healthy
  phiI ~ dunif(0, 1)   # prior survival ill
  psiIH ~ dunif(0, 1)  # prior transition ill -> healthy
  psiHI ~ dunif(0, 1)  # prior transition healthy -> ill
  pH ~ dunif(0, 1)     # prior detection healthy
  pI ~ dunif(0, 1)     # prior detection ill
  pi ~ dunif(0, 1)     # prob init state healthy
  betaH ~ dunif(0,1)
  betaI ~ dunif(0,1)

  # HMM ingredients
  delta[1] <- pi         # Pr(healthy t = 1) = pi
  delta[2] <- 1 - pi     # Pr(ill t = 1) = 1 - pi
  delta[3] <- 0          # Pr(dead t = 1) = 0
  
  gamma[1,1] <- phiH * (1 - psiHI)      # Pr(H t -> H t+1)
  gamma[1,2] <- phiH * psiHI            # Pr(H t -> I t+1)
  gamma[1,3] <- 1 - phiH                # Pr(alive t -> dead t+1)
  gamma[2,1] <- phiI * psiIH            # Pr(I t -> H t+1)
  gamma[2,2] <- phiI * (1 - psiIH)      # Pr(I t -> I t+1)
  gamma[2,3] <- 1 - phiI                # Pr(alive t -> dead t+1)
  gamma[3,1] <- 0                       # Pr(dead t -> alive t+1)
  gamma[3,2] <- 0                       # Pr(dead t -> alive t+1)
  gamma[3,3] <- 1                       # Pr(dead t -> dead t+1)
  
  omega[1,1] <- 1 - pH           # Pr(H t -> non-detected t)
  omega[1,2] <- pH * betaH       # Pr(H t -> detected H t)
  omega[1,3] <- 0                # Pr(H t -> detected I t)
  omega[1,4] <- pH * (1 - betaH) # Pr(H t -> detected U t)
  omega[2,1] <- 1 - pI           # Pr(I t -> non-detected t)
  omega[2,2] <- 0                # Pr(I t -> detected H t)
  omega[2,3] <- pI * betaI       # Pr(I t -> detected I t)
  omega[2,4] <- pI * (1 - betaI) # Pr(I t -> detected U t)
  omega[3,1] <- 1                # Pr(dead t -> non-detected t)
  omega[3,2] <- 0                # Pr(dead t -> detected H t)
  omega[3,3] <- 0                # Pr(dead t -> detected I t)
  omega[3,4] <- 0                # Pr(dead t -> detected U t)
  
  omega.init[1,1] <- 0                # Pr(H t = 1 -> non-detected t = 1)
  omega.init[1,2] <- betaH            # Pr(H t = 1 -> detected H t = 1)
  omega.init[1,3] <- 0                # Pr(H t = 1 -> detected I t = 1)
  omega.init[1,4] <- 1 - betaH        # Pr(H t = 1 -> detected U t = 1)
  omega.init[2,1] <- 0                # Pr(I t = 1 -> non-detected t = 1)
  omega.init[2,2] <- 0                # Pr(I t = 1 -> detected H t = 1)
  omega.init[2,3] <- betaI            # Pr(I t = 1 -> detected I t = 1)
  omega.init[2,4] <- 1 - betaI        # Pr(I t = 1 -> detected U t = 1)
  omega.init[3,1] <- 1                # Pr(dead t = 1 -> non-detected t = 1)
  omega.init[3,2] <- 0                # Pr(dead t = 1 -> detected H t = 1)
  omega.init[3,3] <- 0                # Pr(dead t = 1 -> detected I t = 1)
  omega.init[3,4] <- 0                # Pr(dead t = 1 -> detected U t = 1)

  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:3])
    y[i,first[i]] ~ dcat(omega.init[z[i,first[i]], 1:4])
    for (j in (first[i]+1):K){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:3])
      y[i,j] ~ dcat(omega[z[i,j], 1:4])
    }
  }
})

#' 
#' Get the date of first capture. 
## ------------------------------------------------------------------------------------------
first <- apply(y, 1, function(x) min(which(x !=0)))

#' 
#' Constants in a list. 
## ------------------------------------------------------------------------------------------
my.constants <- list(N = nrow(y), K = ncol(y), first = first)

#' 
#' Data in a list. 
## ------------------------------------------------------------------------------------------
my.data <- list(y = y + 1)

#' 
#' Initial values. 
## ------------------------------------------------------------------------------------------
zinit <- y
zinit[zinit==3] <- sample(c(1,2), sum(zinit==3), replace = TRUE)
for (i in 1:nrow(y)) {
  for (j in 1:ncol(y)) {
    if (j > first[i] & y[i,j]==0) {zinit[i,j] <- which(rmultinom(1, 1, c(1/2,1/2))==1)}
    if (j < first[i]) {zinit[i,j] <- 0}
  }
}
zinit <- as.matrix(zinit)
initial.values <- function() list(phiH = runif(1, 0, 1),
                                  phiI = runif(1, 0, 1),
                                  pH = runif(1, 0, 1),
                                  pI = runif(1, 0, 1),
                                  pi = runif(1, 0, 1),
                                  betaH = runif(1, 0, 1),
                                  betaI = runif(1, 0, 1),
                                  psiHI = runif(1, 0, 1),
                                  psiIH = runif(1, 0, 1),
                                  z = zinit)

#' 
#' Parameters to be monitored. 
## ------------------------------------------------------------------------------------------
parameters.to.save <- c("phiH", 
                        "phiI", 
                        "pH", 
                        "pI", 
                        "pi",
                        "betaH", 
                        "betaI", 
                        "psiHI", 
                        "psiIH")

#' 
#' MCMC details. 
## ------------------------------------------------------------------------------------------
n.iter <- 40000
n.burnin <- 25000
n.chains <- 2

#' 
#' Run nimble. 
## ----eval = FALSE--------------------------------------------------------------------------
out <- nimbleMCMC(code = hmm.disease, 
                           constants = my.constants,
                           data = my.data,              
                           inits = initial.values,
                           monitors = parameters.to.save,
                           niter = n.iter,
                           nburnin = n.burnin, 
                           nchains = n.chains)
save(out, file = "house_finch.RData")

#' 
## ----echo=FALSE----------------------------------------------------------------------------
load("house_finch.RData")

#' 
#' Inspect the results. 
## ------------------------------------------------------------------------------------------
MCMCsummary(out, round = 2)

#' 
#' Let's have a look to the trace plots, and the estimated posterior distribution of the detection probabilities in healthy and ill states. 
## ---- echo = TRUE, fig.width = 7.5, fig.asp = 0.618, dev = "svg"---------------------------
MCMCtrace(out, pdf = FALSE, params = c("pH", "pI"))

#' 
#' What about the assignment probabilities? The non‐infected birds were positively identified, while most infected encounters were classified as unknown.
## ---- echo = TRUE, fig.width = 7.5, fig.asp = 0.618, dev = "svg"---------------------------
MCMCtrace(out, pdf = FALSE, params = c("betaH", "betaI"))

#' 
#' And survival?
## ---- echo = TRUE, fig.width = 7.5, fig.asp = 0.618, dev = "svg"---------------------------
MCMCtrace(out, pdf = FALSE, params = c("phiH", "phiI"))

#' 
#' What about the infection and recovery rates? 
## ---- echo = TRUE, fig.width = 7.5, fig.asp = 0.618, dev = "svg"---------------------------
MCMCtrace(out, pdf = FALSE, params = c("psiHI", "psiIH"))

#' 
#' Last, the proportion of new marked individuals that are in the healthy state. 
## ---- echo = TRUE, fig.width = 7.5, fig.asp = 0.618, dev = "svg"---------------------------
MCMCtrace(out, pdf = FALSE, params = c("pi"))

#' 
#' 
#' ## Individual heterogeneity
#' 
#' For our last example, we will see how to incorporate individual heterogeneity in capture-recapture models using finite mixtures. We will use a case study on gray wolves. The data were kindly provided by the French Office for Biodiversity. 
## ------------------------------------------------------------------------------------------
y <- as.matrix(read.table("wolf.txt"))
# 1 seen
# 0 not seen

#' 
#' Let's write the model. Only thing to notice is that we have two alive states to represent two classes of individuals. We choose to have the detection probabilities dependent on the states and constant survival. 
## ------------------------------------------------------------------------------------------
hmm.phipmix <- nimbleCode({
  
  # priors
  phi ~ dunif(0, 1) # prior survival
  p1 ~ dunif(0, 1) # prior detection
  p2 ~ dunif(0, 1) # prior detection
  pi ~ dunif(0, 1) # prob init state 1

  # HMM ingredients
  gamma[1,1] <- phi      # Pr(alive 1 t -> alive 1 t+1)
  gamma[1,2] <- 0        # Pr(alive 1 t -> alive 2 t+1) // no transition
  gamma[1,3] <- 1 - phi  # Pr(alive t -> dead t+1)
  gamma[2,1] <- 0        # Pr(alive 1 t -> alive 1 t+1) // no transition
  gamma[2,2] <- phi      # Pr(alive 1 t -> alive 2 t+1) 
  gamma[2,3] <- 1 - phi  # Pr(alive t -> dead t+1)
  gamma[3,1] <- 0        # Pr(dead t -> alive 1 t+1)
  gamma[3,2] <- 0        # Pr(dead t -> alive 1 t+1)
  gamma[3,3] <- 1        # Pr(dead t -> dead t+1)
  delta[1] <- pi         # Pr(alive t = 1) = pi
  delta[2] <- 1 - pi     # Pr(alive t = 1) = 1 - pi
  delta[3] <- 0          # Pr(dead t = 1) = 0
  omega[1,1] <- 1 - p1   # Pr(alive state 1 t -> non-detected t)
  omega[1,2] <- p1       # Pr(alive state 1 t -> detected t)
  omega[2,1] <- 1 - p2   # Pr(alive state 2 t -> non-detected t)
  omega[2,2] <- p2       # Pr(alive state 2 t -> detected t)
  omega[3,1] <- 1        # Pr(dead t -> non-detected t)
  omega[3,2] <- 0        # Pr(dead t -> detected t)
  omega.init[1,1] <- 0   # Pr(alive state 1 t -> non-detected t)
  omega.init[1,2] <- 1   # Pr(alive state 1 t -> detected t)
  omega.init[2,1] <- 0   # Pr(alive state 2 t -> non-detected t)
  omega.init[2,2] <- 1   # Pr(alive state 2 t -> detected t)
  omega.init[3,1] <- 1   # Pr(dead t -> non-detected t)
  omega.init[3,2] <- 0   # Pr(dead t -> detected t)
  
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:3])
    y[i,first[i]] ~ dcat(omega.init[z[i,first[i]], 1:2])
    for (j in (first[i]+1):K){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:3])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})

#' 
#' Get the date of first capture. 
## ------------------------------------------------------------------------------------------
first <- apply(y, 1, function(x) min(which(x != 0)))

#' 
#' Constants in a list. 
## ------------------------------------------------------------------------------------------
my.constants <- list(N = nrow(y), 
                     K = ncol(y), 
                     first = first)

#' 
#' Data in a list. 
## ------------------------------------------------------------------------------------------
my.data <- list(y = y + 1)

#' 
#' Initial values. 
## ------------------------------------------------------------------------------------------
zinit <- y
for (i in 1:nrow(y)) {
  for (j in first[i]:ncol(y)) {
    if (j == first[i]) zinit[i,j] <- which(rmultinom(1, 1, c(1/2,1/2))==1) # pick alive state
    if (j > first[i]) zinit[i,j] <- zinit[i,j-1] # because no transition, state remains the same
  }
}
zinit <- as.matrix(zinit)
initial.values <- function() list(phi = runif(1,0,1),
                                  p1 = runif(1,0,1),
                                  p2 = runif(1,0,1),
                                  pi = runif(1,0,1),
                                  z = zinit)

#' 
#' Parameters to be monitored. 
## ------------------------------------------------------------------------------------------
parameters.to.save <- c("phi", "p1", "p2", "pi")

#' 
#' MCMC details. 
## ------------------------------------------------------------------------------------------
n.iter <- 10000
n.burnin <- 5000
n.chains <- 2

#' 
#' Run nimble.
## ------------------------------------------------------------------------------------------
mcmc.phipmix <- nimbleMCMC(code = hmm.phipmix, 
                           constants = my.constants,
                           data = my.data,              
                           inits = initial.values,
                           monitors = parameters.to.save,
                           niter = n.iter,
                           nburnin = n.burnin, 
                           nchains = n.chains)

#' 
#' Numerical summaries. 
## ------------------------------------------------------------------------------------------
MCMCsummary(mcmc.phipmix, round = 2)

#' 
#' Let's have a look to the trace and posterior distribution of the proportion of individuals in each class, as well as the detection probabilities. 
## ---- echo = TRUE, fig.width = 7.5, fig.asp = 0.618, dev = "svg"---------------------------
MCMCtrace(mcmc.phipmix, pdf = FALSE, params = c("pi", "p1", "p2"))

#' 
#' Last, survival. 
## ---- echo = TRUE, fig.width = 7.5, fig.asp = 0.618, dev = "svg"---------------------------
MCMCtrace(mcmc.phipmix, pdf = FALSE, params = c("phi"))

#' 
#' 
#' 
