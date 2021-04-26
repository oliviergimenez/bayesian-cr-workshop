#' ---
#' title: "Class 8 live demo: Skip your coffee break: Speed up MCMC convergence"
#' author: "The team"
#' date: "last updated: `r Sys.Date()`"
#' output: html_document
#' ---
#' 
#' 
## ----setup, include=FALSE, echo=FALSE------------------------------------
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
#' In this demo, we illustrate how nimble can be used to speed up MCMC computations. First, we demonstrate how to change the default samplers. Then we use [the nimbleEcology package](https://cran.r-project.org/web/packages/nimbleEcology/vignettes/Introduction_to_nimbleEcology.html) which implements marginalization (by integrating or sum the likelihood over latent states). Marginalization eliminates the need for MCMC sampling these latent variables, which often provide efficiency gains. We also illustrate how to use nimble functions to write a likelihood that works with a dataset in which we pool individuals with the same encounter histories.
#' 
## ------------------------------------------------------------------------
library(tidyverse)
library(nimbleEcology)
library(MCMCvis)

#' 
#' ## Work with samplers: Use slice sampling
#' 
#' We're going back to the Dipper example. Let's consider a model with wing length and individual random effect on survival.  
#' 
#' 
#' Read in data.
## ----echo = TRUE, message = FALSE, warning=FALSE-------------------------
dipper <- read_csv(here::here("slides", "dat", "dipper.csv"))
y <- dipper %>%
  select(year_1981:year_1987) %>%
  as.matrix()

#' 
#' 
#' Get occasion of first capture. 
## ----echo = TRUE, message = FALSE, warning=FALSE-------------------------
first <- apply(y, 1, function(x) min(which(x !=0)))

#' 
#' 
#' Constants in a list. 
## ----echo = TRUE, message = FALSE, warning=FALSE-------------------------
wing.length.st <- as.vector(scale(dipper$wing_length))
my.constants <- list(N = nrow(y), 
                     T = ncol(y), 
                     first = first,
                     winglength = wing.length.st)

#' 
#' Data in a list.
## ------------------------------------------------------------------------
my.data <- list(y = y + 1)

#' 
#' Initial values. 
## ----echo = TRUE, message = FALSE, warning=FALSE-------------------------
zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(beta = rnorm(2,0,1.5),
                                  sdeps = runif(1,0,3),
                                  p = runif(1,0,1),
                                  z = zinits)

#' 
#' MCMC details. 
## ----echo = TRUE, message = FALSE, warning=FALSE-------------------------
parameters.to.save <- c("beta", "sdeps", "p")
n.iter <- 5000
n.burnin <- 2500
n.chains <- 2

#' 
#' Write the model. 
## ------------------------------------------------------------------------
hmm.phiwlrep <- nimbleCode({
    p ~ dunif(0, 1) # prior detection
    omega[1,1] <- 1 - p    # Pr(alive t -> non-detected t)
    omega[1,2] <- p        # Pr(alive t -> detected t)
    omega[2,1] <- 1        # Pr(dead t -> non-detected t)
    omega[2,2] <- 0        # Pr(dead t -> detected t)
  for (i in 1:N){
    logit(phi[i]) <- beta[1] + beta[2] * winglength[i] + eps[i]
    eps[i] ~ dnorm(mean = 0, sd = sdeps)
    gamma[1,1,i] <- phi[i]      # Pr(alive t -> alive t+1)
    gamma[1,2,i] <- 1 - phi[i]  # Pr(alive t -> dead t+1)
    gamma[2,1,i] <- 0           # Pr(dead t -> alive t+1)
    gamma[2,2,i] <- 1           # Pr(dead t -> dead t+1)
  }
  beta[1] ~ dnorm(mean = 0, sd = 1.5)
  beta[2] ~ dnorm(mean = 0, sd = 1.5)
  sdeps ~ dunif(0, 10)
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, i])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})

#' 
#' Run nimble.
## ----echo = TRUE, eval = TRUE--------------------------------------------
mcmc.phiwlrep <- nimbleMCMC(code = hmm.phiwlrep, 
                            constants = my.constants,
                            data = my.data,              
                            inits = initial.values,
                            monitors = parameters.to.save,
                            niter = n.iter,
                            nburnin = n.burnin, 
                            nchains = n.chains)

#' 
#' 
#' Let's have a look to the trace of the standard deviation of the random effect. Hmmm.  
## ---- echo = TRUE, fig.width = 7.5, fig.asp = 0.718, dev = "svg", message=FALSE, warning=FALSE, fig.align="center"----
MCMCtrace(mcmc.phiwlrep, params = "sdeps", pdf = FALSE)

#' 
#' Let's try and change the default sampler for this parameter. What are the samplers used by default?
## ----echo = TRUE, warning=FALSE, message=FALSE---------------------------
hmm.phiwlrep <- nimbleModel(code = hmm.phiwlrep,
                            constants = my.constants,
                            data = my.data,
                            inits = initial.values())
mcmcConf <- configureMCMC(hmm.phiwlrep)

#' 
#' We remove the default sampler, and use slice sampler instead. 
## ------------------------------------------------------------------------
mcmcConf$removeSamplers('sdeps')
mcmcConf$addSampler(target = 'sdeps',
                    type = "slice")
mcmcConf

#' 
#' Compile model and MCMC.
## ------------------------------------------------------------------------
Rmcmc <- buildMCMC(mcmcConf)
Cmodel <- compileNimble(hmm.phiwlrep)
Cmcmc <- compileNimble(Rmcmc, project=hmm.phiwlrep)

#' 
#' Run nimble. 
## ---- echo = TRUE, fig.width = 7.5, fig.asp = 0.618, dev = "svg", warning=FALSE, message=FALSE, fig.align="center"----
Cmcmc$run(10000)
samples1 <- as.matrix(Cmcmc$mvSamples)
Cmcmc$run(10000)
samples2 <- as.matrix(Cmcmc$mvSamples)

#' 
#' Format results in data.frames. 
## ---- echo = TRUE, fig.width = 7.5, fig.asp = 0.618, dev = "svg", warning=FALSE, message=FALSE, fig.align="center"----
df.sdeps <- data.frame(iter = c(2501:10000, 2501:10000),
                 samples = c(samples1[2501:10000,"sdeps"], samples2[2501:10000,"sdeps"]), 
                 chain = c(rep("chain 1", length(samples1[2501:10000,"sdeps"])),
                           rep("chain 2", length(samples2[2501:10000,"sdeps"]))))
df.beta <- data.frame(iter = c(2501:10000, 2501:10000),
                 beta1 = c(samples1[2501:10000,"beta[1]"], samples2[2501:10000,"beta[1]"]), 
                 beta2 = c(samples1[2501:10000,"beta[2]"], samples2[2501:10000,"beta[2]"]), 
                 chain = c(rep("chain 1", length(samples1[2501:10000,"sdeps"])),
                           rep("chain 2", length(samples2[2501:10000,"sdeps"]))))

#' 
#' Trace plot for the standard deviation of the random effect. 
## ---- echo = FALSE, fig.width = 7.5, fig.asp = 0.618, dev = "svg", warning=FALSE, message=FALSE, fig.align="center"----
df.sdeps %>%
  ggplot() + 
  aes(x = iter, y = samples, group = chain, color = chain) + 
  geom_line() + 
  labs(x = "iterations", y = "random effect standard deviation", color = "")

#' 
#' ## Work with samplers: Use block sampling
#' 
#' High correlation in (regression) parameters may make independent samplers inefficient. Let's have a look to the correlation between the intercept and the slope of the relationship between survival and wing length in the European dipper. 
## ---- echo = FALSE, fig.width = 7.5, fig.asp = 0.418, dev = "svg", fig.align="center"----
df.beta %>%
  ggplot() + 
  aes(x = beta1, y = beta2, group = chain, color = chain) +
  geom_point(alpha = .2) + 
  labs(x = "beta1", y = "beta2", color = "")

#' 
#' There seems to be no correlation. Let's act for a minute as if we had a strong correlation. We're gonna see how block sampling might help. In block sampling, we propose candidate values from a multivariate distribution.
#' 
#' First, we remove and replace independent RW samples by block sampling. Then we proceed as usual.
## ------------------------------------------------------------------------
mcmcConf$removeSamplers(c('beta[1]','beta[2]'))
mcmcConf$addSampler(target = c('beta[1]','beta[2]'),
                    type = "RW_block")
mcmcConf
Rmcmc <- buildMCMC(mcmcConf)
Cmodel <- compileNimble(hmm.phiwlrep)
Cmcmc <- compileNimble(Rmcmc, project=hmm.phiwlrep)

#' 
#' Run nimble and a single chain.
## ------------------------------------------------------------------------
Cmcmc$run(10000)
samples <- as.matrix(Cmcmc$mvSamples)

#' 
#' Summarize.
## ------------------------------------------------------------------------
samples %>%
  as_tibble() %>%
  select(!starts_with("z")) %>% # ignore the latent states z
  summarise(across(everything(), list(mean = mean, sd = sd)))

#' 
#' 
#' ## Marginalization with nimbleEcology
#' 
#' Let's get back to the analysis of the Canada geese data, with 3 sites. 
## ------------------------------------------------------------------------
geese <- read_csv("geese.csv", col_names = TRUE)
y <- as.matrix(geese)

#' 
#' Get the occasion of first capture for each individual.
## ------------------------------------------------------------------------
get.first <- function(x) min(which(x != 0))
first <- apply(y, 1, get.first)

#' 
#' We filter out individuals that are first captured at last occasion. These individuals do not contribute to parameter estimation, and also they cause problems with nimbleEcology. 
## ------------------------------------------------------------------------
mask <- which(first!=ncol(y)) # individuals that are not first encountered at last occasion
y <- y[mask, ]                # keep only these
first <- first[mask]

#' 
## ----echo = FALSE--------------------------------------------------------
dHMM <- nimbleFunction(
  run = function (x = double(1), 
                  init = double(1), 
                  probObs = double(2),
                  probTrans = double(2), 
                  len = double(0),
                  checkRowSums = integer(0, default = 0),  
                  log = integer(0, default = 0)) 
{
    if (length(x) != len) 
        nimStop("In dHMM: Argument len must be length of x or 0.")
    if (nimDim(probObs)[1] != nimDim(probTrans)[1]) 
        nimStop("In dHMM: Length of dimension 1 in probObs must equal length of dimension 1 in probTrans.")
    if (nimDim(probTrans)[1] != nimDim(probTrans)[2]) 
      nimStop("In dHMM: probTrans must be a square matrix.")
    ## There was a strict test for sum(init) == 1.  This could be true in R and false in C++!
    if (abs(sum(init) - 1) > 1e-06) 
        nimStop("In dHMM: Initial probabilities must sum to 1.")
    if (checkRowSums) {
        transCheckPasses <- TRUE
        for (i in 1:nimDim(probTrans)[1]) {
            thisCheckSum <- sum(probTrans[i, ])
            if (abs(thisCheckSum - 1) > 1e-06) {
                nimPrint("In dHMM: Problem with sum(probTrans[i,]) with i = ", 
                  i, ". The sum should be 1 but is ", thisCheckSum)
                transCheckPasses <- FALSE
            }
        }
        obsCheckPasses <- TRUE
        for (i in 1:nimDim(probObs)[1]) {
            thisCheckSum <- sum(probObs[i, ])
            if (abs(thisCheckSum - 1) > 1e-06) {
                nimPrint("In dHMM: Problem with sum(probObs[i,]) with i = ", 
                  i, ". The sum should be 1 but is ", thisCheckSum)
                obsCheckPasses <- FALSE
            }
        }
        if (!(transCheckPasses | obsCheckPasses)) 
            nimStop("In dHMM: probTrans and probObs were not specified correctly.  Probabilities in each row (second dimension) must sum to 1.")
        if (!transCheckPasses) 
            nimStop("In dHMM: probTrans was not specified correctly.  Probabilities in each row (second dimension) must sum to 1.")
        if (!obsCheckPasses) 
            nimStop("In dHMM: probObs was not specified correctly. Probabilities in each row must sum to 1.")
    }
    pi <- init
    logL <- 0
    nObsClasses <- nimDim(probObs)[2]
    for (t in 1:len) {
        if (x[t] > nObsClasses | x[t] < 1) 
            nimStop("In dHMM: Invalid value of x[t].")
        Zpi <- probObs[, x[t]] * pi
        sumZpi <- sum(Zpi)
        logL <- logL + log(sumZpi)
        if (t != len) 
            pi <- ((Zpi %*% probTrans)/sumZpi)[1, ]
    }
    if (log) 
        return(logL)
    return(exp(logL))
    returnType(double())
})

#' 
#' Let's write the model. Note that the likelihood is simpler due to the use of the function `dHMM`.
## ------------------------------------------------------------------------
multisite.marginalized <- nimbleCode({
  
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
  
  # survival priors
  phiA ~ dunif(0, 1)
  phiB ~ dunif(0, 1)
  phiC ~ dunif(0, 1)
  # priors for detection
  pA ~ dunif(0, 1)
  pB ~ dunif(0, 1)
  pC ~ dunif(0, 1)
  # priors for transitions: Dirichlet
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
  # initial state probs
  for(i in 1:N) {
    init[i, 1:4] <- gamma[ y[i, first[i] ] - 1, 1:4 ] # First state propagation
  }
    
  # likelihood 
  for (i in 1:N){
      y[i,(first[i]+1):K] ~ dHMM(init = init[i,1:4],  # count data from first[i] + 1
                                 probObs = omega[1:4,1:4],     # observation matrix
                                 probTrans = gamma[1:4,1:4],   # transition matrix
                                 len = K - first[i],           # nb of occasions
                                 checkRowSums = 0)             # do not check whether elements in a row sum tp 1
      }
})

#' 
#' Data in a list. 
## ------------------------------------------------------------------------
my.data <- list(y = y + 1, 
                alpha = c(1, 1, 1))

#' 
#' Constants in a list. 
## ------------------------------------------------------------------------
my.constants <- list(first = first, 
                     K = ncol(y), 
                     N = nrow(y))

#' 
#' Initial values. Note that we do not need initial values for the latent states anymore. 
## ------------------------------------------------------------------------
initial.values <- function(){list(phiA = runif(1, 0, 1), 
                                  phiB = runif(1, 0, 1), 
                                  phiC = runif(1, 0, 1), 
                                  psiA = rdirch(1, c(1,1,1)),
                                  psiB = rdirch(1, c(1,1,1)),
                                  psiC = rdirch(1, c(1,1,1)),
                                  pA = runif(1, 0, 1), 
                                  pB = runif(1, 0, 1), 
                                  pC = runif(1, 0, 1))}  

#' 
#' Parameters to monitor.
## ------------------------------------------------------------------------
parameters.to.save <- c("phiA", "phiB", "phiC", "psiA", "psiB", "psiC","pA", "pB", "pC")

#' 
#' MCMC settings.
## ------------------------------------------------------------------------
n.iter <- 10000
n.burnin <- 5000
n.chains <- 2

#' 
#' Run nimble.
## ----eval = FALSE--------------------------------------------------------
system.time(multisite.marginalized.out <- nimbleMCMC(code = multisite.marginalized, 
                                                     constants = my.constants,
                                                     data = my.data,              
                                                     inits = initial.values(),
                                                     monitors = parameters.to.save,
                                                     niter = n.iter,
                                                     nburnin = n.burnin, 
                                                     nchains = n.chains))

#' 
## ----echo = FALSE--------------------------------------------------------
load("geese.marginalized.RData")

#' 
#' The run took around 2 minutes. For comparison, the standard formulation (see live demo for class 6 on transition estimation) took around 4 minutes. 
#' 
#' Explore outputs.
## ------------------------------------------------------------------------
MCMCsummary(multisite.marginalized.out, round = 2)

#' 
#' ## Marginalization and weighted likelihood with nimbleEcology
#' 
#' In this section, we're gonna use nimble functions to express the likelihood using pooled encounter histories. We use a vector `mult` that contains the number of individuals with a particular encounter history. We hacked the `dHMM` nimbleEcology function below.  
## ------------------------------------------------------------------------
dHMMweighted <- nimbleFunction(
  run = function (x = double(1), 
                  init = double(1), 
                  probObs = double(2),
                  probTrans = double(2), 
                  len = double(0),
                  mult = double(0), # NEWLY ADDED: argument stating number of occurrences 
                                    # of same encounter history in entire dataset 
                  checkRowSums = integer(0, default = 0),  
                  log = integer(0, default = 0)) 
  {
    if (length(x) != len) 
      nimStop("In dHMM: Argument len must be length of x or 0.")
    if (nimDim(probObs)[1] != nimDim(probTrans)[1]) 
      nimStop("In dHMM: Length of dimension 1 in probObs must equal length of dimension 1 in probTrans.")
    if (nimDim(probTrans)[1] != nimDim(probTrans)[2]) 
      nimStop("In dHMM: probTrans must be a square matrix.")
    ## There was a strict test for sum(init) == 1.  This could be true in R and false in C++!
    if (abs(sum(init) - 1) > 1e-06) 
      nimStop("In dHMM: Initial probabilities must sum to 1.")
    if (checkRowSums) {
      transCheckPasses <- TRUE
      for (i in 1:nimDim(probTrans)[1]) {
        thisCheckSum <- sum(probTrans[i, ])
        if (abs(thisCheckSum - 1) > 1e-06) {
          nimPrint("In dHMM: Problem with sum(probTrans[i,]) with i = ", 
                   i, ". The sum should be 1 but is ", thisCheckSum)
          transCheckPasses <- FALSE
        }
      }
      obsCheckPasses <- TRUE
      for (i in 1:nimDim(probObs)[1]) {
        thisCheckSum <- sum(probObs[i, ])
        if (abs(thisCheckSum - 1) > 1e-06) {
          nimPrint("In dHMM: Problem with sum(probObs[i,]) with i = ", 
                   i, ". The sum should be 1 but is ", thisCheckSum)
          obsCheckPasses <- FALSE
        }
      }
      if (!(transCheckPasses | obsCheckPasses)) 
        nimStop("In dHMM: probTrans and probObs were not specified correctly.  Probabilities in each row (second dimension) must sum to 1.")
      if (!transCheckPasses) 
        nimStop("In dHMM: probTrans was not specified correctly.  Probabilities in each row (second dimension) must sum to 1.")
      if (!obsCheckPasses) 
        nimStop("In dHMM: probObs was not specified correctly. Probabilities in each row must sum to 1.")
    }
    pi <- init
    logL <- 0
    nObsClasses <- nimDim(probObs)[2]
    for (t in 1:len) {
      if (x[t] > nObsClasses | x[t] < 1) 
        nimStop("In dHMM: Invalid value of x[t].")
      Zpi <- probObs[, x[t]] * pi
      sumZpi <- sum(Zpi)
      logL <- logL + log(sumZpi) * mult # NEWLY ADDED
      if (t != len) 
        pi <- ((Zpi %*% probTrans)/sumZpi)[1, ]
    }
    if (log) 
      return(logL)
    return(exp(logL))
    returnType(double())
  })

#' 
#' Write the model.
## ------------------------------------------------------------------------
multisite.marginalized <- nimbleCode({
  
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
  
  # survival priors
  phiA ~ dunif(0, 1)
  phiB ~ dunif(0, 1)
  phiC ~ dunif(0, 1)
  # detection priors
  pA ~ dunif(0, 1)
  pB ~ dunif(0, 1)
  pC ~ dunif(0, 1)
  # transition priors: Dirichlet
  psiA[1:3] ~ ddirch(alpha[1:3])
  psiB[1:3] ~ ddirch(alpha[1:3])
  psiC[1:3] ~ ddirch(alpha[1:3])
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
  
  for(i in 1:N) {
    init[i, 1:4] <- gamma[ y[i, first[i] ] - 1, 1:4 ] # First state propagation
    }
  
  # likelihood 
  for (i in 1:N){
    y[i,(first[i]+1):K] ~ dHMMweighted(init = init[i,1:4], # count data from first[i] + 1
                                       mult = mult[i],
                                       probObs = omega[1:4,1:4],
                                       probTrans = gamma[1:4,1:4],
                                       len = K - first[i],
                                       checkRowSums = 0)
  }
})

#' 
#' We need to pool the individual encounter histories by unique encounter histories, and to record the number of individuals with a particular unique encounter history. 
## ------------------------------------------------------------------------
geese <- read_csv("geese.csv", col_names = TRUE)
y <- as.matrix(geese)
y_weighted <- y %>% 
  as_tibble() %>% 
  group_by_all() %>% 
  summarise(mult = n()) %>% 
  relocate(mult) %>% 
  as.matrix()
head(y_weighted)
mult <- y_weighted[,1] # nb of individuals w/ a particular encounter history
y <- y_weighted[,-1] # pooled data

#' 
#' There are 5 individuals that were detected only once, in site A in year 1989, 8 individuals that were detected only once, in 1989 in site B, and so on. 
#' 
#' Get the occasion of first capture for each history.
## ------------------------------------------------------------------------
get.first <- function(x) min(which(x != 0))
first <- apply(y, 1, get.first)

#' 
#' Filter out individuals that are first captured at last occasion.
## ------------------------------------------------------------------------
mask <- which(first!=ncol(y))
y <- y[mask, ]

#' 
#' Apply filter on occasion of first capture and sample size.
## ------------------------------------------------------------------------
first <- first[mask]
mult <- mult[mask]

#' 
#' Data in a list. 
## ------------------------------------------------------------------------
my.data <- list(y = y + 1, 
                alpha = c(1, 1, 1))

#' 
#' Constants in a list. 
## ------------------------------------------------------------------------
my.constants <- list(first = first, 
                     K = ncol(y), 
                     N = nrow(y),
                     mult = mult)

#' 
#' Initial values. Note that we do not need initial values for the latent states anymore. 
## ------------------------------------------------------------------------
initial.values <- function(){list(phiA = runif(1, 0, 1), 
                                  phiB = runif(1, 0, 1), 
                                  phiC = runif(1, 0, 1), 
                                  psiA = rdirch(1, c(1,1,1)),
                                  psiB = rdirch(1, c(1,1,1)),
                                  psiC = rdirch(1, c(1,1,1)),
                                  pA = runif(1, 0, 1), 
                                  pB = runif(1, 0, 1), 
                                  pC = runif(1, 0, 1))}  

#' 
#' Parameters to monitor.
## ------------------------------------------------------------------------
parameters.to.save <- c("phiA", "phiB", "phiC", "psiA", "psiB", "psiC","pA", "pB", "pC")

#' 
#' MCMC settings.
## ------------------------------------------------------------------------
n.iter <- 10000
n.burnin <- 5000
n.chains <- 2

#' 
#' Run nimble. 
## ----eval = FALSE--------------------------------------------------------
system.time(multisite.marginalized.out <- nimbleMCMC(code = multisite.marginalized, 
                                                     constants = my.constants,
                                                     data = my.data,              
                                                     inits = initial.values(),
                                                     monitors = parameters.to.save,
                                                     niter = n.iter,
                                                     nburnin = n.burnin, 
                                                     nchains = n.chains))

#' 
#' This run won't work cause we need a `rHMMweighted` as well. 
## ----eval = FALSE--------------------------------------------------------
rHMMweighted <- nimbleFunction(
  run = function(n = integer(),    ## Observed capture (state) history
                 init = double(1),
                 probObs = double(2),
                 probTrans = double(2),
                 len = double(0, default = 0),
                 mult = double(0),
                 checkRowSums = double(0, default = 1)) {
    returnType(double(1))
    if (dim(probObs)[1] != dim(probTrans)[1]) stop("In rHMM: Number of cols in probObs must equal number of cols in probTrans.")
    if (dim(probTrans)[1] != dim(probTrans)[2]) stop("In rHMM: probTrans must be a square matrix.")
    if (abs(sum(init) - 1) > 1e-06) stop("In rHMM: Initial probabilities must sum to 1.")
    if (checkRowSums) {
      transCheckPasses <- TRUE
      for (i in 1:dim(probTrans)[1]) {
        thisCheckSum <- sum(probTrans[i,])
        if (abs(thisCheckSum - 1) > 1e-6) {
          ## Compilation doesn't support more than a simple string for stop()
          ## so we provide more detail using a print().
          print("In rHMM: Problem with sum(probTrans[i,]) with i = ", i, ". The sum should be 1 but is ", thisCheckSum)
          transCheckPasses <- FALSE
        }
      }
      obsCheckPasses <- TRUE
      for (i in 1:dim(probObs)[1]) {
        thisCheckSum <- sum(probObs[i,])
        if (abs(thisCheckSum - 1) > 1e-6) {
          print("In rHMM: Problem with sum(probObs[i,]) with i = ", i, ". The sum should be 1 but is ", thisCheckSum)
          obsCheckPasses <- FALSE
        }
      }
      if(!(transCheckPasses | obsCheckPasses))
        stop("In rHMM: probTrans and probObs were not specified correctly.  Probabilities in each row (second dimension) must sum to 1.")
      if(!transCheckPasses)
        stop("In rHMM: probTrans was not specified correctly.  Probabilities in each row (second dimension) must sum to 1.")
      if(!obsCheckPasses)
        stop("In rHMM: probObs was not specified correctly. Probabilities in each row must sum to 1.")
    }
    ans <- numeric(len)
    probInit <- init
    trueInit <- 0
    r <- runif(1, 0, 1)
    j <- 1
    while (r > sum(probInit[1:j])) j <- j + 1
    trueInit <- j
    trueState <- trueInit
    for (i in 1:len) {
      # Transition to a new true state
      r <- runif(1, 0, 1)
      j <- 1
      while (r > sum(probTrans[trueState, 1:j])) j <- j + 1
      trueState <- j
      # Detect based on the true state
      r <- runif(1, 0, 1)
      j <- 1
      while (r > sum(probObs[trueState, 1:j])) j <- j + 1
      ans[i] <- j
    }
    return(ans)
  })

#' 
#' Run nimble again. 
## ----eval = FALSE--------------------------------------------------------
system.time(multisite.marginalized.out <- nimbleMCMC(code = multisite.marginalized, 
                                                     constants = my.constants,
                                                     data = my.data,              
                                                     inits = initial.values(),
                                                     monitors = parameters.to.save,
                                                     niter = n.iter,
                                                     nburnin = n.burnin, 
                                                     nchains = n.chains))

#' 
#' Wow, this run took 1.5 minute, even faster than the marginalized version of the likelihood. The outputs are similar. 
## ----echo = FALSE--------------------------------------------------------
load(here::here("worksheets", "weighted_geese.RData"))

#' 
## ------------------------------------------------------------------------
MCMCsummary(multisite.marginalized.out, round = 2)

#' 
#' ## Analysis of the whole Canada geese dataset
#' 
#' So far we've worked with a subset of the original dataset. Fitting the same model to the whole data would be difficult, if not impossible. Let's try it with the weighted likelihood. We first read in the data. 
## ------------------------------------------------------------------------
geese <- read_csv2("allgeese.csv", col_names = TRUE)
geese <- as.matrix(geese)
y <- geese[,-7]
mult <- geese[,7]

#' 
#' In total, we have `r sum(mult)` banded geese, and only `r length(mult)` unique capture histories.
#' 
#' Get the occasion of first capture for each individual
## ------------------------------------------------------------------------
get.first <- function(x) min(which(x != 0))
first <- apply(y, 1, get.first)

#' 
#' Filter out individuals that are first captured at last occasion
## ------------------------------------------------------------------------
mask <- which(first!=ncol(y))
y <- y[mask, ]

#' 
#' Recalculate occasion of first capture and sample size
## ------------------------------------------------------------------------
first <- first[mask]
mult <- mult[mask]

#' 
#' Data in list. 
## ------------------------------------------------------------------------
my.data <- list(y = y + 1, 
                alpha = c(1, 1, 1))

#' 
#' Constant in a list.
## ------------------------------------------------------------------------
my.constants <- list(first = first, 
                     K = ncol(y), 
                     N = nrow(y),
                     mult = mult)

#' 
#' Initial values.
## ------------------------------------------------------------------------
initial.values <- function(){list(phiA = runif(1, 0, 1), 
                                  phiB = runif(1, 0, 1), 
                                  phiC = runif(1, 0, 1), 
                                  psiA = rdirch(1, c(1,1,1)),
                                  psiB = rdirch(1, c(1,1,1)),
                                  psiC = rdirch(1, c(1,1,1)),
                                  pA = runif(1, 0, 1), 
                                  pB = runif(1, 0, 1), 
                                  pC = runif(1, 0, 1))}

#' 
#' MCMC settings.
## ------------------------------------------------------------------------
parameters.to.save <- c("phiA", "phiB", "phiC", "psiA", "psiB", "psiC","pA", "pB", "pC")
n.iter <- 10000
n.burnin <- 5000
n.chains <- 2

#' 
#' Run nimble.
## ----eval = FALSE--------------------------------------------------------
system.time(multisite.marginalized.out <- nimbleMCMC(code = multisite.marginalized, 
                                                     constants = my.constants,
                                                     data = my.data,              
                                                     inits = initial.values(),
                                                     monitors = parameters.to.save,
                                                     niter = n.iter,
                                                     nburnin = n.burnin, 
                                                     nchains = n.chains))

#' 
#' It would take ages to run the same model on the > 20,000 individual encounter histories. Not even sure it can be run.
#' 
#' Let's inspect the results. 
#' 
## ----echo = FALSE--------------------------------------------------------
load(here::here("worksheets", "allgeese_weighted.RData"))

#' 
## ------------------------------------------------------------------------
MCMCsummary(multisite.marginalized.out, round = 2)

#' 
#' <!-- knitr::purl(here::here("worksheets","7_demo.Rmd"), documentation = 2) -->
