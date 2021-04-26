#' ---
#' title: "Exercises"
#' author: "The team"
#' date: "last updated: `r Sys.Date()`"
#' output: html_document
#' ---
#' 
#' 
## ----setup, include=FALSE, echo=FALSE, message=FALSE, warning=FALSE----
options(htmltools.dir.version = FALSE)
knitr::opts_chunk$set(comment = "")
library(tidyverse)
theme_set(theme_light(base_size = 14))
update_geom_defaults("point", list(size = 2)) 
library(here)
library(nimble)

#' 
#' ## To train yourself
#' 
#' Two questions for you to go further with what we've seen on the first day. We will use the dipper dataset again. 
#' 
#' + Build and fit in nimble a model with an additive effect of sex and wing length on survival. 
#' 
#' + Build and fit in nimble a model with an age effect on survival. For age, consider a two-age covariate, with a survival parameter for the time interval right after first capture, and a parameter for beyond. 
#' 
#' ## Preliminary steps
#' 
#' Load packages. 
## ---------------------------------------------------------------
library(nimble)
library(tidyverse)
library(MCMCvis)

#' 
#' Read in the data. 
## ---------------------------------------------------------------
dipper <- read_csv("dipper.csv")

#' 
#' Format the data.
## ---------------------------------------------------------------
y <- dipper %>%
  select(year_1981:year_1987) %>%
  as.matrix()
head(y)

#' 
#' 
#' ## Effect of sex and wing length
#' 
#' Write the model. We use nested indexing with the sex index that contains 1 if the bird is a male, and 2 otherwise. We have $\logit(\phi_i) = \beta_1 + \beta_3 * \text{winglength}_i$ for males, and $\logit(\phi_i) = \beta_2 + \beta_3 * \text{winglength}_i$ for females. 
## ---------------------------------------------------------------
hmm.phisexwlp <- nimbleCode({
  p ~ dunif(0, 1) # prior detection
  omega[1,1] <- 1 - p    # Pr(alive t -> non-detected t)
  omega[1,2] <- p        # Pr(alive t -> detected t)
  omega[2,1] <- 1        # Pr(dead t -> non-detected t)
  omega[2,2] <- 0        # Pr(dead t -> detected t)
  for (i in 1:N){
    logit(phi[i]) <- beta[sex[i]] + beta[3] * winglength[i]
    gamma[1,1,i] <- phi[i]      # Pr(alive t -> alive t+1)
    gamma[1,2,i] <- 1 - phi[i]  # Pr(alive t -> dead t+1)
    gamma[2,1,i] <- 0           # Pr(dead t -> alive t+1)
    gamma[2,2,i] <- 1           # Pr(dead t -> dead t+1)
  }
  beta[1] ~ dnorm(mean = 0, sd = 1.5) # intercept male
  beta[2] ~ dnorm(mean = 0, sd = 1.5) # intercept female
  beta[3] ~ dnorm(mean = 0, sd = 1.5) # slope wing length
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
#' Constants in a list. Note we standardise wing length. 
## ---------------------------------------------------------------
first <- apply(y, 1, function(x) min(which(x !=0)))
wing.length.st <- as.vector(scale(dipper$wing_length))
my.constants <- list(N = nrow(y), 
                     T = ncol(y), 
                     first = first,
                     winglength = wing.length.st,
                     sex = if_else(dipper$sex == "M", 1, 2))

#' 
#' Data in a list. 
## ---------------------------------------------------------------
my.data <- list(y = y + 1)

#' 
#' Initial values. 
## ---------------------------------------------------------------
zinits <- y
zinits[zinits == 0] <- 1
initial.values <- function() list(beta = rnorm(3,0,5),
                                  p = runif(1,0,1),
                                  z = zinits)

#' 
#' Parameters to be monitored. 
## ---------------------------------------------------------------
parameters.to.save <- c("beta", "p")

#' 
#' MCMC details.
## ---------------------------------------------------------------
n.iter <- 5000
n.burnin <- 2500
n.chains <- 2

#' 
#' Run nimble.
## ---------------------------------------------------------------
mcmc.phisexwlp <- nimbleMCMC(code = hmm.phisexwlp, 
                             constants = my.constants,
                             data = my.data,              
                             inits = initial.values,
                             monitors = parameters.to.save,
                             niter = n.iter,
                             nburnin = n.burnin, 
                             nchains = n.chains)

#' 
#' Display results. 
## ---------------------------------------------------------------
MCMCsummary(mcmc.phisexwlp, round = 2)

#' 
#' Let's visualise survival as a function of wing length for both sexes. 
#' 
#' First we put together the values from the two chains we generated in the posterior distributions of the regression parameters. 
## ---------------------------------------------------------------
beta1 <- c(mcmc.phisexwlp$chain1[,'beta[1]'], mcmc.phisexwlp$chain2[,'beta[1]'])
beta2 <- c(mcmc.phisexwlp$chain1[,'beta[2]'], mcmc.phisexwlp$chain2[,'beta[2]'])
beta3 <- c(mcmc.phisexwlp$chain1[,'beta[3]'], mcmc.phisexwlp$chain2[,'beta[3]'])

#' 
#' We get survival estimates for each MCMC iteration. 
## ---------------------------------------------------------------
predicted_survivalM <- matrix(NA, nrow = length(beta1), ncol = length(my.constants$winglength))
predicted_survivalF <- matrix(NA, nrow = length(beta1), ncol = length(my.constants$winglength))
for (i in 1:length(beta1)){
  for (j in 1:length(my.constants$winglength)){
    predicted_survivalM[i,j] <- plogis(beta1[i] + beta3[i] * my.constants$winglength[j]) # males
    predicted_survivalF[i,j] <- plogis(beta2[i] + beta3[i] * my.constants$winglength[j]) # females
  }
}

#' 
#' From here, we may calculate posterior mean and credible intervals. 
## ---------------------------------------------------------------
mean_survivalM <- apply(predicted_survivalM, 2, mean)
lciM <- apply(predicted_survivalM, 2, quantile, prob = 2.5/100)
uciM <- apply(predicted_survivalM, 2, quantile, prob = 97.5/100)
mean_survivalF <- apply(predicted_survivalF, 2, mean)
lciF <- apply(predicted_survivalF, 2, quantile, prob = 2.5/100)
uciF <- apply(predicted_survivalF, 2, quantile, prob = 97.5/100)
ord <- order(my.constants$winglength)
df <- data.frame(wing_length = c(my.constants$winglength[ord], my.constants$winglength[ord]),
                 survival = c(mean_survivalM[ord], mean_survivalF[ord]),
                 lci = c(lciM[ord],lciF[ord]),
                 uci = c(uciM[ord],uciF[ord]),
                 sex = c(rep("male", length(mean_survivalM)), rep("female", length(mean_survivalF))))

#' 
#' Now on a plot. 
## ---- fig.width = 7.5, fig.asp = 0.618, dev = "svg"-------------
df %>%
  ggplot() + 
  aes(x = wing_length, y = survival, color = sex) + 
  geom_line() + 
  geom_ribbon(aes(ymin = lci, ymax = uci, fill = sex), alpha = 0.5) + 
  ylim(0,1) + 
  labs(x = "wing length", y = "estimated survival", color = "", fill = "")

#' 
#' 
#' ## Incorporating age
#' 
#' Age in capture-recapture has a particular meaning in capture-recapture analyses. It is the time elapsed since first encounter, which is a proxy of true age obviously, but not true age. Of course, if age is known at first encounter, then it is the true age. 
#' 
#' Another important remark is that age is an individual covariate, but in contrast with the wing length covariate we considered in the previous examples, age varies over time. The cool thing is that it has no missing value as age at $t+1$ is just age at $t$ to which we add 1. This suggests a way to code the age effect in nimble as follows. 
#' 
## ---------------------------------------------------------------
hmm.phiage.in <- nimbleCode({
  p ~ dunif(0, 1) # prior detection
  omega[1,1] <- 1 - p    # Pr(alive t -> non-detected t)
  omega[1,2] <- p        # Pr(alive t -> detected t)
  omega[2,1] <- 1        # Pr(dead t -> non-detected t)
  omega[2,2] <- 0        # Pr(dead t -> detected t)
  for (i in 1:N){
    for (t in first[i]:(T-1)){
    logit(phi[i,t]) <- beta[1] + beta[2] * equals(t, first[i]) # phi1 = beta1 + beta2 and phi1+ = beta1
    gamma[1,1,i,t] <- phi[i,t]      # Pr(alive t -> alive t+1)
    gamma[1,2,i,t] <- 1 - phi[i,t]  # Pr(alive t -> dead t+1)
    gamma[2,1,i,t] <- 0           # Pr(dead t -> alive t+1)
    gamma[2,2,i,t] <- 1           # Pr(dead t -> dead t+1)
    }
  }
  beta[1] ~ dnorm(mean = 0, sd = 1.5) # phi1+
  beta[2] ~ dnorm(mean = 0, sd = 1.5) # phi1 - phi1+
  phi1plus <- plogis(beta[1])         # phi1+
  phi1 <- plogis(beta[1] + beta[2])   # phi1
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, i, j-1])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})

#' 
#' Constants in a list.
## ---------------------------------------------------------------
first <- apply(y, 1, function(x) min(which(x !=0)))
my.constants <- list(N = nrow(y), 
                     T = ncol(y), 
                     first = first)

#' 
#' Data in a list. 
## ---------------------------------------------------------------
my.data <- list(y = y + 1)

#' 
#' Initial values. 
## ---------------------------------------------------------------
zinits <- y
zinits[zinits == 0] <- 1
initial.values <- function() list(beta = rnorm(2,0,5),
                                  p = runif(1,0,1),
                                  z = zinits)

#' 
#' Parameters to be monitored. 
## ---------------------------------------------------------------
parameters.to.save <- c("phi1", "phi1plus", "p")

#' 
#' MCMC details.
## ---------------------------------------------------------------
n.iter <- 5000
n.burnin <- 2500
n.chains <- 2

#' 
#' Run nimble.
## ---------------------------------------------------------------
mcmc.phi.age.in <- nimbleMCMC(code = hmm.phiage.in, 
                             constants = my.constants,
                             data = my.data,              
                             inits = initial.values,
                             monitors = parameters.to.save,
                             niter = n.iter,
                             nburnin = n.burnin, 
                             nchains = n.chains)

#' 
#' Display results. 
## ---------------------------------------------------------------
MCMCsummary(mcmc.phi.age.in, round = 2)

#' 
#' 
#' Another method to include an age effect is to create an individual by time covariate and use nested indexing (as in the previous example) to distinguish survival after first detection from survival afterwards. 
## ---------------------------------------------------------------
age <- matrix(NA, nrow = nrow(y), ncol = ncol(y) - 1)
for (i in 1:nrow(age)){
  for (j in 1:ncol(age)){
    if (j == first[i]) age[i,j] <- 1
    if (j > first[i]) age[i,j] <- 2
  }
}

#' 
#' The model. Careful, now survival is no longer defined on the logit scale as in the previous model, so we use uniform priors.  
## ---------------------------------------------------------------
hmm.phiage.out <- nimbleCode({
  p ~ dunif(0, 1) # prior detection
  omega[1,1] <- 1 - p    # Pr(alive t -> non-detected t)
  omega[1,2] <- p        # Pr(alive t -> detected t)
  omega[2,1] <- 1        # Pr(dead t -> non-detected t)
  omega[2,2] <- 0        # Pr(dead t -> detected t)
  for (i in 1:N){
    for (t in first[i]:(T-1)){
    phi[i,t] <- beta[age[i,t]] # beta1 = phi1, beta2 = phi1+
    gamma[1,1,i,t] <- phi[i,t]      # Pr(alive t -> alive t+1)
    gamma[1,2,i,t] <- 1 - phi[i,t]  # Pr(alive t -> dead t+1)
    gamma[2,1,i,t] <- 0           # Pr(dead t -> alive t+1)
    gamma[2,2,i,t] <- 1           # Pr(dead t -> dead t+1)
    }
  }
  beta[1] ~ dunif(0, 1) # phi1
  beta[2] ~ dunif(0, 1) # phi1+
  phi1 <- beta[1]
  phi1plus <- beta[2]
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, i, j-1])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})

#' 
#' Constants in a list, inculding the age matrix covariate.
## ---------------------------------------------------------------
first <- apply(y, 1, function(x) min(which(x !=0)))
my.constants <- list(N = nrow(y), 
                     T = ncol(y), 
                     first = first,
                     age = age)

#' 
#' Data in a list. 
## ---------------------------------------------------------------
my.data <- list(y = y + 1)

#' 
#' Initial values. 
## ---------------------------------------------------------------
zinits <- y
zinits[zinits == 0] <- 1
initial.values <- function() list(beta = runif(2,0,1),
                                  p = runif(1,0,1),
                                  z = zinits)

#' 
#' Parameters to be monitored. 
## ---------------------------------------------------------------
parameters.to.save <- c("phi1", "phi1plus", "p")

#' 
#' MCMC details.
## ---------------------------------------------------------------
n.iter <- 5000
n.burnin <- 2500
n.chains <- 2

#' 
#' Run nimble.
## ---------------------------------------------------------------
mcmc.phi.age.out <- nimbleMCMC(code = hmm.phiage.out, 
                               constants = my.constants,
                               data = my.data,              
                               inits = initial.values,
                               monitors = parameters.to.save,
                               niter = n.iter,
                               nburnin = n.burnin, 
                               nchains = n.chains)

#' 
#' Display results. 
## ---------------------------------------------------------------
MCMCsummary(mcmc.phi.age.out, round = 2)

#' 
#' 
#' 
#' <!-- knitr::purl(here::here("worksheets","exercises.Rmd"), documentation = 2) -->
