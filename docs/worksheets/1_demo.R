#' ---
#' title: "Class 2 live demo: Crash course on Bayesian statistics and MCMC algorithms"
#' author: "The team"
#' date: "last updated: `r Sys.Date()`"
#' output: html_document
#' ---
#' 
#' 
## ----setup, include=FALSE, echo=FALSE-------------------------------------------------------------------------
options(htmltools.dir.version = FALSE)
knitr::opts_chunk$set(comment = "", message = FALSE, warning = FALSE)
library(tidyverse)
theme_set(theme_light(base_size = 14))
update_geom_defaults("point", list(size = 2)) 
library(here)
library(nimble)

#' 
#' ## Introduction
#' 
#' In this demo, we show how to implement the Bayes theorem by brute force with numerical integration, then we use a Metropolis algorithm. We use the survival example with $y = 19$ alive individuals out of $n = 57$ released. A binomial likelihood is assumed, with a beta prior on survival. 
#' 
#' Load the packages we will need. 
## -------------------------------------------------------------------------------------------------------------
library(tidyverse)

#' 
#' The data.
## -------------------------------------------------------------------------------------------------------------
y <- 19 # nb of success
n <- 57 # nb of attempts

#' 
#' ## Brute force via numerical integration
#' 
#' First we defined the prior. 
## ---- echo = TRUE, fig.width = 7.5, fig.asp = 0.618, dev = "svg"----------------------------------------------
a <- 1 # param of the prior
b <- 1 # param of the prior
p <- seq(0,1,.002) # grid
prior <- dbeta(p,a,b) # eval beta prior
dfprior <- data.frame(p = p, prior = prior) 
dfprior %>%
  ggplot() + 
  geom_line(aes(x = p, y = prior), 
            size = 1.5,
            color = wesanderson::wes_palettes$Royal1[1])

#' 
#' Now we form the numerator and denominator of the Bayes formula. 
## ---- echo = TRUE, fig.width = 7.5, fig.asp = 0.618, dev = "svg"----------------------------------------------
numerator <- function(p) dbinom(y,n,p) * dbeta(p,a,b) # likelihood x prior
denominator <- integrate(numerator,0,1)$value # denominator
numerical_posterior <- data.frame(p = p, posterior = numerator(p)/denominator) 
numerical_posterior %>%
  ggplot() + 
  geom_line(aes(x = p, y = posterior), 
            size = 1.5, 
            col = wesanderson::wes_palettes$Royal1[2], 
            alpha = 0.5)

#' 
#' Now we compare the posterior distribution we obtained by numerical integration to the explicit distribution we get through conjugacy (binomial likelihood and beta prior gives beta posterior). 
## ---- echo = TRUE, fig.width = 7.5, fig.asp = 0.618, dev = "svg"----------------------------------------------
explicit_posterior <- dbeta(p, y + a, n - y + b) # beta posterior
dfexpposterior <- data.frame(p = p, explicit_posterior = explicit_posterior)
ggplot() + 
  geom_line(data = numerical_posterior, 
            aes(x = p, y = posterior), 
            size = 1.5, 
            col = wesanderson::wes_palettes$Royal1[2],
            alpha = 0.5) + 
  geom_line(data = dfexpposterior, 
            aes(x = p, y = explicit_posterior),
            size = 1.5, 
            col = wesanderson::wes_palettes$Royal1[3], 
            linetype = "dashed")

#' 
#' Last, we put everything on the same plot. 
## ---- echo = TRUE, fig.width = 7.5, fig.asp = 0.618, dev = "svg"----------------------------------------------
ggplot() + 
  geom_line(data = numerical_posterior, 
            aes(x = p, y = posterior), 
            size = 1.5, 
            col = wesanderson::wes_palettes$Royal1[2], 
            alpha = 0.5) + 
  geom_line(data = dfexpposterior, 
            aes(x = p, y = explicit_posterior),
            col = wesanderson::wes_palettes$Royal1[3], 
            size = 1.5, 
            linetype = "dashed") + 
  geom_line(data = dfprior,
            aes(x = p, y = prior),
            col = wesanderson::wes_palettes$Royal1[1],
            size = 1.5)

#' 
#' 
#' ## Implement the Metropolis algorithm
#' 
#' In this section, we implement the Metropolis algorithm. 
#' 
#' First the data. 
## ---- echo = TRUE, fig.width = 7.5, fig.asp = 0.618, dev = "svg"----------------------------------------------
# survival data, 19 "success" out of 57 "attempts"
survived <- 19
released <- 57

#' 
#' Then we define the log-likelihood function, which is binomial. Note that we use the log. 
## ---- echo = TRUE, fig.width = 7.5, fig.asp = 0.618, dev = "svg"----------------------------------------------
# log-likelihood function
loglikelihood <- function(x, p){
  dbinom(x = x, size = released, prob = p, log = TRUE)
}

#' 
#' We then define the prior, a uniform distribution between 0 and 1, or a beta distribution with parameters 1 and 1.  
## ---- echo = TRUE, fig.width = 7.5, fig.asp = 0.618, dev = "svg"----------------------------------------------
# prior density
logprior <- function(p){
  dunif(x = p, min = 0, max = 1, log = TRUE)
}

#' 
#' Eventually, we form the posterior density distribution, on the log scale. Note that we commented out the denominator because we won't need it in the Metropolis algorithm (it cancels out when we form the acceptance ratio. 
## ---- echo = TRUE, fig.width = 7.5, fig.asp = 0.618, dev = "svg"----------------------------------------------
# posterior density function (log scale)
posterior <- function(x, p){
  loglikelihood(x, p) + logprior(p) # - log(Pr(data))
}

#' 
#' Now we need to propose a candidate value for survival at each iteration of the algorithm. To do so, we work on the logit scale, and use a normal distribution to jump away from the current value. 
## ---- echo = TRUE, fig.width = 7.5, fig.asp = 0.618, dev = "svg"----------------------------------------------
# propose candidate value
move <- function(x, away = .2){ 
  logitx <- log(x / (1 - x))
  logit_candidate <- logitx + rnorm(1, 0, away)
  candidate <- plogis(logit_candidate)
  return(candidate)
}

#' 
#' We pick 100 iterations, pre-allocate some memory to store the results and set the seed for reproducibility.
## ---- echo = TRUE, fig.width = 7.5, fig.asp = 0.618, dev = "svg"----------------------------------------------
# set up the scene
steps <- 100
theta.post <- rep(NA, steps)
set.seed(1234)

#' 
#' To start the algorithm, we also need an initial value for survival. It's going to be 0.5 here. This is our first value stored. When running several chains, pick different initial values. 
## ---- echo = TRUE, fig.width = 7.5, fig.asp = 0.618, dev = "svg"----------------------------------------------
# pick starting value (step 1)
inits <- 0.5
theta.post[1] <- inits

#' 
#' Now we run the Metropolis algorithm. 
## ---- echo = TRUE, fig.width = 7.5, fig.asp = 0.618, dev = "svg"----------------------------------------------
for (t in 2:steps){ # repeat steps 2-4 (step 5)
  
  # propose candidate value for prob of success (step 2)
  theta_star <- move(theta.post[t-1])
  
  # calculate ratio R (step 3)
  pstar <- posterior(survived, p = theta_star) 
  pprev <- posterior(survived, p = theta.post[t-1])
  logR <- pstar - pprev # here we see that the denominator log(Pr(data)) cancels out
                        # logR = loglikelihood(survived, theta_star) + logprior(theta_star) - log(Pr(data)) -          
                        #        (loglikelihood(survived, theta.post[t-1]) + logprior(theta.post[t-1]) - log(Pr(data)))
  R <- exp(logR)
  
  # decide to accept candidate value or to keep current value (step 4)
  accept <- rbinom(1, 1, prob = min(R, 1))
  theta.post[t] <- ifelse(accept == 1, theta_star, theta.post[t-1])
}

#' 
#' Let's have a look to the values we generated from the posterior distribution. 
## ---- echo = TRUE, fig.width = 7.5, fig.asp = 0.618, dev = "svg"----------------------------------------------
head(theta.post)
tail(theta.post)

#' 
#' And now, the trace of the chain.
## ---- echo = TRUE, fig.width = 7.5, fig.asp = 0.618, dev = "svg"----------------------------------------------
df <- data.frame(x = 1:steps, y = theta.post)
df %>%
  ggplot() +
  geom_line(aes(x = x, y = y), size = 1.5, color = wesanderson::wes_palettes$Zissou1[1]) + 
  labs(x = "iterations", y = "values from posterior distribution") + 
  ylim(0.1, 0.6)

#' 
#' Let's run another chain, now with initial value 0.2 for survival, and visualise the trace of the two chains.
## ---- echo = TRUE, fig.width = 7.5, fig.asp = 0.618, dev = "svg"----------------------------------------------
# pick starting value (step 1)
inits <- 0.2
theta.post2 <- rep(NA, steps)
theta.post2[1] <- inits

for (t in 2:steps){ # repeat steps 2-4 (step 5)
  # propose candidate value for prob of success (step 2)
  theta_star <- move(theta.post2[t-1])
  # calculate ratio R (step 3)
  pstar <- posterior(survived, p = theta_star)  
  pprev <- posterior(survived, p = theta.post[t-1])
  logR <- pstar - pprev
  R <- exp(logR)
  
  # decide to accept candidate value or to keep current value (step 4)
  accept <- rbinom(1, 1, prob = min(R, 1))
  theta.post2[t] <- ifelse(accept == 1, theta_star, theta.post2[t-1])
}

df2 <- data.frame(x = 1:steps, y = theta.post2)
  ggplot() +
  geom_line(data = df, aes(x = x, y = y), size = 1.5, color = wesanderson::wes_palettes$Zissou1[1]) + 
  geom_line(data = df2, aes(x = x, y = y), size = 1.5, color = wesanderson::wes_palettes$Zissou1[3]) + 
  labs(x = "iterations", y = "values from posterior distribution") + 
  ylim(0.1, 0.6)

#' 
#' 
#' Last, we re-run these two chains with more iterations, say 5000. We also add the posterior mean in yellow, and the maximum likelihood estimate in red. 
## ---- echo = TRUE, fig.width = 7.5, fig.asp = 0.618, dev = "svg"----------------------------------------------
# set up the scene
steps <- 5000
theta.post <- rep(NA, steps)
set.seed(1234)

# pick starting value (step 1)
inits <- 0.5
theta.post[1] <- inits

for (t in 2:steps){ # repeat steps 2-4 (step 5)
  
  # propose candidate value for prob of success (step 2)
  theta_star <- move(theta.post[t-1])
  
  # calculate ratio R (step 3)
  pstar <- posterior(survived, p = theta_star)  
  pprev <- posterior(survived, p = theta.post[t-1])
  logR <- pstar - pprev
  R <- exp(logR)
  
  # decide to accept candidate value or to keep current value (step 4)
  accept <- rbinom(1, 1, prob = min(R, 1))
  theta.post[t] <- ifelse(accept == 1, theta_star, theta.post[t-1])
}

df <- data.frame(x = 1:steps, y = theta.post)
df %>%
  ggplot() +
  geom_line(aes(x = x, y = y), size = 1, color = wesanderson::wes_palettes$Zissou1[1]) + 
  labs(x = "iterations", y = "values from posterior distribution") + 
  ylim(0.1, 0.6) + 
  geom_hline(aes(yintercept = mean(theta.post)), 
             color = wesanderson::wes_palettes$Zissou1[3],
             size = 1.2) + 
  geom_hline(aes(yintercept = 19/57), 
             color = wesanderson::wes_palettes$Zissou1[5],
             size = 1.2)

