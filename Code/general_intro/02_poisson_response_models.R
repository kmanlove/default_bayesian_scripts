# Basic GLM examples in JAGS
# K. Manlove
# 2024-01-04

# I. Pull in data from glmmTMB
library(glmmTMB)
head(Owls)

# II. Bayesian implementations of the nestling models ----
## A. Naive Poisson model ----
nestlings_jags_data <- list(SiblingNegotiation = Owls$SiblingNegotiation, 
                            FoodTreatment = as.numeric(as.factor(Owls$FoodTreatment)) - 1,
                            # note that we have to convert the FoodTreatment factor to
                            # a number to work appropriately in JAGS. Also, numeric coding needs to start at 0, not 1
                            SexParent = as.numeric(as.factor(Owls$SexParent)) - 1, 
                            # same with SexParent. Also, numeric coding needs to start at 0
                            N_owls = nrow(Owls))

owls_naive <- function(){
  # Priors:
  beta_0 ~ dnorm(0, 0.001)
  beta_ft ~ dnorm(0, 0.001)
  beta_ps ~ dnorm(0, 0.001)
  beta_inter ~ dnorm(0, 0.001)
  
  # Likelihood:
  for(i in 1:N_owls){
    # linear combination of predictors
    lin_comb[i] <- beta_0 + beta_ft * FoodTreatment[i] + beta_ps * SexParent[i] + 
      beta_inter * FoodTreatment[i] * SexParent[i]
    # invert link to get focal parameter
    lambda[i] <- exp(lin_comb[i])
    # specify likelihood prob model
    SiblingNegotiation[i] ~ dpois(lambda[i])
  }
}

nestling_inits <- function(){
  list(beta_0 = rnorm(1, 0, 1), 
       beta_ft = rnorm(1, 0, 1),
       beta_ps = rnorm(1, 0, 1),
       beta_inter = rnorm(1, 0, 1))
}

nestling_params <- c("beta_0", "beta_ft", "beta_ps", "beta_inter")

nestlings_fit_naive <- jags(data = nestlings_jags_data,
                            inits = nestling_inits,
                            parameters.to.save = nestling_params,
                            model.file = owls_naive,
                            n.chains = 5,
                            n.iter = 10000,
                            n.burnin = 2500,
                            n.thin = 1)

nestlings_fit_naive


## B. Model that adds an offset: ----
nestlings_jags_data <- list(SiblingNegotiation = Owls$SiblingNegotiation, 
                            FoodTreatment = as.numeric(as.factor(Owls$FoodTreatment)) - 1,
                            # note that we have to convert the FoodTreatment factor to
                            # a number to work appropriately in JAGS. Also, numeric coding needs to start at 0, not 1
                            SexParent = as.numeric(as.factor(Owls$SexParent)) - 1, 
                            # same with SexParent. Also, numeric coding needs to start at 0
                            N_owls = nrow(Owls),
                            ln_BroodSize = log(Owls$BroodSize))

owls_offset <- function(){
  # Priors:
  beta_0 ~ dnorm(0, 0.001)
  beta_ft ~ dnorm(0, 0.001)
  beta_ps ~ dnorm(0, 0.001)
  beta_inter ~ dnorm(0, 0.001)
  beta_offset ~ dnorm(1, 1000000)
  
  # Likelihood:
  for(i in 1:N_owls){
    # linear combination of predictors
    lin_comb[i] <- beta_0 + beta_ft * FoodTreatment[i] + beta_ps * SexParent[i] + 
      beta_inter * FoodTreatment[i] * SexParent[i] + beta_offset * ln_BroodSize[i]
    # invert link to get focal parameter
    lambda[i] <- exp(lin_comb[i])
    # specify likelihood prob model
    SiblingNegotiation[i] ~ dpois(lambda[i])
  }
}

nestling_inits <- function(){
  list(beta_0 = rnorm(1, 0, 1), 
       beta_ft = rnorm(1, 0, 1),
       beta_ps = rnorm(1, 0, 1),
       beta_inter = rnorm(1, 0, 1),
       beta_offset = rnorm(1, 1, .0000001))
}

nestling_params <- c("beta_0", "beta_ft", "beta_ps", "beta_inter", "beta_offset")

nestlings_fit_offset <- jags(data = nestlings_jags_data,
                             inits = nestling_inits,
                             parameters.to.save = nestling_params,
                             model.file = owls_offset,
                             n.chains = 5,
                             n.iter = 10000,
                             n.burnin = 2500,
                             n.thin = 1)

nestlings_fit_offset

## C. Model that adds a randome effect on nest ----
nestlings_jags_data <- list(SiblingNegotiation = Owls$SiblingNegotiation, 
                            FoodTreatment = as.numeric(as.factor(Owls$FoodTreatment)) - 1,
                            # note that we have to convert the FoodTreatment factor to
                            # a number to work appropriately in JAGS. Also, numeric coding needs to start at 0, not 1
                            SexParent = as.numeric(as.factor(Owls$SexParent)) - 1, 
                            # same with SexParent. Also, numeric coding needs to start at 0
                            N_owls = nrow(Owls),
                            ln_BroodSize = log(Owls$BroodSize),
                            Nest = as.numeric(as.factor(Owls$Nest)),
                            N_nest = length(levels(factor(Owls$Nest))))

owls_re <- function(){
  # Priors:
  beta_0 ~ dnorm(0, 0.001)
  beta_ft ~ dnorm(0, 0.001)
  beta_ps ~ dnorm(0, 0.001)
  beta_inter ~ dnorm(0, 0.001)
  beta_offset ~ dnorm(1, 1000000)
  # Random intercept priors:
  for(n in 1:N_nest){
    alpha[n] ~ dnorm(0, tau_nest)
  }
  sigma_sq_nest ~ dunif(0, 50)
  tau_nest <- pow(sigma_sq_nest, -2)
  
  # Likelihood:
  for(i in 1:N_owls){
    # linear combination of predictors
    lin_comb[i] <- beta_0 + beta_ft * FoodTreatment[i] + beta_ps * SexParent[i] + 
      beta_inter * FoodTreatment[i] * SexParent[i] + beta_offset * ln_BroodSize[i] + 
      alpha[Nest[i]]
    # invert link to get focal parameter
    lambda[i] <- exp(lin_comb[i])
    # specify likelihood prob model
    SiblingNegotiation[i] ~ dpois(lambda[i])
  }
}

nestling_inits <- function(){
  list(beta_0 = rnorm(1, 0, 1), 
       beta_ft = rnorm(1, 0, 1),
       beta_ps = rnorm(1, 0, 1),
       beta_inter = rnorm(1, 0, 1),
       beta_offset = rnorm(1, 1, .0000001),
       alpha = rnorm(n = 27, 0, 1),
       sigma_sq_nest = runif(1, 0, 50))
}

nestling_params <- c("beta_0", "beta_ft", "beta_ps", "beta_inter", "beta_offset", 
                     "alpha", "sigma_sq_nest")

nestlings_fit_re <- jags(data = nestlings_jags_data,
                         inits = nestling_inits,
                         parameters.to.save = nestling_params,
                         model.file = owls_re,
                         n.chains = 5,
                         n.iter = 10000,
                         n.burnin = 2500,
                         n.thin = 1)
nestlings_fit_re

## D. Model that update the probability model to negative binomial ----
nestlings_jags_data <- list(SiblingNegotiation = Owls$SiblingNegotiation, 
                            FoodTreatment = as.numeric(as.factor(Owls$FoodTreatment)) - 1,
                            # note that we have to convert the FoodTreatment factor to
                            # a number to work appropriately in JAGS. Also, numeric coding needs to start at 0, not 1
                            SexParent = as.numeric(as.factor(Owls$SexParent)) - 1, 
                            # same with SexParent. Also, numeric coding needs to start at 0
                            N_owls = nrow(Owls),
                            ln_BroodSize = log(Owls$BroodSize),
                            Nest = as.numeric(as.factor(Owls$Nest)),
                            N_nest = length(levels(factor(Owls$Nest))))

owls_negbin <- function(){
  # Priors:
  beta_0 ~ dnorm(0, 0.001)
  beta_ft ~ dnorm(0, 0.001)
  beta_ps ~ dnorm(0, 0.001)
  beta_inter ~ dnorm(0, 0.001)
  beta_offset ~ dnorm(1, 1000000)
  # Random intercept priors:
  for(n in 1:N_nest){
    alpha[n] ~ dnorm(0, tau_nest)
  }
  sigma_sq_nest ~ dunif(0, 50)
  tau_nest <- pow(sigma_sq_nest, -2)
  r ~ dunif(0, 50)
  
  # Likelihood:
  for(i in 1:N_owls){
    # linear combination of predictors
    lin_comb[i] <- beta_0 + beta_ft * FoodTreatment[i] + beta_ps * SexParent[i] + 
      beta_inter * FoodTreatment[i] * SexParent[i] + beta_offset * ln_BroodSize[i] + 
      alpha[Nest[i]]
    # invert link to get focal parameter
    lambda[i] <- exp(lin_comb[i])
    p[i] <- r / (r + lambda[i])
    # specify likelihood prob model
    SiblingNegotiation[i] ~ dnegbin(p[i], r)
  }
}

nestling_inits <- function(){
  list(beta_0 = rnorm(1, 0, 1), 
       beta_ft = rnorm(1, 0, 1),
       beta_ps = rnorm(1, 0, 1),
       beta_inter = rnorm(1, 0, 1),
       beta_offset = rnorm(1, 1, .0000001),
       alpha = rnorm(n = 27, 0, 1),
       sigma_sq_nest = runif(1, 0, 50),
       r = runif(1, 0, 50))
}

nestling_params <- c("beta_0", "beta_ft", "beta_ps", "beta_inter", "beta_offset", 
                     "alpha", "sigma_sq_nest", "r")

nestlings_fit_negbin <- jags(data = nestlings_jags_data,
                             inits = nestling_inits,
                             parameters.to.save = nestling_params,
                             model.file = owls_negbin,
                             n.chains = 5,
                             n.iter = 10000,
                             n.burnin = 2500,
                             n.thin = 1)

nestlings_fit_negbin

## E. Poisson model with zero inflation ----
nestlings_jags_data <- list(SiblingNegotiation = Owls$SiblingNegotiation, 
                            FoodTreatment = as.numeric(as.factor(Owls$FoodTreatment)) - 1,
                            # note that we have to convert the FoodTreatment factor to
                            # a number to work appropriately in JAGS. Also, numeric coding needs to start at 0, not 1
                            SexParent = as.numeric(as.factor(Owls$SexParent)) - 1, 
                            # same with SexParent. Also, numeric coding needs to start at 0
                            N_owls = nrow(Owls),
                            ln_BroodSize = log(Owls$BroodSize),
                            Nest = as.numeric(as.factor(Owls$Nest)),
                            N_nest = length(levels(factor(Owls$Nest))))

owls_zip <- function(){
  # Priors:
  beta_0_rate ~ dnorm(0, 0.001)
  beta_ft_rate ~ dnorm(0, 0.001)
  beta_ps_rate ~ dnorm(0, 0.001)
  beta_inter_rate ~ dnorm(0, 0.001)
  beta_offset_rate ~ dnorm(1, 1000000)
  
  beta_0_prob ~ dnorm(0, 0.001)
  beta_ft_prob ~ dnorm(0, 0.001)
  beta_ps_prob ~ dnorm(0, 0.001)
  beta_inter_prob ~ dnorm(0, 0.001)
  
  # Random intercept priors:
  for(n in 1:N_nest){
    alpha_rate[n] ~ dnorm(0, tau_nest_rate)
    alpha_prob[n] ~ dnorm(0, tau_nest_prob)
  }
  sigma_sq_nest_rate ~ dunif(0, 50)
  tau_nest_rate <- pow(sigma_sq_nest_rate, -2)
  
  sigma_sq_nest_prob ~ dunif(0, 50)
  tau_nest_prob <- pow(sigma_sq_nest_prob, -2)
  
  # Likelihood:
  for(i in 1:N_owls){
    # linear combination of predictors of the rate for nonstructural zeros
    lin_comb_rate[i] <- beta_0_rate + beta_ft_rate * FoodTreatment[i] + beta_ps_rate * SexParent[i] + 
      beta_inter_rate * FoodTreatment[i] * SexParent[i] + beta_offset_rate * ln_BroodSize[i] + 
      alpha_rate[Nest[i]]
    # linear combination of predictors of the probability for structural zeros
    lin_comb_prob[i] <- beta_0_prob + beta_ft_prob * FoodTreatment[i] + beta_ps_prob * SexParent[i] + 
      beta_inter_prob * FoodTreatment[i] * SexParent[i] + 
      alpha_prob[Nest[i]]
    
    # invert logit to get pi
    pi[i] <- exp(lin_comb_prob[i]) / (1 + exp(lin_comb_prob[i]))
    zero[i] ~ dbern(pi[i])
    
    # invert link to get lambda
    lambda[i] <- exp(lin_comb_rate[i])
    lambda_hacked[i] <- lambda[i]*(1-zero[i]) + 1e-10*zero[i]
    SiblingNegotiation[i] ~ dpois(lambda_hacked[i])
  }
}

nestling_inits <- function(){
  list(beta_0_rate = rnorm(1, 0, 1), 
       beta_ft_rate = rnorm(1, 0, 1),
       beta_ps_rate = rnorm(1, 0, 1),
       beta_inter_rate = rnorm(1, 0, 1),
       beta_offset_rate = rnorm(1, 1, .0000001),
       alpha_rate = rnorm(n = 27, 0, 1),
       sigma_sq_nest_rate = runif(1, 0, 50),
       beta_0_prob = rnorm(1, 0, 1), 
       beta_ft_prob = rnorm(1, 0, 1),
       beta_ps_prob = rnorm(1, 0, 1),
       beta_inter_prob = rnorm(1, 0, 1),
       alpha_prob = rnorm(n = 27, 0, 1),
       sigma_sq_nest_prob = runif(1, 0, 50))
}

nestling_params <- c("beta_0_prob", "beta_ft_prob", "beta_ps_prob", "beta_inter_prob", "beta_offset_rate", 
                     "alpha_prob", "sigma_sq_nest_prob", 
                     "beta_0_rate", "beta_ft_rate", "beta_ps_rate", "beta_inter_rate",
                     "alpha_rate", "sigma_sq_nest_rate")

nestlings_fit_zip <- jags(data = nestlings_jags_data,
                          inits = nestling_inits,
                          parameters.to.save = nestling_params,
                          model.file = owls_zip,
                          n.chains = 5,
                          n.iter = 5000,
                          n.burnin = 2500,
                          n.thin = 1)
nestlings_fit_zip





