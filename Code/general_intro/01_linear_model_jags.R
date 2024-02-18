# installation instructions: https://onlinelibrary.wiley.com/doi/pdf/10.1002/9781119287995.app1
# a useful starter tutorial: http://biometry.github.io/APES//LectureNotes/StatsCafe/Linear_models_jags.html
library(R2jags)

# 00. Big bang ----
library(Sleuth3)
bigbang <- case0701

# plotting the data
plot(bigbang$Distance ~ bigbang$Velocity,
     las = 1, pch = 16,
     xlab = "Velocity",
     ylab = "Distance")

bigbang_jags_data <- list(Distance = bigbang$Distance, 
                          Velocity = bigbang$Velocity,
                          N = nrow(bigbang))

bigbang_jags_model <- function(){
  # Priors:
  beta_0 ~ dnorm(0, 0.01) # intercept
  beta_1 ~ dnorm(0, 0.01) # slope
  sigma ~ dunif(0, 100) # standard deviation
  tau <- 1 / (sigma * sigma) # sigma^2 doesn't work in JAGS
  
  # Likelihood:
  for (i in 1:N){
    lin_comb[i] <- beta_0 + beta_1 * Velocity[i]
    mu[i] <- lin_comb[i]
    Distance[i] ~ dnorm(mu[i], tau) # tau is precision (1 / variance)
  }
}

bigbang_inits <- function(){
  list(beta_0 = rnorm(1), 
       beta_1 = rnorm(1), 
       sigma = runif(1))
}

bigbang_params <- c("beta_0", "beta_1", "sigma")

bigbang_fit <- jags(data = bigbang_jags_data, 
                    inits = bigbang_inits, 
                    parameters.to.save = bigbang_params, 
                    model.file = bigbang_jags_model,
                    n.chains = 3,
                    n.iter = 12000, 
                    n.burnin = 2000, 
                    n.thin = 10, 
                    DIC = F)
bigbang_fit

traceplot(bigbang_fit, mfrow = c(2, 2), ask = F)

bigbang_mcmc <- as.mcmc(bigbang_fit)
plot(bigbang_mcmc)

