
library(jagsUI)

#setup~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#set your working directory to the project root before running this script
#the folder should contain: data/ (with jags.data1.rds) and code/ (for model output)
#setwd("path/to/project/root")

#load jags data
jags.data1 <- readRDS("data/jags.data1.rds")

#write jags model file to the working directory
cat(file="model1.txt", "
model{

#likelihood~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#loop through the clutches (i.e., nClutches) each day (i.e., day)

    for(i in 1:nClutches){
      for(t in init[i]:last[i]){

        #lambda according to a normal distribution
        lambda[i,t] ~ dnorm(mu[i,t], tau)

        #linear predictor
        mu[i,t] <- alpha

        #predict N (i.e.,clutch sizes)
        N[i,t+1] <- N[i,t] + lambda[i,t]

        #observation model (i.e., clutch size corrected due to observation error)
        y[i,t] ~ dnorm(N[i,t], tau.y)
      }
    }

#priors~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#lambda intercept, standard deviation, and variance

    alpha ~ dnorm(0, 0.001)
    sigma ~ dunif(0,10)
    tau <- 1 / (sigma * sigma)

#inital clutch size

    for(i in 1:nClutches){
      N[i,init[i]] ~ dunif(-50,50)
    }

#standard deviation and variance for the observation model

    sigma.y ~ dunif(0,10)
    tau.y <- 1 / (sigma.y * sigma.y)

}
")

#initial values
inits <- function(){list()}

#parameters monitored
parameters <- c("alpha", "sigma", "sigma.y", "lambda", "N")

#mcmc settings
ni <- 80000; nb <- 30000; nc <- 4; nt <- 10

#call jags from r
m <- jags(jags.data1, inits, parameters, "model1.txt", n.chains = nc,
          n.thin = nt, n.iter = ni, n.burnin = nb, parallel = F)
summary(m)



